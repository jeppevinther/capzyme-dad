#!/usr/bin/env python3
"""Independently decode BigWigs and compare every interval to an integer matrix.

Uses disk-backed SQLite and 1-Mb BigWig windows, not the writer/counting engine.
The JSON specification identifies source matrices/dictionaries, ordered samples,
audited depths and raw/RPM/DESeq2 output arms. DESeq2 verification requires a
separate direct-R factor calculation; this checker never estimates size factors.
"""
import argparse
import csv
import itertools
import json
import math
from pathlib import Path
import shutil
import sqlite3
import struct
import tempfile

from core import sha256


def sizes(path):
    result = {}
    for line in Path(path).read_text().splitlines():
        fields = line.split('\t')
        if len(fields) != 2 or fields[0] in result or not fields[0] or not fields[1].isdigit() or int(fields[1]) < 1:
            raise ValueError('Invalid dictionary')
        result[fields[0]] = int(fields[1])
    if not result: raise ValueError('Empty dictionary')
    return result


def intervals(bw, dictionary):
    for chrom, length in dictionary.items():
        for start in range(0, length, 1000000):
            for left, right, value in bw.intervals(chrom, start, min(start+1000000, length)) or ():
                # Reject coverage/extended intervals instead of silently expanding them.
                if right != left+1 or not start <= left < min(start+1000000, length):
                    raise ValueError('Nonsingle-base or overlapping-window interval')
                yield chrom, left, right, value


def _validate_space(spec, workspace):
    import pyBigWig
    matrix, dictionary_path = Path(spec['matrix']), Path(spec['sizes'])
    before = {str(p): sha256(p) for p in (matrix, dictionary_path)}
    if before[str(matrix)] != spec['matrix_sha256']: raise ValueError('Matrix trusted checksum mismatch')
    dictionary = sizes(dictionary_path)
    names = spec['samples']
    if len(names) != len(set(names)) or not names: raise ValueError('Sample order must be explicit/unique')
    dbpath = Path(workspace)/'matrix.sqlite'
    if dbpath.exists(): raise ValueError('New verification workspace required')
    report = {'samples': names, 'source_sha256': before, 'arms': {}, 'features': 0, 'sample_totals': [0]*len(names)}
    with sqlite3.connect(dbpath) as db:
        db.execute('pragma cache_size=-65536')
        db.execute('pragma temp_store=FILE')
        columns = ','.join(f'c{i} integer not null' for i in range(len(names)))
        db.execute(f'create table counts (chrom text, strand text, pos integer, {columns}, primary key(chrom,strand,pos)) without rowid')
        with matrix.open() as f:
            reader = csv.reader(f, delimiter='\t')
            if next(reader) != ['chromosome', 'strand', 'start', 'end', *names]: raise ValueError('Matrix sample order/header mismatch')
            batch = []
            for row in reader:
                if len(row) != 4+len(names) or not all(x.isdigit() for x in row[2:]): raise ValueError('Invalid integer matrix row')
                pos, end = int(row[2]), int(row[3])
                if row[0] not in dictionary or row[1] not in ('+', '-') or pos != end or not 1 <= pos <= dictionary[row[0]]: raise ValueError('Invalid matrix coordinate')
                counts = list(map(int, row[4:]))
                report['sample_totals'] = [a+b for a, b in zip(report['sample_totals'], counts)]
                report['features'] += 1
                batch.append((row[0], row[1], pos-1, *counts))
                if len(batch) == 10000:
                    db.executemany('insert into counts values ('+','.join('?' for _ in batch[0])+')', batch); batch.clear()
            if batch: db.executemany('insert into counts values ('+','.join('?' for _ in batch[0])+')', batch)
        db.commit()
        # One read transaction avoids remote lock acquisition on every contig.
        db.execute('begin')
        if report['sample_totals'] != spec['accepted_depths']: raise ValueError('Matrix totals differ from audited accepted depths')
        for arm in spec['arms']:
            if arm['method'] not in ('raw', 'rpm', 'deseq2'): raise ValueError('Unsupported validation method')
            root = Path(arm['directory'])
            with (root/'normalization_factors.tsv').open() as f: factors = list(csv.DictReader(f, delimiter='\t'))
            if [r['sample_id'] for r in factors] != names: raise ValueError('Factor cohort order mismatch')
            meta = json.loads((root/'provenance.json').read_text())
            if meta['status'] != 'complete' or meta['cohort'] != names or meta['emitted_samples'] != names or meta['normalization'] != arm['method']:
                raise ValueError('Incomplete/wrong cohort provenance')
            independent = None
            if arm['method']=='deseq2':
                source = Path(arm['independent_factors'])
                if sha256(source)!=arm['independent_factors_sha256']: raise ValueError('Independent factors checksum mismatch')
                with source.open() as f: reference_factors=list(csv.DictReader(f, delimiter='\t'))
                if [r['sample_id'] for r in reference_factors]!=names: raise ValueError('Independent factor sample order mismatch')
                independent=[float(r['size_factor']) for r in reference_factors]
                if any(not math.isfinite(s) or s<=0 for s in independent): raise ValueError('Invalid independent size factors')
                if arm['estimator'] not in ('ratio','poscounts') or meta.get('estimator')!=arm['estimator']:
                    raise ValueError('DESeq2 estimator mismatch')
                agreement=[]
                for row, sf in zip(factors, independent):
                    actual=float(row['size_factor'])
                    if not math.isfinite(actual) or actual<=0 or not math.isclose(actual,sf,rel_tol=1e-12,abs_tol=0.):
                        raise ValueError('Independent DESeq2 size factor mismatch: '+row['sample_id'])
                    agreement.append(abs(actual-sf)/sf)
            expected_files, tracks = set(), []
            for i, name in enumerate(names):
                scale = (1./float(factors[i]['size_factor']) if independent is not None else
                         1. if arm['method'] == 'raw' else 1000000/spec['accepted_depths'][i])
                if factors[i]['method'] != arm['method'] or float(factors[i]['multiplier']) != scale: raise ValueError('Factor calculation mismatch')
                if arm['method'] == 'rpm' and int(factors[i]['accepted_records']) != spec['accepted_depths'][i]: raise ValueError('RPM depth mismatch')
                for strand in arm['strands']:
                    filename = name + ('.plus.bw' if strand == '+' else '.minus.bw')
                    expected_files.add(filename)
                    path = root/'tracks'/filename
                    def expected():
                        for chrom in dictionary:
                            for pos, count in db.execute(f'select pos,c{i} from counts where chrom=? and strand=? and c{i}>0 order by pos', (chrom, strand)):
                                exact = count/independent[i] if independent is not None else count*scale
                                rounded = struct.unpack('f', struct.pack('f', exact))[0]
                                yield chrom, pos, pos+1, rounded, exact
                    values = 0; max_error = 0.; max_relative = 0.
                    # A single sequential copy avoids many tiny remote reads.
                    # Cache only one track at a time, never BAMs or raw reads.
                    cached = Path(workspace)/'current.bw'
                    shutil.copyfile(path, cached)
                    cached_hash = sha256(cached)
                    with pyBigWig.open(str(cached)) as bw:
                        if list(bw.chroms().items()) != list(dictionary.items()): raise ValueError('BigWig dictionary/order mismatch')
                        iterator = iter(expected())
                        first = next(iterator, None)
                        empty = first is None
                        expected_rows = [(next(iter(dictionary)), 0, 1, 0., 0.)] if empty else itertools.chain([first], iterator)
                        for a, b in itertools.zip_longest(expected_rows, intervals(bw, dictionary)):
                            if a is None or b is None or a[:4] != b:
                                raise ValueError(f'BigWig interval/value mismatch: {filename}, expected={a}, observed={b}')
                            if not empty:
                                values += 1
                                max_error = max(max_error, abs(b[3]-a[4]))
                                max_relative = max(max_relative, abs((b[3]-a[4])/a[4]))
                    if empty != (filename in meta['empty_tracks']): raise ValueError('Empty-track declaration mismatch')
                    h = sha256(path)
                    if h != cached_hash: raise ValueError('Track changed during cached verification')
                    if meta['outputs_sha256']['tracks/'+filename] != h: raise ValueError('Published track checksum mismatch')
                    tracks.append({'file': filename, 'bytes': path.stat().st_size, 'sha256': h, 'nonzero_values_checked': values,
                                   'empty_zero_sentinel': empty, 'max_float32_absolute_error': max_error, 'max_float32_relative_error': max_relative})
            if {p.name for p in (root/'tracks').iterdir()} != expected_files: raise ValueError('Missing/extra sample or strand tracks')
            report['arms'][arm['method']] = {'tracks': tracks, 'track_count': len(tracks), 'depth_scope': spec['depth_scope']}
            if independent is not None:
                if sha256(source)!=arm['independent_factors_sha256']: raise ValueError('Independent factors changed')
                report['arms'][arm['method']].update(estimator=arm['estimator'],
                    max_factor_relative_difference=max(agreement), factor_relative_tolerance=1e-12,
                    independent_factors_sha256=arm['independent_factors_sha256'],
                    units='DESeq2 size-factor-normalized 5-prime counts (not RPM)')
    if before != {str(p): sha256(p) for p in (matrix, dictionary_path)}: raise ValueError('Source changed during verification')
    report['status'] = 'passed'
    return report


def validate_space(spec, workspace):
    # SQLite scratch is node-local, never a shared NFS database. Only compact
    # verification results are retained on project storage; no source is copied.
    # The production track writer and normalization are unaffected.
    Path(workspace).mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='capzyme-bw-check-', dir='/tmp') as local:
        result = _validate_space(spec, local)
        result['verification_scratch'] = 'node-local temporary SQLite and one-track cache; removed after verification'
        return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--spec', required=True); p.add_argument('--output-dir', required=True)
    a = p.parse_args(); out = Path(a.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=False)
    report = {'status': 'running', 'spaces': {}}
    try:
        for name, spec in json.loads(Path(a.spec).read_text()).items():
            work = out/name; work.mkdir()
            report['spaces'][name] = validate_space(spec, work)
        report['status'] = 'passed'
    except Exception as exc:
        report.update(status='failed', error=str(exc)); raise
    finally:
        (out/'validation.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__': main()
