#!/usr/bin/env python3
"""RNA-reference Bowtie2 mapping with verified aligner-unmapped FASTQ handoff.

Uses the established very-sensitive, end-to-end RNA mapping and MAPQ/flag policy.
This I/O adapter adds read accounting without imposing genome CIGAR/decoy rules.
"""
import argparse
from collections import Counter
from contextlib import closing
import json
from pathlib import Path
import sqlite3
import subprocess
import tempfile

from bam_count import command, policy
from core import sha256
from genome_common import (check_pins, dictionary, executable, index_files,
                           ingest_fastq, new_output, pin, verify_handoff)


def classify(flag, mapq, minimum, flags):
    if flag & (256 | 2048):
        return 'nonprimary'
    if flag & 4:
        return 'unmapped'
    if flag & flags:
        return 'flags'
    return 'low_mapq' if mapq < minimum else 'accepted'


def process_rna_sam(source, db, minimum, flags):
    outcomes, records = Counter(), Counter()
    with open(source) as handle:
        for line in handle:
            if line.startswith('@'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 11:
                raise ValueError('Malformed RNA SAM')
            name, flag, mapq = fields[0], int(fields[1]), int(fields[4])
            if flag & 1 or not 0 <= mapq <= 255:
                raise ValueError('Require single-end SAM with valid MAPQ')
            records['all'] += 1
            reason = classify(flag, mapq, minimum, flags)
            if reason == 'nonprimary':
                records['nonprimary'] += 1
                continue
            row = db.execute('select seen from reads where id=?', (name,)).fetchone()
            if row is None or row[0]:
                raise ValueError('Missing/duplicate RNA primary read: ' + name)
            db.execute('update reads set seen=1,unmapped=? where id=?',
                       (int(reason == 'unmapped'), name))
            outcomes[reason] += 1
    if db.execute('select count(*) from reads where seen=0').fetchone()[0]:
        raise ValueError('RNA aligner omitted reads')
    return {'primary_outcomes': dict(outcomes), 'alignment_records': dict(records),
            'mapq_minimum': minimum, 'excluded_flags': flags}


def run(a):
    flags = policy(a.mapq, a.exclude_duplicates)
    out = new_output(a)
    before = pin([a.fastq, a.reference, *index_files(a.index, 'bowtie2')])
    meta = {'status': 'running', 'parameters': vars(a), 'inputs_sha256': before,
            'code_sha256': {p.name: sha256(p) for p in
                           [Path(__file__), Path(__file__).with_name('genome_common.py')]}}
    try:
        with (out/'commands.jsonl').open('w') as log, tempfile.TemporaryDirectory(prefix='rna-', dir=out) as tmp:
            tmp = Path(tmp)
            with closing(sqlite3.connect(tmp/'reads.sqlite')) as db:
                count = ingest_fastq(a.fastq, db)
                dictionary(a.index, 'bowtie2', out, log, a.reference)
                meta['versions'] = {t: executable(t) for t in ('bowtie2', 'samtools')}
                # Relative read/index aliases also support whitespace in paths.
                for part in index_files(a.index, 'bowtie2'):
                    suffix = str(part)[len(str(Path(a.index).resolve())):]
                    (tmp/('rna'+suffix)).symlink_to(part)
                reads = 'reads.fq.gz' if a.fastq.endswith('.gz') else 'reads.fq'
                (tmp/reads).symlink_to(Path(a.fastq).resolve())
                args = ['bowtie2', '-p', str(a.threads), '--very-sensitive',
                        '--seed', str(a.seed), '-x', 'rna', '-U', reads,
                        '--un', 'unmapped.fastq', '-S', 'all.sam']
                if a.norc:
                    args.append('--norc')
                command(args, log, cwd=tmp)
                qc = process_rna_sam(tmp/'all.sam', db, a.mapq, flags)
                qc['forwarded_unmapped_reads'] = verify_handoff(tmp/'unmapped.fastq', db)
                if sum(qc['primary_outcomes'].values()) != count:
                    raise ValueError('RNA read accounting mismatch')
                (tmp/'unmapped.fastq').replace(out/'unmapped.fastq')
                filtered = tmp/'filtered.bam'
                command(['samtools', 'view', '-b', '-F', str(flags), '-q', str(a.mapq),
                         '-o', str(filtered), str(tmp/'all.sam')], log)
                for source, name in [(tmp/'all.sam', 'alignments'), (filtered, 'accepted')]:
                    command(['samtools', 'sort', '-@', str(a.threads-1), '-m', '256M',
                             '-T', str(tmp/name), '-o', str(out/(name+'.bam')), str(source)], log)
                    command(['samtools', 'quickcheck', '-v', str(out/(name+'.bam'))], log)
                (out/'qc.json').write_text(json.dumps(qc, indent=2)+'\n')
                meta.update(input_reads=count, qc=qc)
                db.commit()
        check_pins(before)
        meta.update(status='complete', inputs_unchanged=True)
        meta['outputs_sha256'] = {p.name: sha256(p) for p in out.iterdir() if p.is_file()}
    except Exception as exc:
        meta.update(status='failed', error=str(exc))
        raise
    finally:
        (out/'provenance.json').write_text(json.dumps(meta, indent=2)+'\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.set_defaults(decoys=None, no_decoys=False, reference_sha256=None)
    p.add_argument('--fastq', required=True)
    p.add_argument('--index', required=True)
    p.add_argument('--reference', required=True)
    p.add_argument('--output-dir', required=True, help='New external output directory')
    p.add_argument('--threads', type=int, default=1, help='Sequential aligner/sort CPU budget')
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--mapq', type=int, default=20)
    p.add_argument('--exclude-duplicates', action='store_true')
    p.add_argument('--norc', action='store_true')
    a = p.parse_args()
    try:
        run(a)
    except (ValueError, OSError, subprocess.SubprocessError) as exc:
        p.exit(1, str(exc)+'\n')


if __name__ == '__main__':
    main()
