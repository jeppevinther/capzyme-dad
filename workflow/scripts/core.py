"""Portable FASTA validation, one-base matrix assembly and thresholding.

Counts are produced by bedtools genomecov -5 -dz, not by a replacement CIGAR
interpretation. The caller supplies accepted, coordinate-sorted BAMs.
"""
import csv
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import tempfile


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def records(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as fh:
        header, sequence = None, []
        for line in fh:
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(sequence).upper()
                header, sequence = line[1:].strip(), []
            elif line.strip():
                if header is None:
                    raise ValueError('Sequence before FASTA header')
                sequence.append(''.join(line.split()))
        if header is not None:
            yield header, ''.join(sequence).upper()


def validate_fasta(path):
    seen = set()
    for header, sequence in records(path):
        name = header.split()[0] if header else ''
        if not name or name in seen or not sequence or not re.fullmatch('[ACGTRYSWKMBDHVN]+', sequence):
            raise ValueError(f'Invalid/duplicate FASTA record: {name!r}')
        seen.add(name)
    if not seen:
        raise ValueError('FASTA contains no records')


def collapse(source, prefix, k, drops, legacy):
    # Use the original implementation directly, with integrity verified first.
    expected = 'cce88700a8b9d1108d4c5aa3758d7c50feeb1958fd35a1f8287edb7a2ac10336'
    if sha256(legacy) != expected:
        raise ValueError('Legacy collapse source hash changed; review compatibility before running')
    spec = importlib.util.spec_from_file_location('capzyme_legacy_prefix', legacy)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    merged, rows, pre, dropped, skipped, stats = module.merge_by_prefix(source, k=k, drop_biotypes=set(drops))
    if not merged:
        raise ValueError('No ncRNA records remain after filtering')
    module.write_fasta(merged, str(prefix) + '.fa')
    module.write_tsv(rows, str(prefix) + '.tsv')
    return dict(k=k, source_sha256=sha256(source), legacy_sha256=expected,
                records=stats, pre_biotypes=pre, dropped=dropped, skipped=skipped, groups=len(rows))


def combine(inputs, output):
    seen = set()
    with open(output, 'w') as out:
        for path in inputs:
            validate_fasta(path)
            for header, sequence in records(path):
                name = header.split()[0]
                if name in seen:
                    raise ValueError(f'Duplicate combined reference ID: {name}')
                seen.add(name)
                out.write(f'>{header}\n')
                for start in range(0, len(sequence), 60):
                    out.write(sequence[start:start + 60] + '\n')


def natural_key(text):
    return tuple((1, int(part)) if part.isdigit() else (0, part) for part in re.split(r'(\d+)', text))


def matrix(bams, names, fai, output, strand='both', log=None):
    if not names or len(names) != len(bams) or len(set(names)) != len(names):
        raise ValueError('Unique sample names must match BAM order')
    sizes = {}
    with open(fai) as fh:
        for line in fh:
            chrom, size, *_ = line.rstrip().split('\t')
            sizes[chrom] = int(size)
    sites = {}
    if log:
        version = subprocess.check_output(['bedtools', '--version'], text=True).strip()
        log.write(json.dumps({'bedtools_version': version}) + '\n'); log.flush()
    strands = ['+', '-'] if strand == 'both' else ['+' if strand == 'plus' else '-']
    # Disk-backed subprocess output keeps pipe handling simple and avoids
    # collecting a second full copy of each genomecov output in memory.
    for column, bam in enumerate(bams):
        for sign in strands:
            cmd = ['bedtools', 'genomecov', '-ibam', str(bam), '-strand', sign, '-5', '-dz']
            if log:
                log.write(json.dumps({'command': cmd}) + '\n'); log.flush()
            with tempfile.TemporaryFile(mode='w+t', dir=Path(output).parent) as data:
                result = subprocess.run(cmd, stdout=data, stderr=log, check=False)
                if log:
                    log.write(f'exit_status={result.returncode}\n'); log.flush()
                result.check_returncode()
                data.seek(0)
                for line in data:
                    chrom, pos, count = line.rstrip().split('\t')
                    pos, count = int(pos), int(count)
                    if chrom not in sizes or not 0 <= pos < sizes[chrom] or count < 1:
                        raise ValueError(f'Invalid genomecov coordinate/count: {line.rstrip()}')
                    key = chrom, sign, pos + 1
                    values = sites.setdefault(key, [0] * len(names))
                    values[column] += count
    with open(output, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t', lineterminator='\n')
        writer.writerow(['chromosome', 'strand', 'start', 'end', *names])
        for chrom, sign, pos in sorted(sites, key=lambda x: (natural_key(x[0]), x[2], x[1])):
            writer.writerow([chrom, sign, pos, pos, *sites[chrom, sign, pos]])
    return dict(samples=names, rows=len(sites), strand=strand,
                accepted_ends_per_sample=[sum(v[i] for v in sites.values()) for i in range(len(names))])


def threshold(source, output, minimum):
    if type(minimum) is not int or minimum < 0:
        raise ValueError('minimum must be a nonnegative integer')
    with open(source, newline='') as inp, open(output, 'w', newline='') as out:
        reader = csv.reader(inp, delimiter='\t')
        writer = csv.writer(out, delimiter='\t', lineterminator='\n')
        writer.writerow(next(reader))
        for row in reader:
            if sum(map(int, row[4:])) >= minimum:
                writer.writerow(row)
