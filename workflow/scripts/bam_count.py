"""Shared read-only BAM preparation for the standalone and Snakemake counters."""
import json
from pathlib import Path
import re
import subprocess
import tempfile
import time

from core import matrix, records, sha256, threshold, validate_fasta

RESERVED = {'chromosome', 'strand', 'start', 'end'}


def sample_names(names, bams):
    if not names or len(names) != len(bams) or len(set(names)) != len(names):
        raise ValueError('Provide one unique sample name per BAM, in BAM order')
    if any(not isinstance(n, str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', n) or n in RESERVED for n in names):
        raise ValueError('Invalid or reserved sample name')
    if len({str(Path(b).resolve()) for b in bams}) != len(bams):
        raise ValueError('Duplicate BAM input paths')


def policy(mapq, exclude_duplicates, strand='both', minimum=0):
    if type(mapq) is not int or not 0 <= mapq <= 255:
        raise ValueError('MAPQ must be an integer from 0 to 255')
    if type(exclude_duplicates) is not bool:
        raise ValueError('exclude_duplicates must be boolean')
    if strand not in ('plus', 'minus', 'both'):
        raise ValueError('strand must be plus, minus or both')
    if type(minimum) is not int or minimum < 0:
        raise ValueError('minimum total count must be a nonnegative integer')
    return 2308 | (1024 if exclude_duplicates else 0)


def command(args, log, **kwargs):
    start = time.monotonic()
    log.write(json.dumps({'command': [str(a) for a in args]}) + '\n'); log.flush()
    result = subprocess.run(args, stderr=log, **kwargs)
    log.write(json.dumps({'exit_status': result.returncode, 'wall_seconds': time.monotonic()-start}) + '\n')
    log.flush()
    result.check_returncode()
    return result


def reference_sizes(fasta, output, expected_sha256=None):
    if Path(fasta).resolve() == Path(output).resolve():
        raise ValueError('Cannot overwrite reference')
    digest = sha256(fasta)
    from reference_hashes import expected_digest
    expected_sha256 = expected_digest(expected_sha256, 'Reference SHA-256')
    if expected_sha256 is not None and digest != expected_sha256:
        raise ValueError('Reference SHA-256 mismatch')
    validate_fasta(fasta)
    sizes = {h.split()[0]: len(s) for h, s in records(fasta)}
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    Path(output).write_text(''.join(f'{name}\t{length}\n' for name, length in sizes.items()))
    return sizes


def read_sizes(path):
    return {row.split('\t')[0]: int(row.split('\t')[1]) for row in Path(path).read_text().splitlines()}


def coordinate_order(bam, dictionary, log):
    """Validate actual alignment order without sorting, indexing or CIGAR parsing."""
    ranks = {name: i for i, name in enumerate(dictionary)}
    args = ['samtools', 'view', str(bam)]
    log.write(json.dumps({'command': args, 'purpose': 'validate actual coordinate order'})+'\n')
    log.flush()
    start, previous, ordered, count = time.monotonic(), None, True, 0
    process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=log, text=True)
    try:
        for line in process.stdout:
            fields = line.split('\t', 4)
            key = (ranks[fields[2]], int(fields[3]))
            if previous is not None and key < previous:
                ordered = False
            previous = key
            count += 1
        status = process.wait()
    except Exception:
        process.terminate()
        process.wait()
        raise
    finally:
        process.stdout.close()
    log.write(json.dumps({'exit_status': status, 'wall_seconds': time.monotonic()-start,
                          'coordinate_order_valid': ordered, 'checked_alignments': count})+'\n')
    log.flush()
    if status:
        raise subprocess.CalledProcessError(status, args)
    return ordered


def prepare_bam(source, output, sizes, mapq=20, exclude_duplicates=False, log=None):
    """Re-filter an output-local copy; sort only if actual coordinates are unordered."""
    flags = policy(mapq, exclude_duplicates)
    source, output = Path(source).resolve(), Path(output).resolve()
    if source == output:
        raise ValueError('Cannot overwrite input BAM')
    if log is None:
        raise ValueError('A command log is required')
    before = source.stat()
    command(['samtools', 'quickcheck', '-v', str(source)], log)
    header = command(['samtools', 'view', '--no-PG', '-H', str(source)], log,
                     stdout=subprocess.PIPE, text=True).stdout
    dictionary = {}
    sort_order = 'unknown'
    for line in header.splitlines():
        tags = dict(field.split(':', 1) for field in line.split('\t')[1:] if ':' in field)
        if line.startswith('@SQ\t'):
            if tags['SN'] in dictionary:
                raise ValueError('Duplicate BAM sequence identifier')
            dictionary[tags['SN']] = int(tags['LN'])
        elif line.startswith('@HD\t'):
            sort_order = tags.get('SO', 'unknown')
    if dictionary != read_sizes(sizes):
        raise ValueError('BAM/reference dictionary mismatch')
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='prepare-', dir=output.parent) as tmp:
        filtered = str(Path(tmp) / 'filtered.bam')
        command(['samtools', 'view', '-b', '-F', str(flags), '-q', str(mapq),
                 '-o', filtered, str(source)], log)
        ordered = coordinate_order(filtered, dictionary, log)
        if ordered:
            Path(filtered).replace(output)
        else:
            command(['samtools', 'sort', '-@', '0', '-m', '256M', '-T', str(Path(tmp)/'sort'),
                     '-o', str(output), filtered], log)
        command(['samtools', 'quickcheck', '-v', str(output)], log)
    after = source.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError('Input BAM changed during preparation')
    return {'source': str(source), 'size_bytes': before.st_size, 'mtime_ns': before.st_mtime_ns,
            'declared_sort_order': sort_order, 'actual_filtered_order_valid': ordered,
            'sort_performed': not ordered, 'sort_policy': 'only unordered output-local copies',
            'mapq': mapq, 'exclude_flags': flags, 'source_unchanged_stat': True}


def count_prepared(bams, names, sizes, output, strand='both', minimum=0, log=None):
    policy(20, False, strand, minimum)
    sample_names(names, bams)
    output = Path(output)
    with tempfile.TemporaryDirectory(prefix='matrix-', dir=output.parent) as tmp:
        raw = Path(tmp) / 'raw.tsv'
        summary = matrix(bams, names, sizes, raw, strand, log)
        threshold(raw, output, minimum)
    return dict(summary, min_total_count=minimum)
