#!/usr/bin/env python3
"""Count single-end aligned 5′ ends using the shared BEDTools matrix engine."""
import argparse
import json
from pathlib import Path
import subprocess
import shutil
import sys
import tempfile

from bam_count import count_prepared, policy, prepare_bam, reference_sizes, sample_names
from core import sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__, epilog=(
        'Requires samtools and bedtools on PATH. Sources are read-only; no source-side indexes. '
        'All BAMs are re-filtered into temporary output-local files; only unordered copies are sorted. '
        'Previously removed alignments cannot be recovered. Default flags: -F 2308; duplicates and QC-fail retained. '
        'Output TSV: chromosome, strand, start, end, ordered samples; start=end is 1-based.'))
    parser.add_argument('bams', nargs='+', help='One or more single-end BAMs, in desired sample order')
    parser.add_argument('-g', '--reference', required=True, help='Reference FASTA matching BAM dictionaries (not modified)')
    parser.add_argument('--reference-sha256', help='Optional required reference SHA-256')
    parser.add_argument('--sample-name', action='append', dest='names', help='Repeat once per BAM; defaults to each BAM filename without .bam')
    parser.add_argument('-o', '--output', required=True, help='Output matrix TSV')
    parser.add_argument('--mapq', type=int, default=20, help='Minimum MAPQ, inclusive (default: 20)')
    parser.add_argument('--exclude-duplicates', action='store_true', help='Also exclude flag 1024 (default: retain duplicates)')
    parser.add_argument('--strand', choices=['plus', 'minus', 'both'], default='both')
    parser.add_argument('--min-total-count', type=int, default=0, help='Inclusive row total across all samples (default: 0)')
    parser.add_argument('--log', help='Command/provenance log (default: OUTPUT.jsonl)')
    args = parser.parse_args()
    try:
        policy(args.mapq, args.exclude_duplicates, args.strand, args.min_total_count)
        names = args.names or [Path(b).stem for b in args.bams]
        sample_names(names, args.bams)
        output = Path(args.output).resolve()
        log_path = Path(args.log or str(output)+'.jsonl').resolve()
        inputs = {Path(b).resolve() for b in args.bams} | {Path(args.reference).resolve()}
        if output == log_path or output in inputs or log_path in inputs:
            raise ValueError('Output/log paths must differ from all inputs and each other')
        # Refuse accidental replacement of a previous result; Snakemake manages its own outputs separately.
        if output.exists() or log_path.exists():
            raise ValueError('Output/log already exists; choose new paths')
        output.parent.mkdir(parents=True, exist_ok=True)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open('x') as log, tempfile.TemporaryDirectory(prefix='count5-', dir=output.parent) as tmp:
            versions = {tool: subprocess.check_output([tool, '--version'], text=True).splitlines()[0]
                        for tool in ('samtools', 'bedtools')}
            sizes = Path(tmp)/'reference.sizes.tsv'
            reference_sizes(args.reference, sizes, args.reference_sha256)
            prepared, inputs_report = [], []
            for i, source in enumerate(args.bams):
                bam = Path(tmp)/f'{i}.accepted.bam'
                inputs_report.append(prepare_bam(source, bam, sizes, args.mapq, args.exclude_duplicates, log))
                prepared.append(bam)
            summary = count_prepared(prepared, names, sizes, output, args.strand, args.min_total_count, log)
            log.write(json.dumps({'exit_status': 0, 'versions': versions, 'inputs': inputs_report,
                                 'executables': {tool: shutil.which(tool) for tool in ('samtools', 'bedtools')},
                                 'python': sys.version,
                                 'reference_sha256': sha256(args.reference), 'matrix': summary,
                                 'output_sha256': sha256(output)})+'\n')
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        parser.exit(1, f'count_5ends: {error}\n')


if __name__ == '__main__':
    main()
