"""Snakemake adapter for the independently executable counting operations."""
import json
from pathlib import Path
import subprocess
import shutil
import sys

sys.path.insert(0, snakemake.scriptdir)
from bam_count import prepare_bam, reference_sizes
from core import matrix, sha256, threshold
from compare_matrices import compare

s = dict(snakemake.params.settings)
for path in [*snakemake.output, *snakemake.log]:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
with open(snakemake.log[0], 'w') as log:
    step = snakemake.params.step
    try:
        if step == 'sizes':
            result = reference_sizes(s['fasta'], snakemake.output[0], s['reference_sha256'])
        elif step == 'prepare':
            result = prepare_bam(snakemake.input.bam, snakemake.output[0], snakemake.input.sizes,
                                 s['mapq'], s['exclude_duplicates'], log)
            result['samtools_version'] = subprocess.check_output(['samtools', '--version'], text=True).splitlines()[0]
        elif step == 'count':
            result = matrix(snakemake.input.bams, list(s['samples']), snakemake.input.sizes,
                            snakemake.output[0], s['strand'], log)
        elif step == 'threshold':
            threshold(snakemake.input.matrix, snakemake.output[0], s['minimum'])
            result = {'minimum': s['minimum']}
        elif step == 'compare':
            result = compare(snakemake.input.expected, snakemake.input.observed, snakemake.input.sizes,
                             s['minimum'], s['comparison'].get('expected_features'),
                             s['comparison']['expected_matrix_sha256'])
            Path(snakemake.output[0]).write_text(json.dumps(result, indent=2)+'\n')
            # Keep the discrepancy report as a log as well: Snakemake removes failed outputs.
            log.write(json.dumps(result)+'\n'); log.flush()
            if not result['equal']:
                raise ValueError('Historical matrix mismatch; inspect comparison log; do not correct silently')
        elif step == 'provenance':
            import snakemake as smk
            result = {'settings': s, 'python': sys.version, 'snakemake': smk.__version__,
                      'reference_sha256': sha256(s['fasta']),
                      'manifest_sha256': sha256(s['manifest']),
                      'code': {str(Path(p).relative_to(Path(snakemake.scriptdir).parents[1])): sha256(p)
                               for p in snakemake.input.code},
                      'preparations': [Path(p).read_text() for p in snakemake.input.tool_logs],
                      'count_log': Path(snakemake.input.count_log).read_text(),
                      'executables': {tool: shutil.which(tool) for tool in ('samtools', 'bedtools')},
                      'output_sha256': sha256(snakemake.input.matrix)}
            Path(snakemake.output[0]).write_text(json.dumps(result, indent=2)+'\n')
        else:
            raise ValueError('Unknown step: '+step)
        log.write(json.dumps({'step': step, 'exit_status': 0, 'result': result})+'\n')
    except Exception as error:
        log.write(json.dumps({'step': step, 'exit_status': 1, 'error': str(error)})+'\n')
        raise
