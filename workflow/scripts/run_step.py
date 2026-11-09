"""Snakemake script adapter. All generated artifacts stay in output_dir."""
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, snakemake.scriptdir)
from core import collapse, combine, matrix, sha256, threshold
from fixed_reference import verify_fixed

for output in snakemake.output:
    Path(output).parent.mkdir(parents=True, exist_ok=True)
for log in snakemake.log:
    Path(log).parent.mkdir(parents=True, exist_ok=True)

with open(snakemake.log[0], 'w') as log:
    step = snakemake.params.step
    try:
        if step == 'verify_fixed':
            summary = verify_fixed(snakemake.params.reference, snakemake.output.sizes, snakemake.output.report)
        elif step == 'collapse':
            summary = collapse(snakemake.input.fasta, snakemake.params.prefix,
                               snakemake.params.k, snakemake.params.drops,
                               snakemake.input.legacy)
        elif step == 'combine':
            combine(snakemake.input.fastas, snakemake.output[0])
            summary = {'sources': [{'path': x, 'sha256': sha256(x)} for x in snakemake.input.fastas]}
        elif step == 'matrix':
            summary = matrix(snakemake.input.bams, snakemake.params.names,
                             snakemake.input.fai, snakemake.output[0],
                             snakemake.params.strand, log)
        elif step == 'threshold':
            threshold(snakemake.input.matrix, snakemake.output[0], snakemake.params.minimum)
            summary = {'min_total_count': snakemake.params.minimum}
        elif step == 'provenance':
            summary = {'settings': dict(snakemake.params.settings), 'versions': {}, 'inputs': {}, 'code': {}}
            for name, path in zip(('cutadapt', 'bowtie2', 'samtools', 'bedtools'), snakemake.input.tool_logs):
                summary['versions'][name] = Path(path).read_text().splitlines()[0]
            import snakemake as smk
            summary['versions']['snakemake'] = smk.__version__
            summary['versions']['python'] = sys.version
            for path in snakemake.params.sources:
                summary['inputs'][path] = sha256(path)
            for path in snakemake.input.code:
                summary['code'][str(Path(path).relative_to(snakemake.params.repo))] = sha256(path)
            def git(*args):
                return subprocess.check_output(['git', '-C', snakemake.params.repo, *args], text=True).strip()
            if (Path(snakemake.params.repo) / '.git').exists():
                summary['git'] = {'commit': git('rev-parse', 'HEAD'), 'branch': git('branch', '--show-current'),
                                  'dirty': bool(git('status', '--porcelain'))}
            else:
                snapshot = Path(snakemake.params.repo) / 'source_snapshot.json'
                if not snapshot.is_file():
                    raise ValueError('Execution snapshot requires source_snapshot.json provenance')
                snapshot_data = json.loads(snapshot.read_text())
                if 'git' in snapshot_data:
                    summary['git'] = snapshot_data['git']
                elif snapshot_data.get('kind') == 'installed-wheel':
                    summary['installed_source'] = snapshot_data
                else:
                    raise ValueError('Unrecognized source snapshot provenance')
                summary['execution_snapshot_sha256'] = sha256(snapshot)
            Path(snakemake.output[0]).write_text(json.dumps(summary, indent=2) + '\n')
        else:
            raise ValueError(f'Unknown step {step}')
        log.write(json.dumps({'step': step, 'exit_status': 0, 'summary': summary}, indent=2) + '\n')
    except Exception:
        log.write(json.dumps({'step': step, 'exit_status': 1}) + '\n')
        raise
