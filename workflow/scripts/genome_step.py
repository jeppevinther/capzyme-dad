"""Snakemake adapter for standalone mapping, shared counting and provenance."""
import json
from pathlib import Path
import sys
import subprocess

sys.path.insert(0, snakemake.scriptdir)
from bam_count import command
from core import matrix, sha256
from sequential import accounting, publish_stage

s = snakemake.params.settings
repo = Path(snakemake.params.repo)
for path in [*snakemake.output, *snakemake.log]:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
with open(snakemake.log[0], 'w') as log:
    step = snakemake.params.step
    # Reserve a CPU slot for the active Python/I/O helper while the aligner or
    # sort runs. Sorting is sequential inside each standalone tool.
    tool_threads = max(1, int(snakemake.threads)-1)
    if step == 'rna':
        args = [sys.executable, str(repo/'workflow/scripts/rna_handoff.py'),
                '--fastq', str(snakemake.input.fastq), '--reference', str(snakemake.input.reference),
                '--index', s['rna_index'], '--threads', str(tool_threads),
                '--mapq', str(s['mapq']), '--seed', str(s['seed'])]
        if s['norc']: args.append('--norc')
        if s['flags'] & 1024: args.append('--exclude-duplicates')
        publish_stage(s['output'], 'rna', snakemake.wildcards.sample, args, log)
    elif step == 'genome':
        g = s['genome']
        if g.get('verification_manifest'):
            from reference_hashes import check_manifest
            check_manifest(g['verification_manifest'],g['reference'],g['index_prefix'],g['aligner'])
        args = [sys.executable, str(repo/'workflow/scripts'/f"{g['aligner']}_genome.py"), 'map',
                '--fastq', str(snakemake.input.fastq), '--index', g['index_prefix'],
                '--reference', g['reference'], '--threads', str(tool_threads),
                '--mapq', str(g['mapq']), '--seed', str(g['seed']), '--uncertain-cigar', 'reject']
        args += ['--decoys', g['decoys']] if g['decoys'] else ['--no-decoys']
        if g['reference_sha256']: args += ['--reference-sha256', g['reference_sha256']]
        if g['exclude_duplicates']: args.append('--exclude-duplicates')
        publish_stage(s['output'], 'genome', snakemake.wildcards.sample, args, log)
    elif step == 'genome_count':
        matrix(snakemake.input.bams, list(s['samples']), snakemake.input.sizes,
               snakemake.output[0], s['genome']['strand'], log)
    elif step == 'accounting':
        accounting(snakemake.params.rna, snakemake.params.genome, Path(snakemake.output[0]), log,
                   snakemake.input.trim)
    elif step == 'provenance':
        result = {'settings': dict(s), 'sample_order': list(s['samples']),
                  'reference_spaces': 'Separate RNA/genome integer matrices; reference-relative strands; 1-based start=end',
                  'inputs_sha256': {p: sha256(p) for p in snakemake.input if Path(p).is_file()},
                  'workflow_code_sha256': {str(Path(p).relative_to(repo)): sha256(p) for p in snakemake.params.code}}
        if (repo/'.git').exists():
            result['git'] = {key: subprocess.check_output(['git','-C',str(repo),*args],text=True).strip()
                             for key,args in [('commit',['rev-parse','HEAD']),('branch',['branch','--show-current']),
                                              ('dirty_status',['status','--porcelain'])]}
        else:
            snapshot = repo/'source_snapshot.json'
            if not snapshot.is_file():
                raise ValueError('Execution snapshot requires source_snapshot.json')
            result['source_snapshot'] = json.loads(snapshot.read_text())
            result['source_snapshot_sha256'] = sha256(snapshot)
        Path(snakemake.output[0]).write_text(json.dumps(result, indent=2)+'\n')
    else:
        raise ValueError('Unknown workflow step '+step)
