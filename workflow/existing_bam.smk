"""Existing RNA-reference BAMs -> filtered/sorted copies -> 5′ matrices only."""
from pathlib import Path
import sys
from snakemake.utils import min_version
min_version('9.9.0')
REPO = Path(workflow.basedir).parent.resolve()
sys.path.insert(0, str(REPO / 'workflow/scripts'))
from bam_configuration import load_bam_settings
S = load_bam_settings(config, REPO)
SETTINGS = S
from workflow_bigwig import plans
TRACKS = plans(config.get('bigwigs'), S, REPO, bam=True)
TRACK_OUTPUTS = [f for p in TRACKS.values() for f in p['files']]
OUT = S['output']
NAMES = list(S['samples'])
SIZES = OUT + '/reference/sizes.tsv'
RAW = OUT + '/counts/5ends.raw.tsv'
FINAL = OUT + '/counts/5ends.min_total_' + str(S['minimum']) + '.tsv'
CODE = [str(REPO/'workflow/scripts'/name) for name in ('bam_count.py','bam_configuration.py','reference_hashes.py','configuration.py','core.py','compare_matrices.py')]

wildcard_constraints:
    sample='[A-Za-z0-9][A-Za-z0-9_.-]*'

ANALYSIS_OUTPUTS = [FINAL, OUT+'/provenance/run.json', *TRACK_OUTPUTS,
    *([OUT+'/comparison/exact.json'] if S['comparison'] else [])]
REPORT_ROUTE = 'bam'
rule all:
    input: ANALYSIS_OUTPUTS, OUT+'/report/report.html', OUT+'/report/sample_summary.tsv', OUT+'/report/report_status.json', OUT+'/report/analysis_context.json'

rule reference_sizes:
    input: fasta=S['fasta'], code=CODE
    output: SIZES
    params: step='sizes', settings=S
    log: OUT + '/logs/reference.jsonl'
    threads: 1
    resources: mem_mb=512
    conda: 'envs/python.yaml'
    script: 'scripts/bam_step.py'

rule prepare_existing_bam:
    input: bam=lambda wc: S['samples'][wc.sample], sizes=SIZES, code=CODE
    output: OUT + '/bam/{sample}.accepted.bam'
    params: step='prepare', settings=S
    log: OUT + '/logs/{sample}.prepare.jsonl'
    threads: 1
    resources: mem_mb=S['resources'].get('filter_mem_mb', 1024)
    conda: 'envs/bam_prepare.yaml'
    script: 'scripts/bam_step.py'

rule count_existing_bams:
    input: bams=expand(OUT+'/bam/{sample}.accepted.bam', sample=NAMES), sizes=SIZES, code=CODE
    output: RAW
    params: step='count', settings=S
    log: OUT + '/logs/count.jsonl'
    threads: S['resources'].get('count_threads', 2)
    resources: mem_mb=S['resources'].get('count_mem_mb', 8192)
    conda: 'envs/counting.yaml'
    script: 'scripts/bam_step.py'

rule threshold_existing_counts:
    input: matrix=RAW, code=CODE
    output: FINAL
    params: step='threshold', settings=S
    log: OUT + '/logs/threshold.jsonl'
    threads: 1
    resources: mem_mb=512
    conda: 'envs/python.yaml'
    script: 'scripts/bam_step.py'

rule existing_provenance:
    input:
        matrix=FINAL, manifest=S['manifest'], reference=S['fasta'], code=CODE,
        tool_logs=expand(OUT+'/logs/{sample}.prepare.jsonl', sample=NAMES),
        count_log=OUT+'/logs/count.jsonl'
    output: OUT+'/provenance/run.json'
    params: step='provenance', settings=S
    log: OUT+'/logs/provenance.jsonl'
    threads: 1
    resources: mem_mb=512
    script: 'scripts/bam_step.py'

if S['comparison']:
    rule compare_historical_matrix:
        input: observed=FINAL, expected=S['comparison']['expected_matrix'], sizes=SIZES, code=CODE
        output: OUT+'/comparison/exact.json'
        params: step='compare', settings=S
        log: OUT+'/logs/compare.jsonl'
        threads: 1
        resources: mem_mb=2048
        conda: 'envs/python.yaml'
        script: 'scripts/bam_step.py'

include: 'rules/tracks.smk'

include: 'rules/report.smk'
