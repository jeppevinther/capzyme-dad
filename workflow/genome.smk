"""Genome-only orchestration for already trimmed FASTQs and an existing index."""
from pathlib import Path
import sys
from snakemake.utils import min_version

min_version('9.9.0')
REPO = Path(workflow.basedir).parent.resolve()
sys.path.insert(0, str(REPO/'workflow/scripts'))
from genome_configuration import genome_only_settings

SETTINGS = genome_only_settings(config, REPO)
from workflow_bigwig import plans
TRACKS = plans(config.get('bigwigs'), SETTINGS, REPO)
TRACK_OUTPUTS = [f for p in TRACKS.values() for f in p['files']]
GENOME = SETTINGS['genome']
GENOME_ONLY = True
OUT = SETTINGS['output']
SAMPLES = list(SETTINGS['samples'])
CODE = [str(REPO/'workflow/scripts'/name) for name in (
    'reference_hashes.py','configuration.py','genome_configuration.py','fixed_reference.py','core.py',
    'bam_count.py','compare_matrices.py','genome_common.py',
    'rna_handoff.py','sequential.py','bowtie2_genome.py','hisat2_genome.py')]

wildcard_constraints:
    sample='[A-Za-z0-9][A-Za-z0-9_.-]*'

ANALYSIS_OUTPUTS = [*TRACK_OUTPUTS, OUT+'/counts/genome.5ends.raw.tsv', OUT+'/provenance/genome.json']
REPORT_ROUTE = 'genome'
rule all:
    input: ANALYSIS_OUTPUTS, OUT+'/report/report.html', OUT+'/report/sample_summary.tsv', OUT+'/report/report_status.json', OUT+'/report/analysis_context.json'

include: 'rules/genome.smk'

include: 'rules/tracks.smk'

include: 'rules/report.smk'
