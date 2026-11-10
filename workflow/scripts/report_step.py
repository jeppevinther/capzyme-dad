"""Final Snakemake adapter; the standalone command uses this same generator."""
import sys
from pathlib import Path
sys.path.insert(0, snakemake.scriptdir)
from run_report import generate_report
code = generate_report(snakemake.params.analysis_dir, context=dict(snakemake.params.context))
Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
Path(snakemake.log[0]).write_text('Record-only reporting exit status: '+str(code)+'\n'
    'This task cannot certify the subsequent Snakemake exit status.\n')
if code:
    # Snakemake may remove declared failed outputs. Keep diagnostic copies as logs.
    import shutil
    root=Path(snakemake.params.analysis_dir)/'report'
    for name in ('report.html','sample_summary.tsv','report_status.json'):
        shutil.copy2(root/name, Path(snakemake.log[0]).with_name('failed_'+name))
    raise ValueError('Required report checks failed; preserved diagnostics beside report.log')
