# Final record-only reporting: no upstream rule depends on these outputs.
from run_report import workflow_context
REPORT_CONTEXT = workflow_context(SETTINGS, config, TRACKS, ANALYSIS_OUTPUTS, REPORT_ROUTE)
rule run_report:
    input:
        analysis=ANALYSIS_OUTPUTS,
        code=[str(REPO/'workflow/scripts/run_report.py')]
    output:
        html=OUT+'/report/report.html',
        summary=OUT+'/report/sample_summary.tsv',
        status=OUT+'/report/report_status.json',
        context=OUT+'/report/analysis_context.json'
    params: context=REPORT_CONTEXT, analysis_dir=OUT
    threads: 1
    resources: mem_mb=512, runtime=5
    log: OUT+'/logs/report.log'
    script: '../scripts/report_step.py'
