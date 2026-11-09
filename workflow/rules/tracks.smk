# Both entry points expose SETTINGS, TRACKS and REPO. Per-plan rules have fixed
# outputs, making cohort/method additions independent of existing track jobs.
for space in ('rna', 'genome'):
    depth_plans = [p for key,p in TRACKS.items() if key.startswith(space+'/') and p['depth']]
    if depth_plans:
        dp = depth_plans[0]
        rule:
            name: 'bigwig_depths_'+space
            input: bams=dp['bams'], code=[str(REPO/'workflow/scripts/bigwig_step.py')]
            output: dp['depth']
            params: step='depths', names=dp['names'], space=space
            log: OUT+'/logs/bigwig.'+space+'.depths.log'
            threads: 1
            resources: mem_mb=512
            script: '../scripts/bigwig_step.py'

for key, plan in TRACKS.items():
    rule:
        name: 'bigwig_'+key.replace('/','_')
        input:
            matrix=plan['args'][1],
            reference=plan['args'][4],
            external=plan['external'],
            depths=([plan['depth']] if plan['depth'] else []),
            code=[str(REPO/f) for f in ('workflow/scripts/bigwig2.py', 'workflow/scripts/workflow_bigwig.py',
                  'workflow/scripts/bigwig_normalization.py', 'workflow/scripts/bam_count.py',
                  'workflow/scripts/core.py', 'workflow/scripts/compare_matrices.py',
                  'workflow/scripts/deseq2_factors.R', 'R/deseq2_normalization.R')]
        output: plan['files']
        params: step='tracks', plan=plan, repo=str(REPO)
        threads: 1
        resources: mem_mb=plan['mem_mb']
        script: '../scripts/bigwig_step.py'
