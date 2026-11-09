"""Thin Snakemake adapter; never counts BAM 5-prime ends."""
import csv
from pathlib import Path
import subprocess
import sys
sys.path.insert(0,snakemake.scriptdir)
from workflow_bigwig import run
from core import sha256

if snakemake.params.step == 'depths':
    output=Path(snakemake.output[0]); output.parent.mkdir(parents=True,exist_ok=True)
    logpath=Path(snakemake.log[0]); logpath.parent.mkdir(parents=True,exist_ok=True)
    with logpath.open('w') as log, output.open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n')
        w.writerow(['sample_id','accepted_records','scope_id'])
        for name,bam in zip(snakemake.params.names,snakemake.input.bams):
            before=sha256(bam)
            cmd=['samtools','view','-c',bam]
            log.write(str(cmd)+'\n');log.flush()
            depth=int(subprocess.check_output(cmd,text=True,stderr=log))
            if sha256(bam)!=before: raise ValueError('Accepted BAM changed during depth audit')
            w.writerow([name,depth,snakemake.params.space+':accepted-all-reference-both-strands'])
            log.write(f'{name}\t{before}\t{depth}\texit_status=0\n')
else:
    run(dict(snakemake.params.plan), Path(snakemake.params.repo))
