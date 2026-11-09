"""Display scaling adapters; statistical size factors are computed only by R/DESeq2."""
import csv
import math
from pathlib import Path
import subprocess


def table(path, required):
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if not reader.fieldnames or not set(required)<=set(reader.fieldnames) or len(set(reader.fieldnames))!=len(reader.fieldnames):
            raise ValueError('Missing/duplicate TSV columns: '+', '.join(required))
        rows=list(reader)
    if any(None in row or any(v is None for v in row.values()) for row in rows):
        raise ValueError('Malformed TSV fields')
    return rows


def supplied_values(path, field):
    values={}
    for row in table(path, ['sample_id',field]):
        name=row['sample_id']
        if not name or name in values: raise ValueError('Duplicate/empty factor sample ID')
        value=float(row[field])
        if not math.isfinite(value) or value<=0: raise ValueError('Supplied factors must be finite and positive')
        values[name]=value
    if not values: raise ValueError('Empty factors file')
    return values


def estimate_r(sites, cohort, controls, method, minimum, output, rscript, repo, log):
    """Serialize original integer data; R owns estimation, including controlGenes."""
    counts=output/'estimation_counts.tsv'
    keys=sorted(sites)
    with counts.open('w',newline='') as fh:
        writer=csv.writer(fh,delimiter='\t',lineterminator='\n')
        writer.writerow(['feature_id',*cohort])
        for i,key in enumerate(keys): writer.writerow([f'f{i+1}',*[sites[key][s] for s in cohort]])
    controls_path='-'
    if controls is not None:
        control_file=output/'control_indices.txt'
        index={key:i+1 for i,key in enumerate(keys)}
        control_file.write_text(''.join(str(index[key])+'\n' for key in controls if key in index))
        controls_path=str(control_file)
    factors=output/'r_size_factors.tsv'; diagnostics=output/'r_diagnostics.tsv'
    args=[rscript,'--vanilla',str(repo/'workflow/scripts/deseq2_factors.R'),
          str(repo/'R/deseq2_normalization.R'),str(counts),controls_path,method,str(minimum),str(factors),str(diagnostics)]
    import json
    log.write(json.dumps({'command':args})+'\n'); log.flush()
    try:
        result=subprocess.run(args,stdout=log,stderr=log)
    except OSError as error:
        log.write(json.dumps({'exit_status':None,'error':str(error)})+'\n'); log.flush()
        raise ValueError('Rscript executable unavailable: '+str(rscript)+
                         '; select --rscript from an isolated R/DESeq2 environment; see commands.log') from error
    log.write(json.dumps({'exit_status':result.returncode})+'\n'); log.flush()
    if result.returncode: raise ValueError('R/DESeq2 estimation failed; see commands.log for dependency or informative-feature errors; verify --rscript and its DESeq2 installation (no estimator fallback)')
    values=supplied_values(factors,'size_factor')
    if set(values)!=set(cohort): raise ValueError('R returned mismatching sample IDs')
    return values,table(diagnostics,['sample_id','method','contributing_features'])
