#!/usr/bin/env python3
"""Create single-base, strand-specific BigWigs with integrated display normalization."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

from bam_count import prepare_bam, reference_sizes, sample_names, policy
from core import matrix, sha256
from compare_matrices import read_matrix
from bigwig_normalization import estimate_r, supplied_values, table

REPO=Path(__file__).resolve().parents[2]


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def sizes_file(path):
    result={}
    for line in Path(path).read_text().splitlines():
        row=line.split('\t')
        if len(row)!=2 or not row[0] or row[0] in result or not row[1].isdigit() or int(row[1])<1:
            raise ValueError('Chromosome sizes require unique identifiers and positive integer lengths')
        result[row[0]]=int(row[1])
    if not result: raise ValueError('Empty reference sizes')
    return result


def column_hash(sites,sample):
    return digest([[*key,sites[key][sample]] for key in sorted(sites) if sites[key][sample]])


def write_bigwig(path,sizes,entries,version_info):
    import pyBigWig
    version_info['pyBigWig']=pyBigWig.__version__
    grouped={}
    for chrom,pos,value in entries: grouped.setdefault(chrom,[]).append((pos,value))
    with pyBigWig.open(str(path),'w') as bw:
        bw.addHeader(list(sizes.items()))
        # Explicit zero sentinel for an otherwise empty strand, recorded in metadata.
        if not entries: grouped={next(iter(sizes)):[(0,0.0)]}
        for chrom in sizes:
            values=grouped.get(chrom,[])
            if values:
                bw.addEntries([chrom]*len(values),[p for p,v in values],
                              ends=[p+1 for p,v in values],values=[float(v) for p,v in values])


def parser():
    p=argparse.ArgumentParser(description=__doc__,epilog=(
        'Matrix values must be original unthresholded integers. RPM uses all accepted alignments '
        'on both strands/all input references. DESeq2/controls use actual R; simple modes never invoke R. '
        'No smoothing, strand negation, background subtraction or coverage extension. New output directory required.'))
    sources=p.add_mutually_exclusive_group(required=True)
    sources.add_argument('--matrix',help='FivePrime Count 2 unthresholded TSV')
    sources.add_argument('--bam',action='append',help='Repeat in sample order; use --sample-name once per BAM')
    sources.add_argument('--bam-manifest',help='TSV with exactly sample_id and bam columns')
    p.add_argument('--sample-name',action='append')
    reference=p.add_mutually_exclusive_group(required=True)
    reference.add_argument('--reference',help='Read-only reference FASTA')
    reference.add_argument('--chrom-sizes',help='Read-only two-column reference sizes TSV')
    p.add_argument('--reference-sha256',help='Optional required FASTA checksum')
    p.add_argument('--normalization',required=True,choices=['raw','rpm','deseq2','controls','supplied'])
    p.add_argument('--cohort',nargs='+',help='Explicit common estimation cohort, in desired order; supplied mode can use saved factor table order')
    p.add_argument('--emit-samples',nargs='+',help='Subset to render; never changes estimated factors')
    p.add_argument('--strands',choices=['plus','minus','both'],default='both')
    p.add_argument('--assert-unthresholded',action='store_true',help='Required for matrix input: caller certifies no features were removed')
    p.add_argument('--matrix-complete-both-strands',action='store_true',help='Certify full-reference both-strand union for matrix-derived RPM depth')
    p.add_argument('--depths',help='Explicit RPM denominator TSV: sample_id, accepted_records, scope_id')
    p.add_argument('--mapq',type=int,default=20)
    p.add_argument('--exclude-duplicates',action='store_true')
    p.add_argument('--deseq2-type',choices=['ratio','poscounts'],help='R estimator (default ratio); no automatic fallback')
    p.add_argument('--controls',help='TSV: chromosome, strand, start, end (1-based single-base controls)')
    p.add_argument('--min-control-features',type=int,help='Controls mode only: minimum informative controls per cohort sample (default 10; at least 2)')
    p.add_argument('--size-factors',help='TSV with sample_id and size_factor or multiplier')
    p.add_argument('--factor-kind',choices=['size-factor','multiplier'])
    p.add_argument('--factor-metadata',help='Saved provenance JSON for verified single-sample regeneration')
    p.add_argument('--rscript',default='Rscript',help='Rscript executable; used only for deseq2/controls')
    p.add_argument('--output-dir',required=True,help='New external directory; existing paths are rejected')
    return p


def validate_options(a):
    """Shared direct/settings option checks; no data processing or estimation."""
    if a.min_control_features is not None and a.normalization!='controls':
        raise ValueError('--min-control-features applies only to controls normalization')
    if a.matrix and not a.assert_unthresholded: raise ValueError('Matrix input requires --assert-unthresholded; thresholded matrices are invalid')
    if bool(a.controls)!=(a.normalization=='controls'): raise ValueError('--controls is required only for controls mode')
    if a.normalization=='controls' and a.min_control_features is not None and a.min_control_features<2: raise ValueError('At least two informative controls required; default minimum is 10')
    if a.normalization=='supplied':
        if not a.size_factors or not a.factor_kind: raise ValueError('supplied requires --size-factors and explicit --factor-kind')
    elif any([a.size_factors,a.factor_kind,a.factor_metadata]): raise ValueError('Factor inputs require supplied mode')
    if a.depths and a.normalization!='rpm': raise ValueError('--depths applies only to rpm')
    if a.deseq2_type and a.normalization not in ('deseq2','controls'): raise ValueError('--deseq2-type applies only to R normalization')
    if a.matrix_complete_both_strands and (not a.matrix or a.normalization!='rpm'):
        raise ValueError('--matrix-complete-both-strands applies only to matrix RPM')
    if a.sample_name and not a.bam: raise ValueError('--sample-name applies only to repeated --bam inputs')
    if a.reference_sha256 and not a.reference: raise ValueError('--reference-sha256 requires FASTA')
    if not a.cohort and a.normalization!='supplied': raise ValueError('Declare one --cohort for all outputs')


def execute(a):
    validate_options(a)
    if a.min_control_features is None: a.min_control_features=10
    policy(a.mapq,a.exclude_duplicates)
    bams,names=[],[]
    if a.bam_manifest:
        rows=table(a.bam_manifest,['sample_id','bam'])
        if rows and set(rows[0])!={'sample_id','bam'}: raise ValueError('BAM manifest requires exactly sample_id and bam')
        bams=[r['bam'] for r in rows]; names=[r['sample_id'] for r in rows]
    elif a.bam:
        bams=a.bam; names=a.sample_name or [Path(b).stem for b in bams]
    if bams: sample_names(names,bams)
    elif not a.matrix: raise ValueError('Empty BAM manifest')
    paths=[a.matrix,a.bam_manifest,a.reference,a.chrom_sizes,a.controls,a.size_factors,a.factor_metadata,a.depths,*bams]
    inputs=[Path(x).expanduser().resolve() for x in paths if x]
    for path in inputs:
        if not path.is_file(): raise ValueError('Input is not a readable file: '+str(path))
    out=Path(a.output_dir).expanduser().resolve()
    if out==REPO or REPO in out.parents: raise ValueError('Analysis outputs must be outside the source checkout')
    if out.exists(): raise ValueError('output-dir already exists; choose a new directory')
    if any(path==out or out in path.parents for path in inputs): raise ValueError('Inputs cannot be inside output-dir')
    before={str(path):sha256(path) for path in inputs}
    out.mkdir(parents=True,exist_ok=False)
    meta={'status':'running','normalization':a.normalization,'parameters':vars(a),'inputs_sha256':before,'versions':{'python':sys.version},
          'units':None,'coordinate_definition':'1-based start=end -> 0-based [start-1,start)',
          'strand_definition':'alignment strand relative to each unchanged reference contig; both tracks nonnegative',
          'min_total_count':0,'warnings':[],
          'code_sha256':{str(p.relative_to(REPO)):sha256(p) for p in [Path(__file__),REPO/'workflow/scripts/bigwig_normalization.py',REPO/'workflow/scripts/bam_count.py',REPO/'workflow/scripts/core.py',REPO/'R/deseq2_normalization.R',REPO/'workflow/scripts/deseq2_factors.R']}}
    try:
        import pyBigWig
        meta['versions']['pyBigWig']=pyBigWig.__version__
        if a.reference: reference_sizes(a.reference,out/'reference.sizes.tsv',a.reference_sha256)
        else: (out/'reference.sizes.tsv').write_bytes(Path(a.chrom_sizes).read_bytes())
        sizes=sizes_file(out/'reference.sizes.tsv')
        meta['reference_dictionary_sha256']=digest(list(sizes.items()))
        with (out/'commands.log').open('w') as log:
            depths={}
            if bams:
                accepted=[]; preparation=[]
                for i,bam in enumerate(bams):
                    dest=out/'accepted'/f'{i}.bam'
                    preparation.append(prepare_bam(bam,dest,out/'reference.sizes.tsv',a.mapq,a.exclude_duplicates,log))
                    cmd=['samtools','view','-c',str(dest)]
                    log.write(json.dumps({'command':cmd})+'\n'); log.flush()
                    depths[names[i]]=int(subprocess.check_output(cmd,text=True,stderr=log))
                    log.write(json.dumps({'exit_status':0})+'\n')
                    accepted.append(dest)
                matrix(accepted,names,out/'reference.sizes.tsv',out/'counts.raw.tsv','both',log)
                source=out/'counts.raw.tsv'
                meta['preparation']=preparation
                meta['versions'].update({t:subprocess.check_output([t,'--version'],text=True).splitlines()[0] for t in ('samtools','bedtools')})
                meta['resolved_executables']={t:shutil.which(t) for t in ('samtools','bedtools')}
            else: source=Path(a.matrix)
            names,sites,totals=read_matrix(source,sizes,0)
            sample_names(names,[f'column_{i}' for i in range(len(names))])
            meta['matrix_sha256']=sha256(source)
            meta['matrix_state']='unthresholded, generated' if bams else 'unthresholded, caller asserted'
            meta['sample_count_hashes']={name:column_hash(sites,name) for name in names}
            if bams and any(depths[n]!=totals[n] for n in names): raise ValueError('Accepted alignment denominator differs from complete 5-prime union')
            supplied={}
            if a.normalization=='supplied': supplied=supplied_values(a.size_factors,'size_factor' if a.factor_kind=='size-factor' else 'multiplier')
            cohort=a.cohort or list(supplied)
            sample_names(cohort,[f'cohort_{i}' for i in range(len(cohort))])
            if a.normalization!='supplied' and any(n not in names for n in cohort): raise ValueError('Cohort sample absent from input matrix')
            emitted=a.emit_samples or [n for n in cohort if n in names]
            if not emitted or len(set(emitted))!=len(emitted) or any(n not in cohort or n not in names for n in emitted):
                raise ValueError('Emitted samples must be unique and present in both cohort and inputs')
            meta.update(cohort=cohort,emitted_samples=emitted,feature_count=len(sites))
            meta['cohort_id']=digest(cohort)
            sf={n:1.0 for n in cohort}; multipliers={n:1.0 for n in cohort}
            if a.normalization=='rpm':
                scope='accepted alignments in all input references, both alignment strands'
                if a.depths:
                    rows=table(a.depths,['sample_id','accepted_records','scope_id']); depths={}; scopes=set()
                    for row in rows:
                        n=row['sample_id']
                        if n in depths or not row['accepted_records'].isdigit() or not row['scope_id']: raise ValueError('Invalid/duplicate depth row')
                        depths[n]=int(row['accepted_records']); scopes.add(row['scope_id'])
                    if len(scopes)!=1: raise ValueError('RPM depths must share one explicitly declared scope_id')
                    scope='externally declared accepted alignment pool: '+next(iter(scopes))
                elif not bams:
                    if not a.matrix_complete_both_strands: raise ValueError('Matrix RPM requires --matrix-complete-both-strands or explicit --depths')
                    depths=totals
                if any(n not in depths or depths[n]<=0 for n in cohort): raise ValueError('RPM needs positive accepted depth for every cohort sample')
                multipliers={n:1e6/depths[n] for n in cohort}
                meta['denominator_definition']=scope
                meta['units']='accepted alignment 5-prime ends per million in declared depth scope'
            elif a.normalization in ('deseq2','controls'):
                controls=None
                if a.normalization=='controls':
                    controls=[]
                    for row in table(a.controls,['chromosome','strand','start','end']):
                        if set(row)!={'chromosome','strand','start','end'}: raise ValueError('Controls require exactly four coordinate columns')
                        if not row['start'].isdigit() or not row['end'].isdigit() or int(row['start'])!=int(row['end']): raise ValueError('Control requires integer 1-based start=end')
                        key=(row['chromosome'],row['strand'],int(row['start']))
                        if key[0] not in sizes or key[1] not in ('+','-') or not 1<=key[2]<=sizes[key[0]]: raise ValueError('Control outside reference/invalid strand')
                        if key in controls: raise ValueError('Duplicate control coordinate')
                        controls.append(key)
                    missing=[k for k in controls if k not in sites]
                    zero=[k for k in controls if k in sites and all(sites[k][n]==0 for n in cohort)]
                    meta['control_diagnostics']={'requested':len(controls),'missing':missing,'all_zero':zero,
                                                 'minimum_informative_per_sample':a.min_control_features}
                    if missing: meta['warnings'].append('Missing controls treated as unobserved/all-zero and excluded; see diagnostics')
                    controls=[k for k in controls if k in sites]
                    if not controls: raise ValueError('No control positions matched the unthresholded matrix')
                work=out/'normalization'; work.mkdir()
                meta['estimator']=a.deseq2_type or 'ratio'
                sf,diagnostics=estimate_r(sites,cohort,controls,a.deseq2_type or 'ratio',
                    a.min_control_features if controls is not None else 1,work,a.rscript,REPO,log)
                meta['normalization_diagnostics']=diagnostics
                meta['estimator']=a.deseq2_type or 'ratio'
                if any(int(d['contributing_features'])<10 for d in diagnostics): meta['warnings'].append('Fewer than ten informative features in at least one sample')
                multipliers={n:1.0/sf[n] for n in cohort}
                meta['units']='DESeq2 size-factor-normalized 5-prime counts (not RPM)'
                meta['resolved_Rscript']=shutil.which(a.rscript) or a.rscript
            elif a.normalization=='supplied':
                if any(n not in supplied for n in cohort): raise ValueError('Missing cohort factor')
                if a.factor_kind=='size-factor':
                    sf={n:supplied[n] for n in cohort}; multipliers={n:1.0/sf[n] for n in cohort}
                else: multipliers={n:supplied[n] for n in cohort}
                meta['units']='externally size-factor-normalized 5-prime counts' if a.factor_kind=='size-factor' else 'external multiplier-scaled 5-prime counts'
                if a.factor_metadata:
                    saved=json.loads(Path(a.factor_metadata).read_text())
                    if saved.get('outputs_sha256',{}).get('normalization_factors.tsv')!=sha256(a.size_factors):
                        raise ValueError('Saved factor file checksum mismatch')
                    if saved['status']!='complete' or saved['reference_dictionary_sha256']!=meta['reference_dictionary_sha256'] or saved['cohort']!=cohort:
                        raise ValueError('Saved factor cohort/reference context mismatch')
                    for n in emitted:
                        if saved['sample_count_hashes'].get(n)!=meta['sample_count_hashes'][n]: raise ValueError('Saved factor input count context mismatch')
                    meta['factor_context_verified']=True
                    meta['units']=saved['units']
                    meta['factor_origin_normalization']=saved['normalization']
                    meta['factor_origin_estimator']=saved.get('estimator')
                else:
                    meta['warnings'].append('Supplied factor origin/context is user-declared, not verified')
            else: meta['units']='raw accepted 5-prime counts'
            if any(not math.isfinite(v) or v<=0 for v in multipliers.values()): raise ValueError('Invalid normalization multiplier')
            with (out/'normalization_factors.tsv').open('w',newline='') as fh:
                writer=csv.writer(fh,delimiter='\t',lineterminator='\n')
                writer.writerow(['sample_id','size_factor','multiplier','accepted_records','method','units','cohort_id'])
                for n in cohort: writer.writerow([n,sf[n] if a.normalization in ('deseq2','controls') or (a.normalization=='supplied' and a.factor_kind=='size-factor') else '',multipliers[n],depths.get(n,''),a.normalization,meta['units'],meta['cohort_id']])
            tracks=out/'tracks'; tracks.mkdir(); meta['empty_tracks']=[]
            for n in emitted:
                for sign in ('+','-') if a.strands=='both' else (('+' if a.strands=='plus' else '-'),):
                    entries=[]
                    for key in sorted(sites):
                        if key[1]==sign and sites[key][n]:
                            value=sites[key][n]*multipliers[n]
                            if not math.isfinite(value) or value>3.4028234663852886e38 or value<1.401298464324817e-45: raise ValueError('Track value outside representable nonzero float32 range')
                            entries.append((key[0],key[2]-1,value))
                    name=n+('.plus.bw' if sign=='+' else '.minus.bw')
                    if not entries: meta['empty_tracks'].append(name)
                    write_bigwig(tracks/name,sizes,entries,meta['versions'])
        if any(sha256(Path(p))!=h for p,h in before.items()): raise ValueError('Input changed during execution')
        meta['normalization']=a.normalization
        meta['outputs_sha256']={str(p.relative_to(out)):sha256(p) for p in out.rglob('*') if p.is_file()}
        meta['status']='complete'
    except Exception as error:
        meta.update(status='failed',error=str(error)); raise
    finally:
        (out/'provenance.json').write_text(json.dumps(meta,indent=2)+'\n')


def main():
    p=parser(); a=p.parse_args()
    try: execute(a)
    except (ValueError,OSError,ImportError,subprocess.CalledProcessError) as error: p.exit(1,f'bigwig2: {error}\n')


if __name__=='__main__': main()
