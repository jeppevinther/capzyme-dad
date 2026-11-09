"""Optional matrix-only workflow coordination. Scientific work stays in bigwig2."""
import importlib.util
import importlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import uuid

from configuration import keys, integer, checked_path
from bigwig2 import parser, validate_options
from bigwig_normalization import table, supplied_values

METHODS = ('raw', 'rpm', 'deseq2', 'controls', 'supplied')


def plans(value, settings, repo, bam=False, dependencies=True):
    if value is None: return {}
    keys(value, ['rna', 'genome'], 'bigwigs')
    result = {}
    for space, block in value.items():
        keys(block, ['methods', 'cohort', 'emit_samples', 'strands', 'mem_mb'], 'bigwigs.'+space)
        if space == 'genome' and not settings.get('genome'):
            raise ValueError('Genome BigWigs require a genome route')
        if space == 'rna' and 'strand' not in settings:
            raise ValueError('RNA BigWigs require an RNA or existing-BAM route')
        methods = block.get('methods', {})
        keys(methods, METHODS, space+'.methods')
        if not methods: raise ValueError('Declare at least one BigWig method')
        names = list(settings['samples'])
        cohort = block.get('cohort', names)
        emitted = block.get('emit_samples', cohort)
        for label, seq, universe in [('cohort', cohort, names), ('emit_samples', emitted, cohort)]:
            if not isinstance(seq, list) or not seq or any(type(n) is not str or n not in universe for n in seq) or len(set(seq)) != len(seq):
                raise ValueError(label+' must be unique configured samples; emitted samples must belong to cohort')
        strand = block.get('strands', settings['genome']['strand'] if space=='genome' else settings['strand'])
        available = settings['genome']['strand'] if space=='genome' else settings['strand']
        if strand not in ('plus','minus','both') or available != 'both' and strand != available:
            raise ValueError('Requested BigWig strands are not present in the unthresholded matrix')
        out = settings['output']
        matrix = out+'/counts/'+('genome.5ends.raw.tsv' if space=='genome' else '5ends.raw.tsv')
        if space=='genome': sizes=out+'/genome/'+names[0]+'/reference.sizes.tsv'
        elif bam: sizes=out+'/reference/sizes.tsv'
        else: sizes=out+'/reference/fixed.sizes.tsv' if settings['reference_mode']=='fixed' else out+'/reference/combined.fa.fai'
        # FASTA .fai has five columns; use FASTA rather than the two-column sizes interface.
        reference = None
        if not bam and space=='rna' and settings['reference_mode']!='fixed':
            reference=out+'/reference/combined.fa'
        bams = [out+'/genome/'+n+'/accepted.bam' if space=='genome' else
                out+'/bam/'+n+'.accepted.bam' if bam or not settings.get('genome') else
                out+'/rna/'+n+'/accepted.bam' for n in names]
        for method, options in methods.items():
            keys(options, ['depths','deseq2_type','controls','min_control_features','size_factors',
                           'factor_kind','factor_metadata','rscript'], space+'.'+method)
            args=['--matrix',matrix,'--assert-unthresholded',
                  '--reference' if reference else '--chrom-sizes', reference or sizes,
                  '--normalization',method,'--cohort',*cohort,'--emit-samples',*emitted,'--strands',strand]
            inputs=[]
            for name, item in options.items():
                if name in ('depths','controls','size_factors','factor_metadata'):
                    item=str(checked_path(item, 'bigwigs.'+name)); inputs.append(item)
                    if Path(out)==Path(item) or Path(out) in Path(item).parents:
                        raise ValueError('BigWig factor/control inputs must be outside output_dir')
                if name=='min_control_features': integer(item, name, 2)
                if not isinstance(item, (str,int)) or isinstance(item,bool): raise ValueError('Invalid BigWig option '+name)
                args += ['--'+name.replace('_','-'),str(item)]
            if method in ('deseq2','controls'):
                if 'deseq2_type' not in options: raise ValueError('Select deseq2_type explicitly; no automatic estimator switch')
                if len(cohort)<2: raise ValueError('DESeq2 requires at least two cohort samples')
            auto_depth = method=='rpm' and not options.get('depths')
            depth_path=out+'/bigwigs/'+space+'/accepted_depths.tsv'
            if auto_depth: args+=['--depths',depth_path]
            final=out+'/bigwigs/'+space+'/'+method
            validate_options(parser().parse_args([*args,'--output-dir',final]))
            if method=='rpm' and options.get('depths'):
                rows=table(options['depths'],['sample_id','accepted_records','scope_id'])
                values={r['sample_id']:r for r in rows}
                if len(values)!=len(rows) or len({r['scope_id'] for r in rows})!=1 or any(n not in values or not values[n]['accepted_records'].isdigit() or int(values[n]['accepted_records'])<=0 or not values[n]['scope_id'] for n in cohort):
                    raise ValueError('RPM needs positive complete accepted-pool depths and one scope_id')
            if method=='supplied':
                values=supplied_values(options['size_factors'], 'size_factor' if options['factor_kind']=='size-factor' else 'multiplier')
                if any(n not in values for n in cohort): raise ValueError('Missing cohort factor')
            if method=='controls':
                rows=table(options['controls'],['chromosome','strand','start','end'])
                coords=[(r['chromosome'],r['strand'],r['start'],r['end']) for r in rows]
                if not rows or len(set(coords))!=len(coords) or any(set(r)!={'chromosome','strand','start','end'} or r['strand'] not in ('+','-') or not r['start'].isdigit() or int(r['start'])<1 or r['start']!=r['end'] for r in rows):
                    raise ValueError('Invalid control single-base positions')
            if dependencies:
                try: importlib.import_module('pyBigWig')
                except ImportError as exc: raise ValueError('Optional BigWigs require working pyBigWig; install the tracks extra') from exc
                if auto_depth and not shutil.which('samtools'): raise ValueError('Automatic accepted depths require samtools')
                if method in ('deseq2','controls'):
                    r=options.get('rscript','Rscript')
                    if not shutil.which(r): raise ValueError('BigWig R normalization requires Rscript and DESeq2; set rscript')
                    probe=subprocess.run([r,'--vanilla','-e', 'if (!requireNamespace("DESeq2", quietly=TRUE)) quit(status=1)'],capture_output=True,text=True,timeout=30,
                        env=dict(os.environ,**{k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')}))
                    if probe.returncode: raise ValueError('BigWig Rscript cannot load DESeq2: '+probe.stderr[-1500:])
            files=['provenance.json','normalization_factors.tsv','reference.sizes.tsv','commands.log','workflow.json']
            files += ['tracks/'+n+'.'+sign+'.bw' for n in emitted for sign in (('plus','minus') if strand=='both' else (strand,))]
            if method in ('deseq2','controls'):
                files += ['normalization/'+f for f in ('estimation_counts.tsv','r_size_factors.tsv','r_diagnostics.tsv','r_diagnostics.tsv.session.txt')]
                if method=='controls': files+=['normalization/control_indices.txt']
            result[space+'/'+method]=dict(args=args, external=inputs, depth=depth_path if auto_depth else None,
                bams=bams, names=names, final=final, files=[final+'/'+f for f in files],
                mem_mb=integer(block.get('mem_mb',8192),'bigwigs.mem_mb',1))
    return result


def run(plan, repo):
    """Preserve every attempt; publish an entirely completed directory by rename."""
    final=Path(plan['final']).resolve()
    attempt=final.parents[2]/'.bigwig_attempts'/final.parent.name/final.name/uuid.uuid4().hex
    attempt.mkdir(parents=True)
    target=attempt/'result'
    command=[sys.executable,str(repo/'workflow/scripts/bigwig2.py'),*plan['args'],'--output-dir',str(target)]
    record={'command':command,'status':'running','attempt':str(attempt)}
    try:
        # Snakemake removes declared outputs after a failed job. Archive a prior
        # publication before starting so that cleanup cannot destroy its files.
        if final.exists(): final.rename(attempt/'previous_publication')
        with (attempt/'driver.log').open('w') as log:
            env=dict(os.environ, **{k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')})
            completed=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env=env)
        record['exit_status']=completed.returncode
        if completed.returncode: raise ValueError('BigWig failed; preserved diagnostics: '+str(attempt/'driver.log'))
        meta=json.loads((target/'provenance.json').read_text())
        if meta['status']!='complete': raise ValueError('Incomplete BigWig provenance')
        record['status']='complete'
        (target/'workflow.json').write_text(json.dumps(record,indent=2)+'\n')
        for filename in plan['files']:
            if not (target/Path(filename).resolve().relative_to(final)).is_file(): raise ValueError('Missing completed output '+filename)
        final.parent.mkdir(parents=True,exist_ok=True)
        target.rename(final)
    except BaseException as exc:
        record.update(status='failed',error=str(exc)); raise
    finally:
        (attempt/'attempt.json').write_text(json.dumps(record,indent=2)+'\n')
