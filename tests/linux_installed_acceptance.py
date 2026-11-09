#!/usr/bin/env python3
"""Synthetic installed-wheel acceptance, independent of a source checkout.

Run with the workflow environment's Python. Supply a separate tracks-only Python
and a Linux-native executable PATH. All inputs/outputs are generated under a new
external directory; commands, statuses and fixtures are retained. No research
input arguments are accepted. This checks existing scientific implementations.
"""
import argparse
import csv
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import time
import traceback


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()


def state(root):
    return {str(p.relative_to(root)):{'sha256':sha(p),'bytes':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns}
            for p in sorted(root.rglob('*')) if p.is_file()}


def check(condition, message):
    if not condition: raise AssertionError(message)


def f32(value):
    return struct.unpack('f',struct.pack('f',value))[0]


class Acceptance:
    def __init__(self,a):
        self.a=a;self.root=Path(a.directory).resolve()
        check(not self.root.exists(),'Acceptance directory must be new')
        check(not any((p/'.git').exists() for p in (self.root,*self.root.parents)), 'Run outside source checkouts')
        self.root.mkdir(parents=True);self.audit=self.root/'audit';self.audit.mkdir()
        self.env={k:v for k,v in os.environ.items() if k not in ('PYTHONPATH','PYTHONHOME')}
        self.env.update(PYTHONDONTWRITEBYTECODE='1',PYTHONNOUSERSITE='1',
            **{k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')})
        self.full=str(Path(a.workflow_python).absolute().with_name('capzyme'))
        self.standalone=str(Path(a.standalone_python).absolute().with_name('capzyme'))
        self.native=a.native_path;self.number=0;self.results={};self.started=time.time()
        self.versions={};self.gaps=[]

    def call(self,args,standalone=False):
        self.number+=1;args=list(map(str,args))
        env=dict(self.env);prefix=Path(self.a.standalone_python if standalone else self.a.workflow_python).absolute().parent
        env['PATH']=str(prefix)+os.pathsep+self.native
        record={'command':args,'cwd':str(self.root),'PATH':env['PATH'],'started_unix':time.time(),
                'log':f'command_{self.number:03d}.log'}
        try:
            r=subprocess.run(args,cwd=self.root,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,errors='replace')
            (self.audit/record['log']).write_text(r.stdout);record['exit_status']=r.returncode
        except BaseException as exc:
            record.update(exit_status=None,error=str(exc));raise
        finally:
            record['elapsed_seconds']=time.time()-record['started_unix']
            with (self.audit/'commands.jsonl').open('a') as f:f.write(json.dumps(record)+'\n')
        check(r.returncode==0,f'Command failed ({r.returncode}): {args}; see {record["log"]}')
        return r.stdout

    def workflow(self,op,cfg,route='rna'):
        return self.call([self.full,'workflow',op,'--route',route,'--config',cfg,
            *(['--cores','4','--memory-mb','16384'] if op!='check' else [])])

    def bw(self,path,expected,dictionary):
        import pyBigWig
        with pyBigWig.open(str(path)) as bw:
            check(list(bw.chroms().items())==list(dictionary.items()),'BigWig reference dictionary/order mismatch')
            for chrom in dictionary:
                actual=list(bw.intervals(chrom) or [])
                wanted=[(p-1,p,f32(v)) for c,p,v in expected if c==chrom]
                check(actual==wanted,f'BigWig interval/value mismatch: {path.name}, {chrom}: {actual} != {wanted}')

    def tracks(self,directory,rows,factors,samples=('B','A'),strands=('plus',),dictionary=None):
        # The generated fixture combines viral records before collapsed ncRNAs.
        dictionary=dictionary or {'virus_plus':140,'virus_minus':140,'RNU2__n2':160}
        expected_names={n+'.'+sign+'.bw' for n in samples for sign in strands}
        check({p.name for p in (directory/'tracks').iterdir()}==expected_names,'Unexpected/missing tracks')
        for n in samples:
            for sign in strands:
                entries=[(chrom,pos,count[n]*factors[n]) for chrom,strand,pos,count in rows
                         if strand==('+' if sign=='plus' else '-') and count[n]]
                self.bw(directory/'tracks'/f'{n}.{sign}.bw',entries,dictionary)
        meta=json.loads((directory/'provenance.json').read_text())
        check(meta['status']=='complete','Incomplete provenance')
        for name,h in meta['outputs_sha256'].items():check(sha(directory/name)==h,'Published artifact checksum mismatch')
        return meta

    def installed(self,python,standalone=False):
        code='''import hashlib,importlib.util,importlib.metadata as m,json,pathlib,shutil,sys
s=importlib.util.find_spec('capzyme_resources');root=pathlib.Path(next(iter(s.submodule_search_locations)))
snapshot=json.loads((root/'source_snapshot.json').read_text())
assert all(hashlib.sha256((root/p).read_bytes()).hexdigest()==h for p,h in snapshot['files_sha256'].items())
print(json.dumps({'python':sys.version,'executable':sys.executable,'package_root':str(root),'manifest_sha256':hashlib.sha256((root/'source_snapshot.json').read_bytes()).hexdigest(),'runtime_files':len(snapshot['files_sha256']),'packages':{n:m.version(n) for n in ['capzyme-dad','pyBigWig','numpy']},'snakemake_module':importlib.util.find_spec('snakemake') is not None,'snakemake_executable':shutil.which('snakemake')}))
'''
        value=json.loads(self.call([python,'-c',code],standalone))
        if standalone:check(not value['snakemake_module'] and value['snakemake_executable'] is None,'Standalone setup exposes Snakemake')
        return value

    def run(self):
        self.versions['workflow']=self.installed(self.a.workflow_python)
        self.versions['standalone']=self.installed(self.a.standalone_python,True)
        check(self.versions['workflow']['manifest_sha256']==self.a.manifest_sha256==self.versions['standalone']['manifest_sha256'],'Installed runtime manifest differs from reviewed wheel')
        for command in (self.full,self.standalone):
            self.call([command,'--help'],command==self.standalone)
            for tool in ('count','bigwig','trim','rna-align','genome-align','collapse','index-record'):
                self.call([command,tool,'--help'],command==self.standalone)
        self.results['installed_help_and_manifest']='passed'
        self.versions['doctor']=self.call([self.full,'workflow','doctor'])
        self.results['dependency_doctor']='passed'
        demo=self.root/'example';self.call([self.full,'workflow','example','--directory',demo,'--with-bigwigs'])
        cfg=demo/'config.yaml';data=json.loads(cfg.read_text());out=demo/'run'
        source_inputs={str(p):{'sha256':sha(p),'mtime_ns':p.stat().st_mtime_ns} for p in demo.iterdir() if p.is_file() and p!=cfg}
        self.workflow('check',cfg);self.workflow('preview',cfg)
        check(not out.exists(),'Preview created outputs')
        self.workflow('run',cfg)
        expected=('chromosome\tstrand\tstart\tend\tB\tA\tempty\nRNU2__n2\t+\t1\t1\t12\t0\t0\n'
                  'RNU2__n2\t+\t11\t11\t0\t3\t0\nvirus_minus\t+\t21\t21\t0\t1\t0\nvirus_plus\t+\t10\t10\t0\t1\t0\n')
        check((out/'counts/5ends.raw.tsv').read_text()==(demo/'expected.tsv').read_text()==expected,'Count matrix/expected fixture mismatch')
        self.results['exact_fixture_count_matrix']={'status':'passed','sample_order':['B','A','empty'],'features':4,'totals':[12,5,0],'sha256':sha(out/'counts/5ends.raw.tsv')}
        rows=[('RNU2__n2','+',1,{'B':12,'A':0}),('RNU2__n2','+',11,{'B':0,'A':3}),
              ('virus_minus','+',21,{'B':0,'A':1}),('virus_plus','+',10,{'B':0,'A':1})]
        depths={n:int(self.call(['samtools','view','-c',out/'bam'/f'{n}.accepted.bam'])) for n in ('B','A','empty')}
        check(depths=={'B':12,'A':5,'empty':0},'Independent accepted depths differ')
        self.tracks(out/'bigwigs/rna/raw',rows,{'B':1,'A':1})
        self.tracks(out/'bigwigs/rna/rpm',rows,{n:1e6/depths[n] for n in ('B','A')})
        for method in ('raw','rpm'):
            with (out/f'bigwigs/rna/{method}/normalization_factors.tsv').open() as f:factor_rows=list(csv.DictReader(f,delimiter='\t'))
            check([r['sample_id'] for r in factor_rows]==['B','A'],'Factor cohort/order differs')
            for row in factor_rows:
                n=row['sample_id'];check(float(row['multiplier'])==(1 if method=='raw' else 1e6/depths[n]),'Factor calculation differs')
                if method=='rpm':check(int(row['accepted_records'])==depths[n],'RPM depth differs')
        self.results['four_raw_rpm_tracks']={'status':'passed','tracks':4,'nonzero_values':8,'independent_depths':depths,'comparison':'exact float32 intervals/values; no extra intervals'}
        before=state(out);text=self.workflow('run',cfg)
        check('Nothing to be done' in text and state(out)==before,'Completed restart changed outputs/scheduled work')
        self.results['completed_noop_restart']={'status':'passed','unchanged_files':len(before)}
        upstream={k:v for k,v in before.items() if not k.startswith(('bigwigs/','.bigwig_attempts/','logs/bigwig.'))}
        old_tracks={k:v for k,v in before.items() if k.startswith('bigwigs/rna/')}
        factors=self.root/'supplied.tsv';factors.write_text('sample_id\tsize_factor\nB\t2\nA\t4\n')
        data['bigwigs']['rna']['methods']['supplied']={'size_factors':str(factors),'factor_kind':'size-factor'}
        cfg.write_text(json.dumps(data,indent=2)+'\n')
        dry=self.workflow('preview',cfg)
        for name in ('trim_fastq','map_filter_sort','count_5ends','threshold_counts'):
            check('rule '+name+':' not in dry,'Adding tracks scheduled upstream work')
        self.workflow('run',cfg);now=state(out)
        check(all(now[k]==v for k,v in upstream.items()),'Adding a method changed upstream files')
        check(all(now[k]==v for k,v in old_tracks.items()),'Adding a method changed existing tracks')
        self.tracks(out/'bigwigs/rna/supplied',rows,{'B':.5,'A':.25})
        self.results['additional_method_without_upstream_work']={'status':'passed','upstream_files_unchanged':len(upstream)}
        preserved={k:v for k,v in now.items() if not k.startswith(('bigwigs/rna/raw/','.bigwig_attempts/rna/raw/'))}
        missing=out/'bigwigs/rna/raw/tracks/A.plus.bw';original_hash=sha(missing);missing.unlink()
        dry=self.workflow('preview',cfg)
        check('bigwig_rna_raw' in dry,'Missing output did not schedule method')
        for name in ('trim_fastq','map_filter_sort','count_5ends','bigwig_rna_rpm','bigwig_rna_supplied'):
            check('rule '+name+':' not in dry,'Missing track scheduled unrelated jobs')
        self.workflow('run',cfg);after=state(out)
        check(sha(missing)==original_hash,'Regenerated track differs')
        check(all(after[k]==v for k,v in preserved.items()),'Recovery changed unaffected files')
        check(list((out/'.bigwig_attempts/rna/raw').glob('*/previous_publication')),'Previous publication not retained')
        self.tracks(out/'bigwigs/rna/raw',rows,{'B':1,'A':1})
        self.results['missing_track_recovery']='passed; affected method regenerated, previous publication retained'
        for path,s in source_inputs.items():check(sha(path)==s['sha256'] and Path(path).stat().st_mtime_ns==s['mtime_ns'],'Synthetic input changed')
        self.standalone_check()
        if self.a.rscript:self.r_check()
        else:self.gaps.append('R/DESeq2 integration not run: '+self.a.r_gap)
        self.results['final_installed_integrity']={'workflow':self.installed(self.a.workflow_python),'standalone':self.installed(self.a.standalone_python,True)}
        self.results['outputs']=state(out)

    def standalone_check(self):
        p=self.root/'standalone';p.mkdir();ref=p/'ref.fa';ref.write_text('>rna\n'+'A'*100+'\n')
        sam=p/'reads.sam';sam.write_text('@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:rna\tLN:100\n'
            'forward\t0\trna\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n'
            'reverse\t16\trna\t91\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n')
        bam=p/'reads.bam';self.call(['samtools','view','-b','-o',bam,sam],True)
        before={str(x):(sha(x),x.stat().st_mtime_ns) for x in (ref,bam)}
        self.call([self.standalone,'count','--bam',bam,'--sample-name','A','--reference',ref,'--output-dir',p/'counts'],True)
        check((p/'counts/counts.tsv').read_text()=='chromosome\tstrand\tstart\tend\tA\nrna\t+\t1\t1\t1\nrna\t-\t100\t100\t1\n','Standalone counts differ')
        for method,mult in [('raw',1),('rpm',500000)]:
            common=['--reference',ref,'--cohort','A','--normalization',method,'--rscript','deliberately-unavailable-Rscript']
            self.call([self.standalone,'bigwig','--matrix',p/'counts/counts.tsv','--assert-unthresholded',
                *(['--matrix-complete-both-strands'] if method=='rpm' else []),*common,'--output-dir',p/(method+'-matrix')],True)
            self.call([self.standalone,'bigwig','--bam',bam,'--sample-name','A',*common,'--output-dir',p/(method+'-bam')],True)
            for strand,pos in [('plus',1),('minus',100)]:
                a=p/(method+'-matrix')/'tracks'/f'A.{strand}.bw';b=p/(method+'-bam')/'tracks'/a.name
                check(a.read_bytes()==b.read_bytes(),'Matrix/BAM BigWigs differ')
                self.bw(a,[('rna',pos,mult)],{'rna':100})
        for path,v in before.items():check((sha(path),Path(path).stat().st_mtime_ns)==v,'Standalone input changed')
        check(not Path(str(bam)+'.bai').exists() and not Path(str(ref)+'.fai').exists(),'Input sidecar created')
        self.results['standalone_without_snakemake']={'status':'passed','count_features':2,'matrix_bam_tracks_compared':4,'decoded_tracks':4,'inputs_unchanged':True}

    def r_check(self):
        p=self.root/'r-integration';p.mkdir();ref=p/'ref.fa';ref.write_text('>rna\n'+'A'*100+'\n')
        bams=[]
        for name,scale in [('B',1),('A',2)]:
            sam=p/(name+'.sam');text='@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:rna\tLN:100\n'
            for sign,flag,pos,num in [('+',0,1,5*scale),('-',16,91,10*scale)]:
                for i in range(num):text+=f'{name}_{flag}_{i}\t{flag}\trna\t{pos}\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n'
            sam.write_text(text);bam=p/(name+'.bam');self.call(['samtools','view','-b','-o',bam,sam]);bams.append(bam)
        before={str(x):(sha(x),x.stat().st_mtime_ns) for x in [ref,*bams]}
        sheet=p/'samples.tsv';sheet.write_text('sample_id\tbam\n'+''.join(f'{n}\t{b}\n' for n,b in zip(['B','A'],bams)))
        controls=p/'controls.tsv';controls.write_text('chromosome\tstrand\tstart\tend\nrna\t+\t1\t1\nrna\t-\t100\t100\n')
        cfg=p/'config.json';out=p/'run'
        methods={'deseq2':{'deseq2_type':'ratio','rscript':self.a.rscript},
                 'controls':{'deseq2_type':'poscounts','rscript':self.a.rscript,'controls':str(controls),'min_control_features':2}}
        cfg.write_text(json.dumps({'schema_version':1,'input_mode':'existing_bam','samples':str(sheet),'output_dir':str(out),
            'reference':{'fasta':str(ref),'fasta_sha256':sha(ref)},'counting':{'strand':'both','min_total_count':100},
            'bigwigs':{'rna':{'cohort':['B','A'],'emit_samples':['A'],'strands':'both','methods':methods}}},indent=2)+'\n')
        for op in ('check','preview','run'):self.workflow(op,cfg,'bam')
        check((out/'counts/5ends.raw.tsv').read_text()=='chromosome\tstrand\tstart\tend\tB\tA\nrna\t+\t1\t1\t5\t10\nrna\t-\t100\t100\t10\t20\n','R fixture raw counts differ')
        independent=p/'independent.R';independent.write_text('''a <- commandArgs(TRUE)
x <- read.delim(a[1], check.names=FALSE, quote="", comment.char="")
m <- as.matrix(x[,c("B","A")]); storage.mode(m) <- "integer"
d <- DESeq2::DESeqDataSetFromMatrix(m, data.frame(row.names=c("B","A"), cohort=c("x","x")), ~1)
if (a[2]=="controls") d <- DESeq2::estimateSizeFactors(d,type="poscounts",controlGenes=seq_len(nrow(m))) else d <- DESeq2::estimateSizeFactors(d,type="ratio")
options(digits=17)
write.table(data.frame(sample_id=names(DESeq2::sizeFactors(d)),size_factor=unname(DESeq2::sizeFactors(d))),a[3],sep="\\t",quote=FALSE,row.names=FALSE)
cat(R.version.string,"\\nDESeq2",as.character(packageVersion("DESeq2")),"\\n")
''')
        factors={};results={}
        rows=[('rna','+',1,{'B':5,'A':10}),('rna','-',100,{'B':10,'A':20})]
        for method in methods:
            native=p/(method+'.native.tsv')
            version=self.call([self.a.rscript,'--vanilla',independent,out/'counts/5ends.raw.tsv',method,native])
            self.versions['R_DESeq2']=version
            with native.open() as f:factors={r['sample_id']:float(r['size_factor']) for r in csv.DictReader(f,delimiter='\t')}
            folder=out/'bigwigs/rna'/method
            with (folder/'normalization_factors.tsv').open() as f:got=list(csv.DictReader(f,delimiter='\t'))
            check([r['sample_id'] for r in got]==['B','A'],'R factor cohort differs')
            for r in got:check(abs(float(r['size_factor'])/factors[r['sample_id']]-1)<=1e-12,'Independent R factors differ')
            meta=self.tracks(folder,rows,{n:1/factors[n] for n in factors},samples=('A',),strands=('plus','minus'),dictionary={'rna':100})
            check(meta['cohort']==['B','A'] and meta['emitted_samples']==['A'],'R cohort/emitted subset conflated')
            args=[self.full,'bigwig','--matrix',out/'counts/5ends.raw.tsv','--assert-unthresholded','--reference',ref,
                  '--cohort','B','A','--emit-samples','A','--strand','both','--normalization',method,
                  '--deseq2-type',methods[method]['deseq2_type'],'--rscript',self.a.rscript,'--output-dir',p/('direct-'+method)]
            if method=='controls':args+=['--controls',controls,'--min-control-features','2']
            self.call(args)
            for track in (folder/'tracks').iterdir():check(track.read_bytes()==(p/('direct-'+method)/'tracks'/track.name).read_bytes(),'Standalone/workflow R tracks differ')
            results[method]={'estimator':methods[method]['deseq2_type'],'native_factors':factors,'tracks':2,'independent_R_factor_and_decoded_track_agreement':'passed'}
        for path,v in before.items():check((sha(path),Path(path).stat().st_mtime_ns)==v,'R fixture input changed')
        self.results['R_integration']=results


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',required=True)
    parser.add_argument('--workflow-python',required=True)
    parser.add_argument('--standalone-python',required=True)
    parser.add_argument('--native-path',required=True,help='Linux executable PATH; exclude existing Snakemake directories')
    parser.add_argument('--manifest-sha256',required=True)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--rscript',help='Existing compatible isolated R/DESeq2 only; never installed by this driver')
    group.add_argument('--r-gap',help='Explicit explanation if no compatible existing R/DESeq2 was available')
    a=parser.parse_args();runner=Acceptance(a);status='failed';error=None
    try:runner.run();status='passed' if not runner.gaps else 'passed_with_coverage_gap'
    except BaseException as exc:error=traceback.format_exc();print(error,file=sys.stderr);raise
    finally:
        result={'status':status,'error':error,'elapsed_seconds':time.time()-runner.started,'checks':runner.results,
                'versions':runner.versions,'coverage_gaps':runner.gaps,'arguments':vars(a)}
        (runner.root/'acceptance.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps({'status':status,'receipt':str(runner.root/'acceptance.json')}),flush=True)


if __name__=='__main__':main()
