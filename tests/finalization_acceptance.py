#!/usr/bin/env python3
"""Compact installed-package acceptance of the approved public workflow templates."""
import argparse
import gzip
import json
import random
from pathlib import Path
import sys
import time
import traceback
import yaml

from conda_install_acceptance import CondaAcceptance
from linux_installed_acceptance import check,sha,state


class FinalizationAcceptance(CondaAcceptance):
    def trim_off(self):
        p=self.root/'trim-modes';p.mkdir();seq='ACCTAACCTAACCTAACCTAA'+'G'*12
        fq=p/'reads.fq';fq.write_text('@high\n'+seq+'\n+\n'+'I'*len(seq)+'\n@low\n'+seq+'\n+\n'+'I'*(len(seq)-4)+'!!!!\n')
        for name,options in [('off',['--no-nextseq-trim']),('on',['--nextseq-trim','20']),('zero',['--nextseq-trim','0'])]:
            self.call([self.full,'trim','--input-file',fq,'--adapter','TTTTTTTTTTTTTTTT','--quality','20','--minimum-length','1','--output-dir',p/name,*options])
        def sequences(name):
            with gzip.open(p/name/'trimmed.fastq.gz','rt') as f:return f.read().splitlines()[1::4]
        check(sequences('off')==[seq,seq[:-4]],'NextSeq off or ordinary quality trimming failed')
        check(sequences('on')[0]==sequences('zero')[0]==seq[:-12],'Numeric NextSeq behavior changed')
        cfg=p/'settings.json';cfg.write_text(json.dumps({'schema_version':2,'input':{'type':'raw_fastq'},'trimming':{'adapter':'TTTTTTTTTTTTTTTT','minimum_length':1,'nextseq_trim':None}}))
        self.call([self.full,'trim','--settings',cfg,'--input-file',fq,'--output-dir',p/'null'])
        check(sequences('null')==sequences('off'),'Explicit null not preserved in settings')
        self.results['NextSeq_off']='CLI/null omit extra treatment; quality trimming active; 0/20 retain established behavior'

    def examples(self):
        p=self.root/'public-examples';p.mkdir();inputs=p/'inputs';inputs.mkdir()
        rng=random.Random(991);rna=''.join(rng.choice('ACGT') for _ in range(48))+'G'*12+''.join(rng.choice('ACGT') for _ in range(60))
        ref=inputs/'rna.fa';ref.write_text('>rna\n'+rna+'\n');idx=inputs/'rna-index'
        self.call(['bowtie2-build',ref,idx])
        genome=self.root/'mapping/genome.fa';chrom=''.join(genome.read_text().splitlines()[1:]);rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
        plus=chrom[600:660];minus=rc(chrom[700:760]);adapter='AGATCGGAAGAGCACACGTCT'
        for name in ['A','B']:
            records=[(name+'_rna',rna[:60]+adapter,'I'*(60+len(adapter))),(name+'_plus',plus+adapter,'I'*(60+len(adapter))),(name+'_minus',minus+adapter,'I'*(60+len(adapter)))]
            if name=='A':records.append(('A_low',rna[:60],'I'*56+'!!!!'))
            (inputs/(name+'.raw.fq')).write_text(''.join('@'+n+'\n'+s+'\n+\n'+q+'\n' for n,s,q in records))
            (inputs/(name+'.trimmed.fq')).write_text(''.join('@'+name+'_'+n+'\n'+s+'\n+\n'+'I'*len(s)+'\n' for n,s in [('plus',plus),('minus',minus)]))
        before={str(x):(sha(x),x.stat().st_mtime_ns) for x in inputs.iterdir()}
        for route,stem in [('sequential','sequential'),('genome','genome')]:
            work=p/stem;work.mkdir()
            cfg=work/(stem+'.settings.yaml')
            data=yaml.safe_load((Path(self.a.examples)/(stem+'.settings.yaml')).read_text())
            data['genome'].update(reference=str(genome),index_prefix=str(self.root/'mapping'/('hisat2-index' if route=='sequential' else 'bowtie2-index')),decoys=None,no_decoys=True)
            if route=='sequential':
                data['reference'].update(fasta=str(ref),index_prefix=str(idx))
                data['preprocessing']['adapter']=adapter
            sheet=work/(stem+'.samples.tsv')
            # Preserve the approved sample-sheet header and order, replacing only paths.
            template=(Path(self.a.examples)/(stem+'.samples.tsv')).read_text().splitlines()
            check(template[0]=='sample_id\tfastq' and [x.split('\t')[0] for x in template[1:]]==['A','B'],'Public sample sheet changed')
            sheet.write_text('sample_id\tfastq\n'+''.join(n+'\t'+str(inputs/(n+('.raw.fq' if route=='sequential' else '.trimmed.fq')))+'\n' for n in ['A','B']))
            cfg.write_text(yaml.safe_dump(data,sort_keys=False))
            for op in ['check','preview','run']:self.workflow(op,cfg,route)
            out=work/data['output_dir'];matrix=out/'counts/genome.5ends.raw.tsv'
            check(matrix.read_text()=='chromosome\tstrand\tstart\tend\tA\tB\nchr\t+\t601\t601\t1\t1\nchr\t-\t760\t760\t1\t1\n','Public genome matrix mismatch')
            if route=='sequential':
                check((out/'counts/5ends.raw.tsv').read_text()=='chromosome\tstrand\tstart\tend\tA\tB\nrna\t+\t1\t1\t2\t1\n','Public RNA matrix mismatch')
                observed=json.loads((out/'reference/fixed_verified.json').read_text())
                check(len(observed['inputs_sha256'])==7 and observed['index_matches_fasta'],'Automatic fixed-reference provenance missing')
                with gzip.open(out/'trimmed/A.fastq.gz','rt') as f:lines=f.read().splitlines()
                reads={lines[i][1:]:lines[i+1] for i in range(0,len(lines),4)}
                check(reads['A_rna']==rna[:60] and reads['A_low']==rna[:56],'Workflow NextSeq-off/quality behavior differs')
            spaces=['rna','genome'] if route=='sequential' else ['genome']
            for space in spaces:
                sizes={'rna':120} if space=='rna' else {'chr':1200}
                for method in ['raw','rpm']:
                    for n in ['A','B']:
                        for strand in ['plus','minus']:
                            if space=='rna':
                                count=2 if n=='A' else 1
                                expected=[('rna',1,count if method=='raw' else 1e6)] if strand=='plus' else [('rna',1,0.0)]
                            else:expected=[('chr',601 if strand=='plus' else 760,1 if method=='raw' else 500000)]
                            self.bw(out/'bigwigs'/space/method/'tracks'/f'{n}.{strand}.bw',expected,sizes)
            saved=state(out);text=self.workflow('run',cfg,route)
            check('Nothing to be done' in text and state(out)==saved,'Public example restart changed outputs')
            upstream={k:v for k,v in saved.items() if not k.startswith(('bigwigs/','.bigwig_attempts/','logs/bigwig.'))}
            factors=work/'supplied.tsv';factors.write_text('sample_id\tsize_factor\nA\t2\nB\t4\n')
            data['bigwigs']['genome']['methods']['supplied']={'size_factors':str(factors),'factor_kind':'size-factor'}
            cfg.write_text(yaml.safe_dump(data,sort_keys=False));dry=self.workflow('preview',cfg,route)
            for rule in ['trim_fastq','map_rna_handoff','map_genome','count_genome','count_5ends']:
                check('rule '+rule+':' not in dry,'Additional tracks scheduled upstream jobs')
            self.workflow('run',cfg,route);after=state(out)
            check(all(after[k]==v for k,v in saved.items()),'Additional tracks changed existing files')
            self.results['public_example_'+route]={'status':'passed','sample_order':['A','B'],'upstream_files_preserved':len(upstream),
                'raw_rpm_tracks':16 if route=='sequential' else 8,'restart_unchanged':True,'hash_expectations':'omitted'}
        for path,v in before.items():check((sha(path),Path(path).stat().st_mtime_ns)==v,'Public example synthetic input/index changed')

    def run(self):
        super().run();self.trim_off();self.examples()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',required=True);parser.add_argument('--examples',required=True)
    parser.add_argument('--platform',required=True,choices=['linux-64','osx-arm64'])
    parser.add_argument('--manifest-sha256',required=True)
    a=parser.parse_args();a.workflow_python=sys.executable;a.standalone_python=sys.executable
    a.native_path='/usr/bin:/bin';a.rscript=None;a.r_gap='Optional R environment is unchanged and outside finalization scope'
    r=FinalizationAcceptance(a);status='failed';error=None
    try:r.run();status='passed'
    except BaseException:error=traceback.format_exc();print(error,file=sys.stderr);raise
    finally:
        (r.root/'acceptance.json').write_text(json.dumps({'status':status,'error':error,'elapsed_seconds':time.time()-r.started,
          'checks':r.results,'versions':r.versions,'arguments':vars(a)},indent=2)+'\n')
        print(json.dumps({'status':status,'receipt':str(r.root/'acceptance.json')}),flush=True)

if __name__=='__main__':main()
