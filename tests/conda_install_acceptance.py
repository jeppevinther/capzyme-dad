#!/usr/bin/env python3
"""Installed-only acceptance for the complete native Conda convenience environment.

No source-tree scientific imports. Reuses the installed-wheel checker, adding
strict executable/package origin checks and tiny Bowtie2/HISAT2 mapping fixtures.
R is intentionally outside this environment and is never invoked.
"""
import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import random
import shutil
import sys
import time
import traceback

from linux_installed_acceptance import Acceptance, check, sha

TOOLS=('python','capzyme','cutadapt','bowtie2','bowtie2-build','bowtie2-inspect',
       'bowtie2-align-s','bowtie2-build-s','bowtie2-inspect-s','hisat2','hisat2-build',
       'hisat2-inspect','hisat2_extract_splice_sites.py','hisat2_extract_exons.py','hisat2-align-s','hisat2-build-s','hisat2-inspect-s','samtools','bedtools','snakemake')


class CondaAcceptance(Acceptance):
    def installed(self,python,standalone=False):
        # The recommended convenience environment intentionally contains Snakemake.
        return super().installed(python,standalone=False)

    def origins(self):
        prefix=Path(sys.prefix).resolve();records=[json.loads(p.read_text()) for p in (prefix/'conda-meta').glob('*.json')]
        expected=self.a.platform
        check(all(r.get('subdir') in (expected,'noarch') for r in records),'Foreign architecture package found')
        check(all(r['name'] not in ('r-base','bioconductor-deseq2') for r in records),'R was unexpectedly installed')
        owners={f:r['name'] for r in records for f in r.get('files',[])}
        result={}
        for tool in TOOLS:
            path=shutil.which(tool,path=str(prefix/'bin')+':/usr/bin:/bin')
            check(path is not None,'Missing required executable '+tool)
            resolved=Path(path).resolve();check(resolved.is_relative_to(prefix),'Executable outside new Conda environment: '+tool)
            rel=str(Path(path).relative_to(prefix))
            owner='pip:capzyme-dad' if tool=='capzyme' else owners.get(rel)
            check(owner is not None,'Executable is not owned by a resolved environment package: '+tool)
            result[tool]={'path':path,'resolved_path':str(resolved),'owner':owner,'sha256':sha(resolved)}
        self.versions['executable_origins']=result
        self.versions['conda_packages']=[{k:r.get(k) for k in ('name','version','build','channel','subdir')} for r in records]
        self.results['native_environment_origins']='passed; required executables and aligner helpers belong to this prefix'
        # Source installation must be a regular wheel install, never editable.
        dist=importlib.metadata.distribution('capzyme-dad')
        direct=json.loads(dist.read_text('direct_url.json'))
        check(not direct.get('dir_info',{}).get('editable',False),'Editable installation is not allowed')

    def standalone_check(self):
        super().standalone_check()
        value=self.results.pop('standalone_without_snakemake')
        self.results['standalone_in_convenience_environment']=dict(value,snakemake_present=True,
            note='Standalone operations tested directly; package dependency separation is unchanged and independently tested previously')
        factors=self.root/'standalone-factors.tsv';factors.write_text('sample_id\tsize_factor\nA\t2\n')
        out=self.root/'standalone/supplied'
        self.call([self.full,'bigwig','--matrix',self.root/'standalone/counts/counts.tsv','--assert-unthresholded',
            '--reference',self.root/'standalone/ref.fa','--cohort','A','--normalization','supplied',
            '--size-factors',factors,'--factor-kind','size-factor','--rscript','deliberately-unavailable-Rscript','--output-dir',out])
        for strand,pos in [('plus',1),('minus',100)]:self.bw(out/'tracks'/f'A.{strand}.bw',[('rna',pos,.5)],{'rna':100})
        self.results['supplied_without_R']='passed'

    def mapping(self):
        p=self.root/'mapping';p.mkdir();rng=random.Random(253)
        chrom=list(''.join(rng.choice('ACGT') for _ in range(1200)))
        chrom[220:222]='GT';chrom[418:420]='AG';chrom=''.join(chrom)
        ref=p/'genome.fa';ref.write_text('>chr\n'+chrom+'\n')
        gtf=p/'genes.gtf';gtf.write_text(''.join(f'chr\tsynthetic\texon\t{a}\t{b}\t.\t+\t.\tgene_id "g"; transcript_id "t";\n' for a,b in [(101,220),(421,540)]))
        rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
        junction=chrom[170:220]+chrom[420:470]
        for kind in ('bowtie2','hisat2'):
            prefix=p/(kind+'-index')
            args=[kind+'-build']
            if kind=='hisat2':
                ss=p/'splice.tsv';exons=p/'exons.tsv'
                ss.write_text(self.call(['hisat2_extract_splice_sites.py',gtf]))
                exons.write_text(self.call(['hisat2_extract_exons.py',gtf]))
                check(ss.stat().st_size>0 and exons.stat().st_size>0,'Empty native annotation tables')
                args+=['--ss',ss,'--exon',exons]
            self.call([*args,ref,prefix])
            self.call([self.full,'index-record','--reference',ref,'--index',prefix,'--aligner',kind,'--output-dir',p/(kind+'-record')])
            reads=[('plus',chrom[600:660]),('minus',rc(chrom[700:760])),('unmapped','N'*60)]
            if kind=='hisat2':reads += [('splice_plus',junction),('splice_minus',rc(junction))]
            fastq=p/(kind+'.fq');fastq.write_text(''.join('@'+n+'\n'+s+'\n+\n'+'I'*len(s)+'\n' for n,s in reads))
            source={str(x):(sha(x),x.stat().st_mtime_ns) for x in [ref,fastq,*p.glob(kind+'-index.*')]}
            out=p/(kind+'-map')
            self.call([self.full,'genome-align','--aligner',kind,'--input-file',fastq,'--reference',ref,
                '--index',prefix,'--no-decoys','--threads','2','--output-dir',out])
            self.call([self.full,'count','--bam',out/'accepted.bam','--sample-name','test','--reference',ref,'--output-dir',p/(kind+'-counts')])
            import csv
            with (p/(kind+'-counts')/'counts.tsv').open() as f:
                observed={(r['chromosome'],r['strand'],int(r['start'])):int(r['test']) for r in csv.DictReader(f,delimiter='\t')}
            wanted={('chr','+',601):1,('chr','-',760):1}
            if kind=='hisat2':wanted.update({('chr','+',171):1,('chr','-',470):1})
            check(observed==wanted,'Incorrect mapping/5-prime coordinates for '+kind)
            sam=self.call(['samtools','view',out/'accepted.bam'])
            if kind=='hisat2':
                splice={f[0]:f[5] for line in sam.splitlines() if (f:=line.split('\t'))[0].startswith('splice')}
                check(splice=={'splice_plus':'50M200N50M','splice_minus':'50M200N50M'},'Spliced CIGAR mismatch')
            check((out/'unmapped.fastq').read_text().splitlines()[0]=='@unmapped','Unmapped handoff differs')
            for path,v in source.items():check((sha(path),Path(path).stat().st_mtime_ns)==v,'Synthetic mapping input/index changed')
            self.results[kind+'_synthetic_mapping']={'status':'passed','accepted':len(wanted),'unmapped':1,'exact_coordinates':True}

    def run(self):
        self.origins()
        for kind in ('bowtie2','hisat2'):
            self.call([self.full,'workflow','doctor','--route','genome','--aligner',kind])
        super().run()
        # Optional R absence is intentional in this convenience environment.
        self.gaps=[]
        self.results['optional_R']='Not installed or tested; raw/RPM/supplied passed with unavailable Rscript'
        self.mapping()
        self.origins()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',required=True)
    parser.add_argument('--platform',required=True,choices=['linux-64','osx-arm64'])
    parser.add_argument('--manifest-sha256',required=True,help='Reviewed manifest identity for the installed artifact')
    a=parser.parse_args()
    a.workflow_python=sys.executable;a.standalone_python=sys.executable
    a.native_path='/usr/bin:/bin'
    a.rscript=None;a.r_gap='R intentionally excluded from the recommended Conda environment'
    runner=CondaAcceptance(a);status='failed';error=None
    try:runner.run();status='passed'
    except BaseException:error=traceback.format_exc();print(error,file=sys.stderr);raise
    finally:
        result={'status':status,'elapsed_seconds':time.time()-runner.started,'error':error,'checks':runner.results,
                'versions':runner.versions,'arguments':vars(a),'coverage_gaps':runner.gaps}
        (runner.root/'acceptance.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps({'status':status,'receipt':str(runner.root/'acceptance.json')}),flush=True)

if __name__=='__main__':main()
