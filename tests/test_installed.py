"""Fresh-install acceptance; set the two installed command paths to exercise it."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

BASE=os.environ.get('CAPZYME_INSTALLED_BASE')
FULL=os.environ.get('CAPZYME_INSTALLED_FULL')

@unittest.skipUnless(BASE and FULL,'Set installed base+tracks and workflow command paths')
class InstalledTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix='capzyme-installed-');self.addCleanup(self.tmp.cleanup)
        self.p=Path(self.tmp.name).resolve()
        self.env={k:v for k,v in os.environ.items() if k!='PYTHONPATH'}
        self.env['PYTHONDONTWRITEBYTECODE']='1'
        # Put the selected environment before pre-existing tools/Python commands.
        self.env['PATH']=str(Path(FULL).parent)+os.pathsep+self.env['PATH']

    def call(self,cli,*args,ok=True):
        env=dict(self.env)
        if cli==BASE or cli==str(Path(BASE).with_name('python')):
            env['PATH']=str(Path(BASE).parent)+':/opt/homebrew/bin:/usr/bin:/bin'
        r=subprocess.run([cli,*map(str,args)],cwd=self.p,env=env,text=True,capture_output=True)
        self.assertEqual(r.returncode==0,ok,r.stdout+r.stderr)
        return r

    def test_standalone_commands_without_snakemake_and_matrix_bam_agreement(self):
        py=str(Path(BASE).with_name('python'))
        self.call(py,'-c',"import importlib.util; assert importlib.util.find_spec('snakemake') is None; assert importlib.util.find_spec('yaml') is None; import shutil; assert shutil.which('snakemake') is None")
        self.call(BASE,'--help')
        for tool in ('count','bigwig','trim','rna-align','genome-align','collapse','index-record'):
            self.call(BASE,tool,'--help')
        reference=self.p/'ref.fa';reference.write_text('>rna\n'+'A'*100+'\n')
        sam=self.p/'reads.sam';sam.write_text('@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:rna\tLN:100\n'
            'forward\t0\trna\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n'
            'reverse\t16\trna\t91\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n')
        bam=self.p/'reads.bam'
        subprocess.run(['samtools','view','-b','-o',str(bam),str(sam)],check=True,env=self.env)
        before=bam.read_bytes()
        self.call(BASE,'count','--bam',bam,'--sample-name','A','--reference',reference,'--output-dir',self.p/'counts')
        for method in ('raw','rpm'):
            options=['--reference',reference,'--cohort','A','--normalization',method,'--rscript','nonexistent-R']
            self.call(BASE,'bigwig','--matrix',self.p/'counts/counts.tsv','--assert-unthresholded',
                *(['--matrix-complete-both-strands'] if method=='rpm' else []),*options,'--output-dir',self.p/(method+'-matrix'))
            self.call(BASE,'bigwig','--bam',bam,'--sample-name','A',*options,'--output-dir',self.p/(method+'-bam'))
            for strand in ('plus','minus'):
                self.assertEqual((self.p/(method+'-matrix')/'tracks'/('A.'+strand+'.bw')).read_bytes(),
                                 (self.p/(method+'-bam')/'tracks'/('A.'+strand+'.bw')).read_bytes())
        self.assertEqual(before,bam.read_bytes());self.assertFalse(Path(str(bam)+'.bai').exists())
        self.call(BASE,'collapse','--reference',reference,'--output-dir',self.p/'collapsed')
        self.assertTrue((self.p/'collapsed/collapsed.fa').exists())
        # Exercise every installed operation, not only its help parser.
        import random
        rng=random.Random(936)
        seq=''.join(rng.choice('ACGT') for _ in range(160))
        mapping_ref=self.p/'mapping.fa';mapping_ref.write_text('>unique\n'+seq+'\n')
        adapter='AGATCGGAAGAGCACACGTCT'
        reads=self.p/'raw.fq';read=seq[30:75]+adapter
        reads.write_text('@test\n'+read+'\n+\n'+'I'*len(read)+'\n')
        self.call(BASE,'trim','--input-file',reads,'--adapter',adapter,'--output-dir',self.p/'trim')
        for kind in ('bowtie2','hisat2'):
            prefix=self.p/(kind+'-idx')
            subprocess.run([kind+'-build',str(mapping_ref),str(prefix)],check=True,env=self.env,capture_output=True)
            self.call(BASE,'index-record','--reference',mapping_ref,'--index',prefix,'--aligner',kind,'--output-dir',self.p/(kind+'-record'))
            if kind=='bowtie2':
                self.call(BASE,'rna-align','--input-file',self.p/'trim/trimmed.fastq.gz','--reference',mapping_ref,
                    '--index',prefix,'--output-dir',self.p/'rna')
            self.call(BASE,'genome-align','--aligner',kind,'--no-decoys','--input-file',self.p/'trim/trimmed.fastq.gz',
                '--reference',mapping_ref,'--index',prefix,'--output-dir',self.p/(kind+'-genome'))
        for directory in ('rna','bowtie2-genome','hisat2-genome'):
            records=subprocess.check_output(['samtools','view','-c',str(self.p/directory/'accepted.bam')],text=True,env=self.env)
            self.assertEqual(records.strip(),'1')


    def test_readme_workflow_example_and_normal_restart(self):
        demo=self.p/'example'
        self.call(FULL,'workflow','example','--directory',demo,'--with-bigwigs')
        config=demo/'config.yaml'
        self.call(FULL,'workflow','doctor')
        for op in ('check','preview','run'):
            self.call(FULL,'workflow',op,'--config',config)
        self.assertEqual((demo/'run/counts/5ends.raw.tsv').read_bytes(),(demo/'expected.tsv').read_bytes())
        tracks=list((demo/'run/bigwigs').rglob('*.bw'));self.assertEqual(len(tracks),4)
        before={p:(p.read_bytes(),p.stat().st_mtime_ns) for p in tracks}
        result=self.call(FULL,'workflow','run','--config',config)
        self.assertIn('Nothing to be done',result.stdout+result.stderr)
        for p,v in before.items():self.assertEqual((p.read_bytes(),p.stat().st_mtime_ns),v)
        # Native-tool provenance in a wheel identifies shipped code by checksum.
        package=Path(FULL).parent.parent/'lib/python3.12/site-packages/capzyme_resources'
        snapshot=json.loads((package/'source_snapshot.json').read_text())
        self.assertEqual(snapshot['kind'],'installed-wheel')
        import hashlib
        for rel,h in snapshot['files_sha256'].items():self.assertEqual(hashlib.sha256((package/rel).read_bytes()).hexdigest(),h)
