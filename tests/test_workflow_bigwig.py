"""Synthetic workflow/standalone equality, dependency isolation and safe retries."""
import copy
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'workflow/scripts'))
from core import sha256
from workflow_bigwig import plans, run
import test_core as fixtures

TOOLS=all(shutil.which(t) for t in ('samtools','bedtools','snakemake')) and importlib.util.find_spec('pyBigWig')
R=os.environ.get('CAPZYME_BIGWIG_RSCRIPT','Rscript')

@unittest.skipUnless(TOOLS,'Counting, workflow and track dependencies required')
class WorkflowBigwigTests(unittest.TestCase):
    def setUp(self):
        self.f=fixtures.CoordinateTests('test_empty_matrix_and_total_threshold');self.f.setUp();self.addCleanup(self.f.doCleanups)
        self.p=self.f.p; self.bam=self.f.bam()
        self.before={p:(sha256(p),p.stat().st_mtime_ns) for p in (self.bam,self.f.fasta)}
        sheet=self.p/'samples.tsv';sheet.write_text(f'sample_id\tbam\nB\t{self.bam}\nA\t{self.bam}\n')
        # Two files are required by the existing input contract.
        other=self.p/'A.bam';shutil.copyfile(self.bam,other)
        sheet.write_text(f'sample_id\tbam\nB\t{self.bam}\nA\t{other}\n')
        self.cfg={'schema_version':1,'input_mode':'existing_bam','samples':str(sheet),'output_dir':str(self.p/'results'),
            'reference':{'fasta':str(self.f.fasta),'fasta_sha256':sha256(self.f.fasta)},
            'counting':{'strand':'both','min_total_count':100}}
        self.path=self.p/'config.json'

    def workflow(self,dry=False,ok=True):
        self.path.write_text(json.dumps(self.cfg))
        r=subprocess.run(['snakemake','--snakefile',str(ROOT/'workflow/existing_bam.smk'),'--configfile',str(self.path),
            '--directory',str(self.p),'--cores','4','--resources','mem_mb=16384',*(['--dry-run'] if dry else [])],capture_output=True,text=True)
        self.assertEqual(r.returncode==0,ok,r.stdout+r.stderr)
        return r.stdout+r.stderr

    def block(self,methods):
        self.cfg['bigwigs']={'rna':{'methods':methods,'cohort':['B','A'],'emit_samples':['A'],'strands':'both'}}

    def test_raw_rpm_addition_no_upstream_rerun_missing_track_and_failure_recovery(self):
        self.workflow()
        root=self.p/'results'
        upstream={p:(sha256(p),p.stat().st_mtime_ns) for p in root.rglob('*') if p.is_file()}
        self.block({'raw':{},'rpm':{}})
        dry=self.workflow(True)
        for name in ('prepare_existing_bam','count_existing_bams','threshold_existing_counts'):
            self.assertNotIn('rule '+name+':',dry)
        self.workflow()
        for p,v in upstream.items(): self.assertEqual((sha256(p),p.stat().st_mtime_ns),v)
        for method in ('raw','rpm'):
            target=root/'bigwigs/rna'/method
            meta=json.loads((target/'provenance.json').read_text())
            self.assertEqual(meta['cohort'],['B','A']);self.assertEqual(meta['emitted_samples'],['A'])
            self.assertEqual(meta['matrix_sha256'],sha256(root/'counts/5ends.raw.tsv'))
            self.assertNotEqual(meta['matrix_sha256'],sha256(root/'counts/5ends.min_total_100.tsv'))
            cmd=[sys.executable,str(ROOT/'toolbox.py'),'bigwig','--matrix',str(root/'counts/5ends.raw.tsv'),
                '--assert-unthresholded','--chrom-sizes',str(root/'reference/sizes.tsv'),'--normalization',method,
                '--cohort','B','A','--emit-samples','A','--output-dir',str(self.p/method)]
            if method=='rpm':cmd+=['--depths',str(root/'bigwigs/rna/accepted_depths.tsv')]
            subprocess.run(cmd,check=True,capture_output=True)
            for track in (target/'tracks').glob('*.bw'):
                self.assertEqual(track.read_bytes(),(self.p/method/'tracks'/track.name).read_bytes())
        self.assertIn('Nothing to be done',self.workflow())
        # Missing one declared output safely republishes a complete directory.
        lost=root/'bigwigs/rna/raw/tracks/A.minus.bw';expected=lost.read_bytes();lost.unlink()
        self.workflow();self.assertEqual(lost.read_bytes(),expected)
        self.assertTrue(list((root/'.bigwig_attempts/rna/raw').glob('*/previous_publication')))
        # Simulate an interrupted/failed first attempt; no final output is published.
        from bam_configuration import load_bam_settings
        settings=load_bam_settings(self.cfg,ROOT)
        plan=plans(self.cfg['bigwigs'],settings,ROOT,bam=True)['rna/raw']
        plan=copy.deepcopy(plan);plan['final']=str(root/'bigwigs/rna/failure')
        plan['files']=[p.replace('/rna/raw/','/rna/failure/') for p in plan['files']]
        plan['args'][1]=str(self.p/'missing-matrix')
        with self.assertRaisesRegex(ValueError,'preserved diagnostics'):run(plan,ROOT)
        self.assertFalse(Path(plan['final']).exists())
        self.assertTrue(list((root/'.bigwig_attempts/rna/failure').glob('*/attempt.json')))
        plan['args'][1]=str(root/'counts/5ends.raw.tsv')
        with patch('workflow_bigwig.subprocess.run',side_effect=KeyboardInterrupt):
            with self.assertRaises(KeyboardInterrupt): run(plan,ROOT)
        self.assertFalse(Path(plan['final']).exists())
        run(plan,ROOT)
        self.assertTrue(Path(plan['final']).is_dir())
        # A failed replacement must preserve the previous published tracks even
        # when the coordinator removes declared failed outputs.
        previous=(Path(plan['final'])/'tracks/A.plus.bw').read_bytes()
        plan['args'][1]=str(self.p/'missing-again')
        with self.assertRaises(ValueError): run(plan,ROOT)
        self.assertFalse(Path(plan['final']).exists())
        archives=list((root/'.bigwig_attempts/rna/failure').glob('*/previous_publication/tracks/A.plus.bw'))
        self.assertTrue(archives);self.assertEqual(archives[-1].read_bytes(),previous)
        plan['args'][1]=str(root/'counts/5ends.raw.tsv');run(plan,ROOT)
        for p,v in self.before.items():self.assertEqual((sha256(p),p.stat().st_mtime_ns),v)

    def test_method_inputs_and_dependency_fail_before_work(self):
        for method,options,error in [('deseq2',{},'explicitly'),('controls',{'deseq2_type':'ratio'},'controls'),
                                     ('supplied',{},'supplied'),('raw',{'depths':str(self.p/'absent')},'does not exist')]:
            self.block({method:options})
            result=self.workflow(True,False);self.assertIn(error,result)
            self.assertFalse((self.p/'results').exists())
        self.block({'raw':{}});self.cfg['bigwigs']['rna']['emit_samples']=['absent']
        self.assertIn('emit_samples',self.workflow(True,False))
        from bam_configuration import load_bam_settings
        settings=load_bam_settings(self.cfg,ROOT)
        self.block({'raw':{}})
        with patch('workflow_bigwig.importlib.import_module',side_effect=ImportError('missing')):
            with self.assertRaisesRegex(ValueError,'pyBigWig'): plans(self.cfg['bigwigs'],settings,ROOT,bam=True)
        self.block({'deseq2':{'deseq2_type':'ratio','rscript':'definitely-missing-Rscript'}})
        with self.assertRaisesRegex(ValueError,'Rscript'): plans(self.cfg['bigwigs'],settings,ROOT,bam=True)


    def test_supplied_and_actual_r_modes_exact_standalone(self):
        try: available=subprocess.run([R,'--vanilla','-e','quit(status=if(requireNamespace("DESeq2",quietly=TRUE)) 0 else 1)'],capture_output=True).returncode==0
        except OSError: available=False
        factors=self.p/'factors.tsv';factors.write_text('sample_id\tsize_factor\nB\t2\nA\t4\n')
        self.block({'supplied':{'size_factors':str(factors),'factor_kind':'size-factor'}})
        if available:
            controls=self.p/'controls.tsv';controls.write_text('chromosome\tstrand\tstart\tend\nrna\t+\t1\t1\nrna\t+\t6\t6\n')
            # Fixture has rna positions 1 and 6 on plus strand.
            self.cfg['bigwigs']['rna']['methods'].update({'deseq2':{'deseq2_type':'ratio','rscript':R},
                'controls':{'deseq2_type':'poscounts','controls':str(controls),'min_control_features':2,'rscript':R}})
        self.workflow()
        from bam_configuration import load_bam_settings
        resolved=plans(self.cfg['bigwigs'],load_bam_settings(self.cfg,ROOT),ROOT,bam=True)
        for key,plan in resolved.items():
            direct=self.p/('direct-'+key.split('/')[1])
            subprocess.run([sys.executable,str(ROOT/'workflow/scripts/bigwig2.py'),*plan['args'],'--output-dir',str(direct)],check=True,capture_output=True)
            for t in (Path(plan['final'])/'tracks').glob('*.bw'):self.assertEqual(t.read_bytes(),(direct/'tracks'/t.name).read_bytes())
            self.assertEqual((Path(plan['final'])/'normalization_factors.tsv').read_bytes(),(direct/'normalization_factors.tsv').read_bytes())
        if not available:self.skipTest('Supplied integration passed; R/DESeq2 integration unavailable')

    def test_separate_reference_pools_default_cohort_and_shared_strand_factor(self):
        import pyBigWig
        out=self.p/'pools'
        (out/'counts').mkdir(parents=True)
        (out/'reference').mkdir()
        (out/'reference/fixed.sizes.tsv').write_text('rna\t100\n')
        (out/'genome/B').mkdir(parents=True)
        (out/'genome/B/reference.sizes.tsv').write_text('rna\t100\n')
        contents='chromosome\tstrand\tstart\tend\tB\tA\nrna\t+\t1\t1\t4\t8\nrna\t-\t100\t100\t4\t8\n'
        (out/'counts/5ends.raw.tsv').write_text(contents)
        (out/'counts/genome.5ends.raw.tsv').write_text(contents)
        rna=self.p/'rna-factors.tsv';rna.write_text('sample_id\tsize_factor\nB\t2\nA\t4\n')
        genome=self.p/'genome-factors.tsv';genome.write_text('sample_id\tsize_factor\nA\t8\nB\t4\n')
        settings={'output':str(out),'samples':{'B':'b.bam','A':'a.bam'},'strand':'both',
                  'reference_mode':'fixed','genome':{'strand':'both'}}
        cfg={'rna':{'methods':{'supplied':{'size_factors':str(rna),'factor_kind':'size-factor'}},'emit_samples':['A']},
             'genome':{'cohort':['A','B'],'methods':{'supplied':{'size_factors':str(genome),'factor_kind':'size-factor'}}}}
        resolved=plans(cfg,settings,ROOT)
        for key,plan in resolved.items():
            run(plan,ROOT)
            meta=json.loads((Path(plan['final'])/'provenance.json').read_text())
            self.assertEqual(meta['cohort'],['B','A'] if key.startswith('rna') else ['A','B'])
            for strand,pos in [('plus',0),('minus',99)]:
                with pyBigWig.open(str(Path(plan['final'])/'tracks'/('A.'+strand+'.bw'))) as bw:
                    self.assertEqual(bw.values('rna',pos,pos+1),[2.0] if key.startswith('rna') else [1.0])
