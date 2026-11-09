"""Actual BigWig and R/DESeq2 synthetic tests. No experimental fixtures."""
import csv
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'workflow/scripts'))
from core import sha256
import test_core as fixtures

HAS_BIGWIG=importlib.util.find_spec('pyBigWig') is not None
RSCRIPT=os.environ.get('CAPZYME_BIGWIG_RSCRIPT',shutil.which('Rscript') or 'Rscript')
try:
    HAS_R=subprocess.run([RSCRIPT,'--vanilla','-e','quit(status=if(requireNamespace("DESeq2",quietly=TRUE)) 0 else 1)'],capture_output=True).returncode==0
except OSError: HAS_R=False


@unittest.skipUnless(HAS_BIGWIG,'pyBigWig required')
class BigWigTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix='bigwig test '); self.addCleanup(self.tmp.cleanup)
        self.p=Path(self.tmp.name)
        self.sizes=self.p/'sizes.tsv'; self.sizes.write_text('rna\t100\nvirus_minus\t50\n')
        self.matrix=self.p/'counts.tsv'
        self.header='chromosome\tstrand\tstart\tend\tB\tA\n'
        self.rows=[['rna','+',i,i,i*5,i*10] for i in range(1,11)]
        self.rows += [['rna','-',100,100,3,6],['virus_minus','+',1,1,7,14],['rna','+',20,20,100,1]]
        self.write_matrix(self.rows)
        self.controls=self.p/'controls.tsv'
        self.controls.write_text('chromosome\tstrand\tstart\tend\n'+''.join(f'rna\t+\t{i}\t{i}\n' for i in range(1,11)))
        self.before={p:sha256(p) for p in (self.sizes,self.matrix,self.controls)}
        self.counter=0

    def write_matrix(self,rows):
        self.matrix.write_text(self.header+''.join('\t'.join(map(str,r))+'\n' for r in rows))

    def invoke(self,method,*extra,fail=False,source=True):
        self.counter+=1; out=self.p/f'run{self.counter}'
        args=[sys.executable,str(ROOT/'workflow/scripts/bigwig2.py'),'--normalization',method,
              '--output-dir',str(out),'--rscript',RSCRIPT]
        if source: args += ['--matrix',str(self.matrix),'--assert-unthresholded','--chrom-sizes',str(self.sizes),'--cohort','B','A']
        args+=list(map(str,extra))
        r=subprocess.run(args,text=True,capture_output=True)
        if fail: self.assertNotEqual(r.returncode,0,r.stdout+r.stderr)
        else: self.assertEqual(r.returncode,0,r.stdout+r.stderr)
        return out,r

    def values(self,out,name,contig='rna'):
        import pyBigWig
        with pyBigWig.open(str(out/'tracks'/name)) as bw: return bw.values(contig,0,bw.chroms(contig))

    def test_raw_rpm_coordinates_strands_and_no_R(self):
        raw,_=self.invoke('raw','--rscript',self.p/'nonexistent-R')
        self.assertEqual(self.values(raw,'B.plus.bw')[0],5)
        self.assertEqual(self.values(raw,'B.minus.bw')[99],3)
        self.assertEqual(self.values(raw,'B.plus.bw','virus_minus')[0],7)
        self.assertTrue(math.isnan(self.values(raw,'B.plus.bw')[99]))
        rpm,_=self.invoke('rpm','--matrix-complete-both-strands','--strands','plus','--rscript',self.p/'nonexistent-R')
        self.assertFalse((rpm/'tracks/B.minus.bw').exists())
        self.assertAlmostEqual(self.values(rpm,'B.plus.bw')[0],5e6/(275+3+7+100),delta=0.001)
        meta=json.loads((rpm/'provenance.json').read_text())
        self.assertEqual(meta['cohort'],['B','A']); self.assertNotIn('DESeq2',meta['versions'])
        for p,h in self.before.items(): self.assertEqual(sha256(p),h)

    def test_matrix_rpm_requires_explicit_scope_and_zero_rejection(self):
        _,result=self.invoke('rpm',fail=True)
        self.assertIn('both-strands',result.stderr)
        depth=self.p/'depths.tsv'; depth.write_text('sample_id\taccepted_records\tscope_id\nB\t1000\tfull_library\nA\t2000\tfull_library\n')
        out,_=self.invoke('rpm','--depths',depth)
        self.assertEqual(self.values(out,'B.plus.bw')[0],5000)
        self.write_matrix([[*r[:4],0,0] for r in self.rows])
        _,result=self.invoke('rpm','--matrix-complete-both-strands',fail=True)
        self.assertIn('positive',result.stderr)
        raw,_=self.invoke('raw')
        self.assertEqual(self.values(raw,'B.plus.bw')[0],0)
        self.assertIn('B.plus.bw',json.loads((raw/'provenance.json').read_text())['empty_tracks'])

    def test_supplied_regeneration_order_context_and_errors(self):
        factors=self.p/'factors.tsv'; factors.write_text('sample_id\tsize_factor\nB\t2\nA\t4\n')
        first,_=self.invoke('supplied','--size-factors',factors,'--factor-kind','size-factor')
        self.assertEqual(self.values(first,'B.plus.bw')[0],2.5)
        self.matrix.write_text('chromosome\tstrand\tstart\tend\tB\n'+''.join('\t'.join(map(str,r[:5]))+'\n' for r in self.rows))
        second,_=self.invoke('supplied','--matrix',self.matrix,'--assert-unthresholded','--chrom-sizes',self.sizes,
            '--size-factors',first/'normalization_factors.tsv','--factor-kind','multiplier',
            '--factor-metadata',first/'provenance.json','--rscript',self.p/'no-R',source=False)
        self.assertEqual((first/'tracks/B.plus.bw').read_bytes(),(second/'tracks/B.plus.bw').read_bytes())
        self.assertEqual(json.loads((second/'provenance.json').read_text())['cohort'],['B','A'])
        saved_factors=first/'normalization_factors.tsv'
        saved_factors.write_text(saved_factors.read_text().replace('0.5','0.6'))
        _,bad=self.invoke('supplied','--matrix',self.matrix,'--assert-unthresholded','--chrom-sizes',self.sizes,
            '--size-factors',saved_factors,'--factor-kind','multiplier','--factor-metadata',first/'provenance.json',source=False,fail=True)
        self.assertIn('checksum mismatch',bad.stderr)
        factors.write_text('sample_id\tsize_factor\nB\tNaN\nA\t2\n')
        _,r=self.invoke('supplied','--size-factors',factors,'--factor-kind','size-factor',fail=True)
        self.assertIn('finite',r.stderr)

    def test_invalid_matrix_controls_options_and_output_paths(self):
        _,r=self.invoke('legacy-rpm',fail=True); self.assertIn('invalid choice',r.stderr)
        self.write_matrix([self.rows[0],self.rows[0]])
        _,r=self.invoke('raw',fail=True); self.assertIn('Duplicate',r.stderr)
        self.write_matrix([['rna','+',0,0,1,1]])
        _,r=self.invoke('raw',fail=True); self.assertIn('coordinate',r.stderr)
        self.write_matrix(self.rows)
        self.controls.write_text('chromosome\tstrand\tstart\tend\nrna\t+\t101\t101\n')
        _,r=self.invoke('controls','--controls',self.controls,fail=True); self.assertIn('outside reference',r.stderr)
        self.controls.write_text('chromosome\tstrand\tstart\tend\nrna\t+\t1\t1\nrna\t+\t1\t1\n')
        _,r=self.invoke('controls','--controls',self.controls,fail=True); self.assertIn('Duplicate control',r.stderr)

    def test_method_specific_rejection_and_missing_R_actionable(self):
        for method,extra in [('raw',['--min-control-features','10']),
                             ('raw',['--matrix-complete-both-strands']),
                             ('rpm',['--deseq2-type','ratio']),
                             ('deseq2',['--depths',str(self.p/'absent-depths')])]:
            _,r=self.invoke(method,*extra,fail=True)
            self.assertIn('applies only',r.stderr)
            args=[sys.executable,str(ROOT/'toolbox.py'),'bigwig','--matrix',str(self.matrix),
                  '--assert-unthresholded','--chrom-sizes',str(self.sizes),'--cohort','B','A',
                  '--normalization',method,'--output-dir',str(self.p/'preview'),*extra,'--preview']
            r=subprocess.run(args,text=True,capture_output=True)
            self.assertNotEqual(r.returncode,0); self.assertIn('applies only',r.stderr)
        out,r=self.invoke('deseq2','--deseq2-type','ratio','--rscript',self.p/'missing-Rscript',fail=True)
        self.assertIn('Rscript executable unavailable',r.stderr)
        self.assertIn('--rscript',r.stderr)
        self.assertEqual(json.loads((out/'provenance.json').read_text())['status'],'failed')

    @unittest.skipUnless(HAS_R,'actual DESeq2 required')
    def test_standalone_settings_deseq_ordered_cohort_subset_and_strand_scaling(self):
        direct,_=self.invoke('deseq2','--deseq2-type','ratio','--cohort','A','B','--emit-samples','B')
        cfg=self.p/'tools.json'; out=self.p/'settings-result'
        cfg.write_text(json.dumps({'schema_version':2,'input':{'type':'counts','matrix':str(self.matrix),'assert_unthresholded':True},
            'output_dir':str(out),'references':{'rna':{'chrom_sizes':str(self.sizes)}},
            'bigwig':{'rna':{'normalization':'deseq2','deseq2_type':'ratio','cohort':['A','B'],'emit_samples':['B'],'strands':'both'}}}))
        standalone_env=dict(os.environ,PATH=str(Path(RSCRIPT).parent)+os.pathsep+'/usr/bin:/bin')
        r=subprocess.run([sys.executable,str(ROOT/'toolbox.py'),'bigwig','--settings',str(cfg),'--rscript',RSCRIPT],text=True,capture_output=True,env=standalone_env)
        self.assertEqual(r.returncode,0,r.stdout+r.stderr)
        self.assertEqual((direct/'normalization_factors.tsv').read_bytes(),(out/'normalization_factors.tsv').read_bytes())
        self.assertEqual(sorted(p.name for p in (out/'tracks').iterdir()),['B.minus.bw','B.plus.bw'])
        for name in ('B.plus.bw','B.minus.bw'): self.assertEqual((direct/'tracks'/name).read_bytes(),(out/'tracks'/name).read_bytes())
        meta=json.loads((out/'provenance.json').read_text())
        self.assertEqual(meta['cohort'],['A','B']); self.assertEqual(meta['emitted_samples'],['B']); self.assertEqual(meta['estimator'],'ratio')
        with (out/'normalization_factors.tsv').open() as fh:
            factors=next(r for r in csv.DictReader(fh,delimiter='\t') if r['sample_id']=='B')
        self.assertAlmostEqual(self.values(out,'B.plus.bw')[0],5/float(factors['size_factor']),delta=1e-6)
        self.assertAlmostEqual(self.values(out,'B.minus.bw')[99],3/float(factors['size_factor']),delta=1e-6)

    @unittest.skipUnless(HAS_R,'actual DESeq2 required')
    def test_actual_R_ratio_poscounts_controls_and_direct_DESeq2_agreement(self):
        for mode,kind in [('deseq2','ratio'),('deseq2','poscounts'),('controls','ratio'),('controls','poscounts')]:
            out,_=self.invoke(mode,'--deseq2-type',kind,*(['--controls',self.controls] if mode=='controls' else []))
            # Independent direct call to DESeq2, not a Python estimator or the shared wrapper.
            script=f'''suppressPackageStartupMessages(library(DESeq2)); x <- as.matrix(read.delim({json.dumps(str(out/'normalization/estimation_counts.tsv'))},row.names=1,check.names=FALSE));
            d <- DESeqDataSetFromMatrix(x,data.frame(row.names=colnames(x),g=c("a","a")),~1);
            idx <- {'scan('+json.dumps(str(out/'normalization/control_indices.txt'))+',quiet=TRUE)' if mode=='controls' else 'NULL'};
            d <- if (is.null(idx)) estimateSizeFactors(d,type={json.dumps(kind)}) else estimateSizeFactors(d,type={json.dumps(kind)},controlGenes=idx);
            write.table(data.frame(sample_id=colnames(x),size_factor=sizeFactors(d)),{json.dumps(str(self.p/'direct.tsv'))},sep="\\t",row.names=FALSE,quote=FALSE)'''
            r=subprocess.run([RSCRIPT,'--vanilla','-e',script],text=True,capture_output=True)
            self.assertEqual(r.returncode,0,r.stderr)
            with (self.p/'direct.tsv').open() as fh: direct=list(csv.DictReader(fh,delimiter='\t'))
            with (out/'normalization_factors.tsv').open() as fh: actual=list(csv.DictReader(fh,delimiter='\t'))
            for a,b in zip(actual,direct): self.assertAlmostEqual(float(a['size_factor']),float(b['size_factor']),places=12)
            meta=json.loads((out/'provenance.json').read_text())
            self.assertEqual(meta['estimator'],kind); self.assertIn('not RPM',meta['units'])
            if mode=='controls':
                self.assertAlmostEqual(float(actual[1]['size_factor'])/float(actual[0]['size_factor']),2,places=12)
                self.assertAlmostEqual(self.values(out,'B.plus.bw')[19],100/float(actual[0]['size_factor']),delta=0.0001)
                self.assertTrue(all(int(d['contributing_features'])==10 for d in meta['normalization_diagnostics']))

    @unittest.skipUnless(HAS_R,'actual DESeq2 required')
    def test_sparse_ratio_failure_poscounts_success_and_insufficient_controls(self):
        sparse=[['rna','+',i,i,(i+1)*2 if i%2 else 0,(i+1)*3 if not i%2 else 0] for i in range(1,21)]
        self.write_matrix(sparse)
        failed,r=self.invoke('deseq2','--deseq2-type','ratio',fail=True)
        self.assertIn('no estimator fallback',r.stderr)
        self.assertIn('Insufficient informative', (failed/'commands.log').read_text())
        out,_=self.invoke('deseq2','--deseq2-type','poscounts')
        self.assertTrue(all(int(d['contributing_features'])==10 for d in json.loads((out/'provenance.json').read_text())['normalization_diagnostics']))
        self.controls.write_text('chromosome\tstrand\tstart\tend\nrna\t+\t1\t1\nrna\t+\t99\t99\n')
        failed,r=self.invoke('controls','--controls',self.controls,'--deseq2-type','poscounts',fail=True)
        meta=json.loads((failed/'provenance.json').read_text())
        self.assertEqual(len(meta['control_diagnostics']['missing']),1)
        self.assertEqual(meta['status'],'failed')

    @unittest.skipUnless(HAS_R,'actual DESeq2 required')
    def test_control_missing_zero_diagnostics_with_sufficient_valid_controls(self):
        self.write_matrix([*self.rows,['rna','+',99,99,0,0]])
        self.controls.write_text(self.controls.read_text()+'rna\t+\t98\t98\nrna\t+\t99\t99\n')
        out,_=self.invoke('controls','--controls',self.controls,'--deseq2-type','ratio')
        meta=json.loads((out/'provenance.json').read_text())
        self.assertEqual(meta['control_diagnostics']['requested'],12)
        self.assertEqual(len(meta['control_diagnostics']['missing']),1)
        self.assertEqual(len(meta['control_diagnostics']['all_zero']),1)
        self.assertTrue(all(int(d['contributing_features'])==10 for d in meta['normalization_diagnostics']))

    @unittest.skipUnless(shutil.which('samtools') and shutil.which('bedtools'),'BAM tools required')
    def test_bam_route_reuses_counter_and_preserves_sources(self):
        f=fixtures.CoordinateTests('test_empty_matrix_and_total_threshold'); f.setUp(); self.addCleanup(f.doCleanups)
        source=f.bam(mapq=0); before=sha256(source); refbefore=sha256(f.fasta)
        out,_=self.invoke('rpm','--bam',source,'--sample-name','B','--reference',f.fasta,'--cohort','B',source=False)
        self.assertAlmostEqual(self.values(out,'B.plus.bw')[0],2e6/9,delta=0.02)
        self.assertEqual(sha256(source),before); self.assertEqual(sha256(f.fasta),refbefore)
        self.assertFalse(Path(str(source)+'.bai').exists())
        self.assertFalse(json.loads((out/'provenance.json').read_text())['preparation'][0]['sort_performed'])
