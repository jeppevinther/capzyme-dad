"""Sequential Snakemake tests with real small aligners and the shared counter."""
import copy
import csv
import io
import json
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import test_genome_mapping as fixture
from synthetic import ADAPTER

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'workflow/scripts'))
from core import sha256
from genome_common import fastq_records
from rna_handoff import process_rna_sam
from sequential import publish_stage

TOOLS = all(shutil.which(t) for t in ('snakemake','cutadapt','bowtie2','bowtie2-build','hisat2',
                                     'hisat2-build','samtools','bedtools','hisat2_extract_splice_sites.py',
                                     'hisat2_extract_exons.py'))


@unittest.skipUnless(TOOLS, 'Full sequential workflow dependencies required')
class SequentialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        fixture.GenomeMappingTests.setUpClass()
        cls.g = fixture.GenomeMappingTests

    @classmethod
    def tearDownClass(cls):
        cls.g.tearDownClass()

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='sequential fixture ')
        self.addCleanup(self.temp.cleanup)
        self.p = Path(self.temp.name)
        self.rna = self.p/'rna.fa'
        self.rna.write_text('>rna\n'+self.g.chrom[600:660]+'\n>repeat1\n'+self.g.repeat+'\n>repeat2\n'+self.g.repeat+'\n')
        r = self.p/'rna_index'
        self.g.run_cli('bowtie2', 'prepare', '--reference', self.rna, '--output-dir', r)
        names = [('rna',self.g.chrom[600:660]), ('low',self.g.repeat[40:100]),
                 ('plus',self.g.chrom[50:110]), ('minus',fixture.rc(self.g.chrom[700:760])),
                 ('splice_plus',self.g.chrom[170:220]+self.g.chrom[420:470]),
                 ('splice_minus',fixture.rc(self.g.chrom[170:220]+self.g.chrom[420:470])),
                 ('decoy',self.g.decoy[120:180]), ('unknown','N'*60)]
        for sample, reads in [('B',names),('A',[('plasmid',self.g.plasmid[40:100])])]:
            fixture.write_reads(self.p/(sample+'.fq'), [(n,seq+ADAPTER) for n,seq in reads])
        self.manifest=self.p/'samples.tsv'
        self.manifest.write_text('sample_id\tfastq\n'+''.join(f'{s}\t{self.p/(s+".fq")}\n' for s in ['B','A']))
        from fixed_reference import fixed_inputs
        prefix=str(r/'index/genome')
        parts=[str(p) for p in (r/'index').glob('genome.*bt2')]
        self.cfg={'schema_version':1,'samples':str(self.manifest),'output_dir':str(self.p/'out'),
            'reference':{'mode':'fixed','fasta':str(self.rna),'fasta_sha256':sha256(self.rna),
                         'index_prefix':prefix,'index_format':'bt2',
                         'index_sha256':{p[len(prefix):]:sha256(p) for p in parts}},
            'preprocessing':{'adapter':ADAPTER,'quality':20,'minimum_length':15,'nextseq_trim':20},
            'mapping':{'norc':True,'threads':4},'counting':{'strand':'both','mapq':20},
            'genome':{'aligner':'hisat2','reference':str(self.g.fasta),'index_prefix':self.g.ht,
                      'reference_sha256':sha256(self.g.fasta),'decoys':str(self.g.decoys),
                      'policy_files':[str(self.g.gtf)],'threads':4,'mem_mb':12000,'count_mem_mb':4000}}
        self.config=self.p/'config.json'
        self.write_config()

    def write_config(self):
        self.config.write_text(json.dumps(self.cfg))

    def run_workflow(self, *extra, snake='Snakefile', fail=False):
        cmd=['snakemake','-s',str(ROOT/'workflow'/snake),'--configfile',str(self.config),
             '--directory',str(self.p),'--cores','8','--resources','mem_mb=24000','--printshellcmds',*extra]
        r=subprocess.run(cmd,text=True,capture_output=True)
        diagnostics = '\n'.join(p.read_text()[-5000:] for p in (Path(self.cfg['output_dir'])/'logs').glob('*.jsonl')) if r.returncode else ''
        if fail:self.assertNotEqual(r.returncode,0,r.stdout+r.stderr)
        else:self.assertEqual(r.returncode,0,r.stdout+r.stderr+diagnostics)
        return r.stdout+r.stderr

    def test_exact_handoff_counts_policy_and_restart(self):
        before={p:sha256(p) for p in self.p.rglob('*') if p.is_file()}
        dry=self.run_workflow('--dry-run')
        self.assertIn('map_rna_handoff',dry);self.assertIn('map_genome',dry)
        self.assertIn('mem_mb=12000',dry)
        self.assertIn(self.g.ht+'.8.ht2',dry)
        self.run_workflow()
        out=self.p/'out'
        rna=out/'rna/B'; genome=out/'genome/B'
        self.assertEqual({n for n,s,q in fastq_records(rna/'unmapped.fastq')},
                         {'plus','minus','splice_plus','splice_minus','decoy','unknown'})
        self.assertNotIn('low',(rna/'unmapped.fastq').read_text())
        self.assertEqual(json.loads((rna/'qc.json').read_text())['primary_outcomes'],
                         {'accepted':1,'low_mapq':1,'unmapped':6})
        report=json.loads((out/'qc/B.sequential.json').read_text())
        self.assertEqual(report['accepted_read_id_overlap'],0)
        self.assertEqual(sum(report['exclusive_outcomes'].values()),8)
        self.assertEqual(sum(report['exclusive_raw_outcomes'].values()),8)
        self.assertEqual(report['raw_input_reads'],8)
        self.assertEqual(report['exclusive_outcomes']['genome_decoy'],1)
        counts={(r['chromosome'],r['strand'],int(r['start'])):(int(r['B']),int(r['A']))
                for r in fixture.tsv(out/'counts/genome.5ends.raw.tsv')}
        self.assertEqual(counts,{('chr','+',51):(1,0),('chr','-',760):(1,0),
                                 ('chr','+',171):(1,0),('chr','-',470):(1,0),('plasmid','+',41):(0,1)})
        rows=fixture.tsv(out/'counts/5ends.raw.tsv')
        self.assertEqual([(r['chromosome'],r['start'],r['B'],r['A']) for r in rows],[('rna','1','1','0')])
        cli=self.p/'standalone.tsv'
        subprocess.run([sys.executable,str(ROOT/'workflow/scripts/count_5ends.py'),
                        str(genome/'accepted.bam'),str(out/'genome/A/accepted.bam'),
                        '--sample-name','B','--sample-name','A','-g',str(self.g.fasta),'-o',str(cli)],
                       check=True,capture_output=True)
        self.assertEqual(cli.read_bytes(),(out/'counts/genome.5ends.raw.tsv').read_bytes())
        self.assertIn('Nothing to be done',self.run_workflow('--dry-run'))
        self.run_workflow('--forcerun','map_genome')
        self.assertTrue(any((out/'.history/genome/B').iterdir()))
        self.assertEqual(cli.read_bytes(),(out/'counts/genome.5ends.raw.tsv').read_bytes())
        for p,digest in before.items():self.assertEqual(sha256(p),digest,str(p))

    def test_failure_preserves_attempts_then_restart_and_complete_index_dependencies(self):
        decoys=self.p/'bad_decoys.txt';decoys.write_text('absent\n')
        self.cfg['genome']['decoys']=str(decoys);self.write_config()
        self.run_workflow(fail=True)
        attempts=list((self.p/'out/.attempts/genome').rglob('provenance.json'))
        self.assertTrue(attempts)
        self.assertTrue(any(json.loads(p.read_text())['status']=='failed' for p in attempts))
        self.cfg['genome']['decoys']=str(self.g.decoys);self.write_config()
        self.run_workflow('--rerun-incomplete')
        self.assertTrue((self.p/'out/provenance/run.json').is_file())
        # Missing final family part must fail during configuration, before scheduling.
        from genome_configuration import genome_settings
        from unittest.mock import patch
        import genome_common
        with patch.object(genome_common.Path,'is_file',lambda p: False if str(p).endswith('.8.ht2') else Path.exists(p)):
            with self.assertRaisesRegex(ValueError,'complete'):
                genome_settings(self.cfg['genome'],ROOT,self.p/'new')

    def test_bowtie_genome_only_orchestration_and_rna_only_compatibility(self):
        fixture.write_reads(self.p/'B.fq',[('plus',self.g.chrom[50:110]),('minus',fixture.rc(self.g.chrom[700:760]))])
        fixture.write_reads(self.p/'A.fq',[('plasmid',self.g.plasmid[40:100])])
        self.cfg={k:v for k,v in self.cfg.items() if k in ('schema_version','samples','output_dir','genome')}
        self.cfg['genome'].update(aligner='bowtie2',index_prefix=self.g.bt)
        self.cfg['output_dir']=str(self.p/'genome_only')
        self.write_config()
        dry=self.run_workflow('--dry-run',snake='genome.smk')
        self.assertNotIn('trim_fastq',dry);self.assertNotIn('map_rna_handoff',dry)
        self.run_workflow(snake='genome.smk')
        self.assertTrue((self.p/'genome_only/counts/genome.5ends.raw.tsv').is_file())
        counts={(r['chromosome'],r['strand'],int(r['start'])):(int(r['B']),int(r['A']))
                for r in fixture.tsv(self.p/'genome_only/counts/genome.5ends.raw.tsv')}
        self.assertEqual(counts,{('chr','+',51):(1,0),('chr','-',760):(1,0),('plasmid','+',41):(0,1)})
        self.assertIn('Nothing to be done',self.run_workflow('--dry-run',snake='genome.smk'))

    def test_empty_handoff_generated_reference_preserves_established_rna_matrix(self):
        from synthetic import generate
        path=generate(self.p/'generated')
        self.cfg=json.loads(path.read_text())
        self.cfg['genome']={'aligner':'hisat2','reference':str(self.g.fasta),
                            'index_prefix':self.g.ht,'decoys':str(self.g.decoys),'threads':4}
        self.write_config()
        self.run_workflow()
        out=Path(self.cfg['output_dir'])
        self.assertEqual((out/'counts/5ends.raw.tsv').read_bytes(),(self.p/'generated/expected.tsv').read_bytes())
        self.assertEqual((out/'counts/genome.5ends.raw.tsv').read_text(),
                         'chromosome\tstrand\tstart\tend\tB\tA\tempty\n')
        for name in ('B','A','empty'):
            report=json.loads((out/f'qc/{name}.sequential.json').read_text())
            self.assertEqual(report['forwarded_unmapped_reads'],0)
            self.assertEqual(report['exclusive_outcomes'].get('genome_accepted',0),0)


class SequentialSafetyTests(unittest.TestCase):
    def test_rna_flags_qcfail_duplicates_and_nonprimary_do_not_enter_handoff(self):
        from genome_common import ingest_fastq
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);fq=p/'reads.fq'
            fixture.write_reads(fq,[(n,'A'*30) for n in ('qc','dup','low','flag','unmap')])
            sam=p/'records.sam'
            sam.write_text('@SQ\tSN:rna\tLN:100\n'+''.join(
                f'{n}\t{flag}\trna\t1\t{q}\t30M\t*\t0\t0\t'+('A'*30)+'\t'+('I'*30)+'\n'
                for n,flag,q in [('qc',512,20),('dup',1024,20),('low',0,19),('flag',256,60),('flag',0,0),('unmap',4,0)]))
            with sqlite3.connect(':memory:') as db:
                ingest_fastq(fq,db)
                qc=process_rna_sam(sam,db,20,2308)
                self.assertEqual(qc['primary_outcomes'],{'accepted':2,'low_mapq':2,'unmapped':1})
                self.assertEqual(db.execute('select id from reads where unmapped=1').fetchall(),[('unmap',)])

    def test_publish_refuses_unmanaged_results_and_preserves_failed_attempts(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);dest=p/'genome/A';dest.mkdir(parents=True);(dest/'user.txt').write_text('keep')
            def fake(args, log):
                out=Path(args[-1]);out.mkdir()
                (out/'provenance.json').write_text(json.dumps({'status':'complete','inputs_unchanged':True}))
            with patch('sequential.command',side_effect=fake):
                with self.assertRaisesRegex(ValueError,'unmanaged'):
                    publish_stage(p,'genome','A',['fake'],io.StringIO())
            self.assertEqual((dest/'user.txt').read_text(),'keep')


if __name__=='__main__':unittest.main()
