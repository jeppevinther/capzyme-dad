"""Actual tiny-genome aligner/counter integration, without experimental inputs."""
import csv
import io
import json
from pathlib import Path
import random
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'workflow/scripts'))
from core import sha256
from genome_common import certain_cigar, fastq_records, ingest_fastq, process_sam
from types import SimpleNamespace


def rc(seq):return seq.translate(str.maketrans('ACGT','TGCA'))[::-1]
def write_reads(path,reads):path.write_text(''.join(f'@{n}\n{s}\n+\n'+('I'*len(s))+'\n' for n,s in reads))
def tsv(path):
    with open(path) as f:return list(csv.DictReader(f,delimiter='\t'))


@unittest.skipUnless(all(shutil.which(x) for x in ['bowtie2','bowtie2-build','bowtie2-inspect','hisat2','hisat2-build','hisat2-inspect','samtools','bedtools','hisat2_extract_splice_sites.py','hisat2_extract_exons.py']),'Actual aligners/builders/counter dependencies required')
class GenomeMappingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory(prefix='genome fixture ');cls.p=Path(cls.temp.name)
        rng=random.Random(253);dna=lambda n:''.join(rng.choice('ACGT') for _ in range(n))
        chrom=list(dna(1200));chrom[220:222]='GT';chrom[418:420]='AG';cls.chrom=''.join(chrom);cls.plasmid=dna(500);cls.decoy=dna(500);cls.repeat=dna(300)
        cls.fasta=cls.p/'genome.fa';cls.fasta.write_text(''.join(f'>{n} AC:synthetic description\n{s}\n' for n,s in [('chr',cls.chrom),('plasmid',cls.plasmid),('hs38d1',cls.decoy),('repeatA',cls.repeat),('repeatB',cls.repeat)]))
        cls.decoys=cls.p/'decoys.txt';cls.decoys.write_text('hs38d1\n')
        cls.gtf=cls.p/'genes.gtf';cls.gtf.write_text(''.join(f'chr\tsynthetic\texon\t{a}\t{b}\t.\t+\t.\tgene_id "g"; transcript_id "t";\n' for a,b in [(101,220),(421,540)]))
        cls.run_cli('bowtie2','prepare','--reference',cls.fasta,'--output-dir',cls.p/'btindex')
        cls.run_cli('hisat2','prepare','--reference',cls.fasta,'--gtf',cls.gtf,'--decoys',cls.decoys,'--output-dir',cls.p/'htindex')
        cls.bt=str(cls.p/'btindex/index/genome');cls.ht=str(cls.p/'htindex/index/genome')
        cls.original={p:sha256(p) for p in cls.p.rglob('*') if p.is_file()}

    @classmethod
    def tearDownClass(cls):cls.temp.cleanup()

    @classmethod
    def run_cli(cls,kind,*args,fail=False):
        cmd=[sys.executable,str(ROOT/'workflow/scripts'/f'{kind}_genome.py'),*map(str,args)]
        r=subprocess.run(cmd,capture_output=True,text=True)
        if fail:
            if r.returncode==0:raise AssertionError('Unexpected success: '+str(cmd))
        elif r.returncode:raise AssertionError(r.stdout+r.stderr)
        return r

    def setUp(self):
        self.t=tempfile.TemporaryDirectory(prefix='genome map test ');self.addCleanup(self.t.cleanup);self.d=Path(self.t.name)

    def reads(self,splice=False):
        reads=[('plus',self.chrom[600:660]),('minus',rc(self.chrom[700:760])),('left_boundary',self.chrom[:60]),('right_boundary',rc(self.chrom[-60:])),('plasmid',self.plasmid[40:100]),('decoy',self.decoy[120:180]),('multi',self.repeat[40:100]),('unmapped','N'*60)]
        if splice:
            junction=self.chrom[170:220]+self.chrom[420:470]
            reads += [('splice_plus',junction),('splice_minus',rc(junction))]
        return reads

    def mapping(self,kind,reads,reference=True,**extra):
        fq=self.d/(kind+str(len(list(self.d.iterdir())))+'.fq');write_reads(fq,reads);out=fq.with_suffix('.out')
        args=['map','--fastq',fq,'--index',self.bt if kind=='bowtie2' else self.ht,'--decoys',self.decoys,'--output-dir',out]
        if reference:args+=['--reference',self.fasta]
        for k,v in extra.items():args+=['--'+k.replace('_','-'),v]
        self.run_cli(kind,*args);return out,fq

    def matrix(self,out):
        dest=out/'counts.tsv'
        r=subprocess.run([sys.executable,str(ROOT/'workflow/scripts/count_5ends.py'),str(out/'accepted.bam'),'-g',str(self.fasta),'--sample-name','sample','--output',str(dest)],text=True,capture_output=True)
        self.assertEqual(r.returncode,0,r.stderr)
        return {(r['chromosome'],r['strand'],int(r['start'])):int(r['sample']) for r in tsv(dest)}

    def sam(self,path):return subprocess.check_output(['samtools','view',str(path)],text=True).splitlines()

    def test_bowtie_both_strands_multiple_contigs_boundaries_counter_and_handoff(self):
        out,fq=self.mapping('bowtie2',self.reads())
        self.assertEqual(self.matrix(out),{('chr','+',601):1,('chr','-',760):1,('chr','+',1):1,('chr','-',1200):1,('plasmid','+',41):1})
        qc=json.loads((out/'qc.json').read_text());self.assertEqual(qc['primary_outcomes'],{'accepted':5,'decoy':1,'low_mapq':1,'unmapped':1});self.assertEqual(qc['decoy_mapped_primary'],1)
        self.assertEqual([n for n,s,q in fastq_records(out/'unmapped.fastq')],['unmapped'])
        args=[r['command'] for r in map(json.loads,[l for l in (out/'commands.jsonl').read_text().splitlines() if l.startswith('{"command":')]) if 'command' in r]
        align=next(a for a in args if a[0]=='bowtie2');self.assertIn('--end-to-end',align);self.assertNotIn('--norc',align)
        for line in self.sam(out/'alignments.bam'):
            f=line.split('\t')
            if not int(f[1])&4:self.assertNotRegex(f[5],'[SH]')
        self.assertEqual(sha256(fq),json.loads((out/'provenance.json').read_text())['inputs_sha256'][str(fq.resolve())])

    def test_hisat_splice_plus_minus_and_actual_counter(self):
        out,fq=self.mapping('hisat2',self.reads(True))
        counts=self.matrix(out)
        expected={('chr','+',601):1,('chr','-',760):1,('chr','+',1):1,('chr','-',1200):1,('plasmid','+',41):1,('chr','+',171):1,('chr','-',470):1}
        self.assertEqual(counts,expected)
        mapped=[l.split('\t') for l in self.sam(out/'accepted.bam')]
        self.assertEqual({f[0]:f[5] for f in mapped if f[0].startswith('splice')},{'splice_plus':'50M200N50M','splice_minus':'50M200N50M'})
        self.assertTrue(all('S' not in f[5] and 'H' not in f[5] for f in mapped))
        qc=json.loads((out/'qc.json').read_text());self.assertEqual(qc['primary_outcomes']['low_mapq'],1);self.assertEqual(qc['primary_outcomes']['decoy'],1)
        self.assertGreater(qc['alignment_records'].get('nonprimary',0),0)

    def test_bacterial_genome_only_fasta_build_and_index_only_map(self):
        fq=self.d/'bacterial.fq';write_reads(fq,[('one',self.plasmid[10:70])]);out=self.d/'genome_only'
        self.run_cli('bowtie2','map','--reference',self.fasta,'--fastq',fq,'--output-dir',out)
        self.assertTrue((out/'index/genome.1.bt2').exists());self.assertEqual(self.matrix(out),{('plasmid','+',11):1})
        index_only,_=self.mapping('bowtie2',[('one',self.chrom[600:660])],reference=False)
        self.assertEqual(len(tsv(index_only/'mapping_by_reference.tsv')),5)

    def test_sequential_ncrna_unmapped_only_not_low_mapq(self):
        rna=self.d/'rna.fa';rna.write_text('>r1\n'+self.chrom[600:660]+'\n>r2\n'+self.repeat+'\n>r3\n'+self.repeat+'\n')
        fq=self.d/'first.fq';write_reads(fq,[('rna',self.chrom[600:660]),('ambiguous',self.repeat[40:100]),('genome',self.chrom[700:760])])
        first=self.d/'first';self.run_cli('bowtie2','map','--reference',rna,'--fastq',fq,'--output-dir',first)
        self.assertEqual([r[0] for r in fastq_records(first/'unmapped.fastq')],['genome'])
        second=self.d/'second';self.run_cli('hisat2','map','--index',self.ht,'--reference',self.fasta,'--fastq',first/'unmapped.fastq','--decoys',self.decoys,'--output-dir',second)
        accepted1={l.split('\t')[0] for l in self.sam(first/'accepted.bam')};accepted2={l.split('\t')[0] for l in self.sam(second/'accepted.bam')}
        self.assertEqual(accepted1,{'rna'});self.assertEqual(accepted2,{'genome'});self.assertFalse(accepted1&accepted2)
        self.assertEqual(json.loads((first/'qc.json').read_text())['primary_outcomes']['low_mapq'],1)

    def test_clipping_restrictions_actual_reads(self):
        reads=[('clean',self.chrom[600:660]),('adapter','T'*35+self.chrom[600:660])]
        for kind in ['bowtie2','hisat2']:
            out,_=self.mapping(kind,reads)
            self.assertEqual([r[0] for r in fastq_records(out/'unmapped.fastq')],['adapter'])
            self.assertEqual({f.split('\t')[0] for f in self.sam(out/'accepted.bam')},{'clean'})
            self.assertEqual(self.matrix(out),{('chr','+',601):1})

    def test_input_immutability_invalid_reference_decoys_and_outputs(self):
        out,fq=self.mapping('hisat2',[('one',self.chrom[600:660])])
        for p,h in self.original.items():self.assertEqual(sha256(p),h,str(p))
        self.assertFalse(Path(str(self.fasta)+'.fai').exists())
        for extra in [['--decoys',self.d/'missing.txt'],['--reference',self.fasta,'--reference-sha256','0'*64]]:
            self.run_cli('hisat2','map','--fastq',fq,'--index',self.ht,'--output-dir',self.d/('bad'+str(len(list(self.d.iterdir())))),*extra,fail=True)
        self.run_cli('bowtie2','map','--fastq',fq,'--index',self.bt,'--output-dir',out,fail=True)
        bad=self.d/'bad.fa';bad.write_text('>wrong\nACGT\n')
        self.run_cli('bowtie2','map','--reference',bad,'--fastq',fq,'--index',self.bt,'--output-dir',self.d/'badref',fail=True)

    @unittest.skipUnless(__import__('importlib.util',fromlist=['find_spec']).find_spec('pyBigWig'),'pyBigWig required')
    def test_bigwig_existing_bam_compatibility(self):
        import pyBigWig
        out,_=self.mapping('hisat2',[('plus',self.chrom[600:660]),('minus',rc(self.chrom[700:760]))])
        bw=self.d/'bw'
        r=subprocess.run([sys.executable,str(ROOT/'workflow/scripts/bigwig2.py'),'--bam',str(out/'accepted.bam'),'--sample-name','G','--reference',str(self.fasta),'--normalization','raw','--cohort','G','--output-dir',str(bw)],capture_output=True,text=True)
        self.assertEqual(r.returncode,0,r.stderr)
        with pyBigWig.open(str(bw/'tracks/G.plus.bw')) as plus,pyBigWig.open(str(bw/'tracks/G.minus.bw')) as minus:
            self.assertEqual(plus.values('chr',600,601),[1.0]);self.assertEqual(minus.values('chr',759,760),[1.0])

    def test_gzip_fastq_and_explicit_no_decoys(self):
        import gzip
        fq=self.d/'reads.fq.gz'
        with gzip.open(fq,'wt') as f:f.write('@decoy\n'+self.decoy[120:180]+'\n+\n'+'I'*60+'\n')
        out=self.d/'no_decoys'
        self.run_cli('hisat2','map','--index',self.ht,'--fastq',fq,'--reference',self.fasta,'--no-decoys','--output-dir',out)
        self.assertEqual(self.matrix(out),{('hs38d1','+',121):1})

    def test_explicit_annotation_tables_and_single_exon_gtf(self):
        out=self.d/'tables'
        self.run_cli('hisat2','prepare','--reference',self.fasta,'--splice-sites',self.p/'htindex/splice_sites.tsv','--exons',self.p/'htindex/exons.tsv','--decoys',self.decoys,'--output-dir',out)
        self.assertEqual(json.loads((out/'provenance.json').read_text())['status'],'complete')
        gtf=self.d/'single.gtf';gtf.write_text('chr\ts\texon\t101\t220\t.\t+\t.\tgene_id "g"; transcript_id "t";\n')
        single=self.d/'single';r=self.run_cli('hisat2','prepare','--reference',self.fasta,'--gtf',gtf,'--no-decoys','--output-dir',single,fail=True)
        self.assertIn('Exon-only/no-junction',r.stderr)

    def test_hisat_preparation_annotation_mismatch_and_missing_decoy_declaration(self):
        bad=self.d/'bad.gtf';bad.write_text('unknown\ts\texon\t1\t20\t.\t+\t.\tgene_id "g"; transcript_id "t";\n')
        self.run_cli('hisat2','prepare','--reference',self.fasta,'--gtf',bad,'--no-decoys','--output-dir',self.d/'badannot',fail=True)
        fq=self.d/'x.fq';write_reads(fq,[('x',self.chrom[600:660])])
        self.run_cli('hisat2','map','--index',self.ht,'--fastq',fq,'--output-dir',self.d/'missing_policy',fail=True)
        self.run_cli('hisat2','map','--reference',self.fasta,'--fastq',fq,'--no-decoys','--output-dir',self.d/'missing_index',fail=True)


class GenomeInputPolicyTests(unittest.TestCase):
    def test_inspector_description_is_not_part_of_reference_identifier(self):
        from unittest.mock import patch
        from genome_common import dictionary
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);fasta=p/'g.fa';fasta.write_text('>chr1 AC:example rl:Chromosome\nACGT\n')
            result=SimpleNamespace(stdout='Sequence-1\tchr1 AC:example rl:Chromosome\t4\n')
            with patch('genome_common.command',return_value=result):
                self.assertEqual(dictionary('unused','bowtie2',p,io.StringIO(),fasta),{'chr1':4})
            self.assertEqual((p/'reference.sizes.tsv').read_text(),'chr1\t4\n')

    def test_complete_index_families_and_ambiguity(self):
        from genome_common import index_files
        with tempfile.TemporaryDirectory() as tmp:
            prefix=str(Path(tmp)/'genome')
            for i in range(1,9):Path(prefix+f'.{i}.ht2l').touch()
            self.assertEqual(len(index_files(prefix,'hisat2')),8)
            Path(prefix+'.1.ht2').touch()
            with self.assertRaisesRegex(ValueError,'unambiguous'):index_files(prefix,'hisat2')
            b=str(Path(tmp)/'bowtie')
            for e in ['.1','.2','.3','.4','.rev.1','.rev.2']:Path(b+e+'.bt2l').touch()
            self.assertEqual(len(index_files(b,'bowtie2')),6)
            Path(b+'.2.bt2l').unlink()
            with self.assertRaisesRegex(ValueError,'complete'):index_files(b,'bowtie2')

    def test_cigar_certainty_and_bounds(self):
        sizes={'chr':1000}
        for cigar,seq in [('30M100N30M','A'*60),('30M2D30M','A'*60),('30M2I30M','A'*62)]:self.assertTrue(certain_cigar(cigar,1,sizes,'chr',seq))
        for cigar,seq in [('2S58M','A'*60),('60M2H','A'*60),('2I58M','A'*60),('2D60M','A'*60),('60M1P','A'*60),('1N60M','A'*60)]:self.assertFalse(certain_cigar(cigar,1,sizes,'chr',seq))
        for cigar,pos in [('0M',1),('60M',950),('3Z',1)]:
            with self.assertRaises(ValueError):certain_cigar(cigar,pos,sizes,'chr','A'*60)

    def test_duplicate_names_and_malformed_fastq(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'x.fq';write_reads(p,[('x','A'*20),('x','C'*20)])
            with sqlite3.connect(':memory:') as db:
                with self.assertRaisesRegex(ValueError,'unique'):ingest_fastq(p,db)
            p.write_text('@x\nAAAA\n+\nII\n')
            with self.assertRaises(ValueError):list(fastq_records(p))

    def test_filtering_duplicates_qcfail_uncertainty_and_fail_fast(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);reads=[('qc','A'*60),('dup','A'*60),('low','A'*60),('clip','A'*60)]
            fq=p/'x.fq';write_reads(fq,reads)
            sam=p/'x.sam';sam.write_text('@HD\tVN:1.6\n@SQ\tSN:chr\tLN:1000\n'+''.join(f'{n}\t{flag}\tchr\t10\t{mapq}\t{cigar}\t*\t0\t0\t'+('A'*60)+'\t'+('I'*60)+'\n' for n,flag,mapq,cigar in [('qc',512,20,'60M'),('dup',1024,20,'60M'),('low',0,19,'60M'),('clip',0,60,'2S58M')]))
            for exclude,uncertain,expected in [(False,'reject',2),(True,'reject',1),(False,'error',None)]:
                with sqlite3.connect(':memory:') as db:
                    ingest_fastq(fq,db);a=SimpleNamespace(mapq=20,exclude_duplicates=exclude,uncertain_cigar=uncertain)
                    if expected is None:
                        with self.assertRaisesRegex(ValueError,'uncertain'):process_sam(sam,p/'out.sam',{'chr':1000},set(),a,db)
                    else:self.assertEqual(process_sam(sam,p/'out.sam',{'chr':1000},set(),a,db)['primary_outcomes']['accepted'],expected)

if __name__=='__main__':unittest.main()
