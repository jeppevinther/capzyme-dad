"""Record-only reporting: no experimental data, BAM readers or scientific tools."""
import copy
import csv
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest
from urllib.parse import unquote
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'workflow/scripts'))
from run_report import generate_report, workflow_context, FAILED, RECORDED


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def put(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2)+'\n')


class ReportTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='report fixture ');self.addCleanup(self.temp.cleanup)
        self.p=Path(self.temp.name).resolve();self.root=self.p/'analysis';self.root.mkdir();self.names=['02','01']
        s={'output':str(self.root),'samples':dict.fromkeys(self.names,'external-input'), 'pre':{'quality':20,'adapter':'<script>alert("x")</script>'},
            'genome':{'strand':'both'},'strand':'plus','minimum':0,'mapq':20,'flags':2308,'reference_mode':'fixed','fixed':{'fasta':'custom.fa'}}
        self.context=workflow_context(s,{}, {}, [''+str(self.root/'counts/5ends.raw.tsv'),str(self.root/'counts/genome.5ends.raw.tsv'),str(self.root/'provenance/run.json')], 'sequential')
        put(self.root/'provenance/run.json',{'settings':s,'sample_order':self.names,'git':{'commit':'reviewed-source'}})
        put(self.root/'provenance/genome.json',{'settings':s,'sample_order':self.names})
        (self.root/'counts').mkdir()
        for name in ['5ends.raw.tsv','genome.5ends.raw.tsv']:(self.root/'counts'/name).write_text('chromosome\tstrand\tstart\tend\t02\t01\nchr\t+\t1\t1\t3\t3\n')
        for name in self.names:self.library(name)
        (self.root/'logs').mkdir(exist_ok=True)
        put(self.root/'logs/count_5ends.json',{'step':'matrix','exit_status':0,'summary':{'samples':self.names,'accepted_ends_per_sample':[3,3]}})
        (self.root/'logs/genome.count.jsonl').write_text('exit_status=0\n')

    def stage(self,space,name,outcomes,total):
        d=self.root/space/name;d.mkdir(parents=True)
        (d/'accepted.bam').write_bytes(b'not-a-BAM-never-open-this')
        (d/'alignments.bam').write_bytes(b'no-scans')
        (d/'unmapped.fastq').write_text('not-scanned')
        q={'primary_outcomes':outcomes,'alignment_records':{'all':total+1,'nonprimary':1}}
        put(d/'qc.json',q)
        fastq=str(self.root/'rna'/name/'unmapped.fastq') if space=='genome' else 'external.trimmed.fq'
        m={'status':'complete','inputs_unchanged':True,'input_reads':total,'qc':q,'parameters':{'output_dir':str(d),'fastq':fastq,'reference':'ref.fa'},'inputs_sha256':{fastq:'a'*64},
           'outputs_sha256':{'accepted.bam':'b'*64,'alignments.bam':'c'*64,'unmapped.fastq':'a'*64,'qc.json':sha(d/'qc.json')}}
        put(d/'provenance.json',m)
        put(d/'publication.json',{'standalone_output_dir':str(d),'published_output_dir':str(d),'standalone_provenance_sha256':sha(d/'provenance.json')})
        return m

    def library(self,name):
        self.stage('rna',name,{'accepted':3,'low_mapq':1,'unmapped':4},8)
        self.stage('genome',name,{'accepted':2,'decoy':1,'unmapped':1},4)
        put(self.root/'qc'/f'{name}.cutadapt.json',{'read_counts':{'input':10,'output':8}})
        pins={space:sha(self.root/space/name/'provenance.json') for space in ('rna','genome')}
        x={'status':'complete','trimmed_input_reads':8,'raw_input_reads':10,'trimming_removed_reads':2,'forwarded_unmapped_reads':4,
            'exclusive_outcomes':{'rna_accepted':3,'rna_low_mapq':1,'genome_accepted':2,'genome_decoy':1,'genome_unmapped':1},
            'exclusive_raw_outcomes':{'rna_accepted':3,'rna_low_mapq':1,'genome_accepted':2,'genome_decoy':1,'genome_unmapped':1,'trimming_removed':2},
            'accepted_read_id_overlap':0,'handoff_sequence_quality_verified':True,'stage_provenance_sha256':pins,
            'handoff_integrity':{'rna_output_sha256':'a'*64,'genome_input_sha256':'a'*64,'current_handoff_sha256':'a'*64,'stage_provenance_sha256':pins},
            'trim_report_sha256':sha(self.root/'qc'/f'{name}.cutadapt.json')}
        put(self.root/'qc'/f'{name}.sequential.json',x)

    def run_report(self,context=True):
        out=self.p/'new report #1';code=generate_report(self.root,out,self.context if context else None)
        return code,json.loads((out/'report_status.json').read_text()),(out/'report.html').read_text(),list(csv.DictReader((out/'sample_summary.tsv').read_text().splitlines(),delimiter='\t'))

    def test_order_zero_ids_denominators_recorded_vs_fresh_no_large_scan(self):
        read=Path.read_bytes
        def safe(p):
            self.assertNotIn(p.suffix,('.bam','.bw','.fq','.fastq'));return read(p)
        with patch.object(Path,'read_bytes',safe):code,status,page,rows=self.run_report()
        self.assertEqual(code,0);self.assertEqual([r['sample_id'] for r in rows],self.names)
        self.assertEqual(rows[0]['rna_accepted_reads_pct_of_raw_reads'],'30.0')
        self.assertEqual(rows[0]['rna_accepted_reads_pct_of_rna_input_reads'],'37.5')
        self.assertEqual(rows[0]['genome_mapped_reads_pct_of_genome_input_reads'],'75.0')
        self.assertEqual(rows[0]['genome_nonprimary_records'],'1')
        self.assertNotIn('treatment',rows[0]);self.assertIn(RECORDED,page)
        self.assertEqual(status['workflow_exit_status'],'not assessed; a final workflow task cannot certify subsequent coordinator exit')

    def test_disabled_features_versus_missing_requested_outputs(self):
        code,_,page,_=self.run_report();self.assertEqual(code,0);self.assertIn('not enabled',page)
        self.context['requested_outputs'].append('bigwigs/rna/raw/tracks/02.plus.bw')
        code,status,_,_=self.run_report();self.assertEqual(code,1)
        self.assertTrue(any(c['status']==FAILED and '02.plus.bw' in c['check'] for c in status['checks']))

    def test_partial_run_informative_default_outputs_nonzero(self):
        p=self.p/'partial';p.mkdir()
        self.assertEqual(generate_report(p),1)
        self.assertIn('Required checks failed',(p/'report/report.html').read_text())
        (self.root/'genome/01/provenance.json').unlink()
        code,status,_,_=self.run_report();self.assertEqual(code,1)
        self.assertTrue((self.p/'new report #1/sample_summary.tsv').exists())

    def test_inconsistent_ids_totals_handoff_overlap_and_completion(self):
        for file,field,value in [('provenance/run.json','sample_order',['01','02']),('qc/02.sequential.json','trimmed_input_reads',9),
                                  ('qc/02.sequential.json','accepted_read_id_overlap',1),('qc/02.sequential.json','handoff_sequence_quality_verified',False),
                                  ('genome/02/provenance.json','status','failed')]:
            p=self.root/file;old=p.read_text();d=json.loads(old);d[field]=value;put(p,d)
            with self.subTest(file=file,field=field):self.assertEqual(self.run_report()[0],1)
            p.write_text(old)

    def test_genome_only_and_bam_nonapplicable_fields(self):
        for route in ('genome','bam'):
            ctx=copy.deepcopy(self.context);ctx['route']=route;s=ctx['settings'];s.pop('pre');s.pop('reference_mode');s.pop('fixed')
            if route=='genome':s.pop('strand');put(self.root/'provenance/genome.json',{'settings':s,'sample_order':self.names});(self.root/'provenance/run.json').unlink()
            else:
                s.pop('genome');s['fasta']='external.fa';put(self.root/'provenance/run.json',{'settings':s})
                for name in self.names:
                    dest=self.root/'bam'/f'{name}.accepted.bam';dest.parent.mkdir(exist_ok=True);dest.write_bytes(b'none')
                    put(self.root/'logs'/f'{name}.prepare.jsonl',{'step':'prepare','exit_status':0,'result':{'source_unchanged_stat':True}})
                put(self.root/'logs/provenance.jsonl',{'step':'provenance','exit_status':0,'result':{}})
                put(self.root/'logs/count.jsonl',{'step':'count','exit_status':0,'result':{'samples':self.names,'accepted_ends_per_sample':[3,3]}})
                (self.root/'counts/5ends.min_total_0.tsv').write_text((self.root/'counts/5ends.raw.tsv').read_text())
            ctx['requested_outputs']=[]
            out=self.p/route;self.assertEqual(generate_report(self.root,out,ctx),0)
            rows=list(csv.DictReader((out/'sample_summary.tsv').read_text().splitlines(),delimiter='\t'))
            self.assertEqual(rows[0]['raw_reads'],'');self.assertEqual(rows[0]['trimming_removed_reads'],'')
            self.assertIn('not applicable',(out/'report.html').read_text())

    def test_html_escaping_and_links_resolve_from_explicit_directory(self):
        code,_,page,_=self.run_report();self.assertEqual(code,0)
        self.assertIn('&lt;script&gt;',page);self.assertNotIn('<script>',page)
        for target in re.findall(r'href="([^"]+)"',page):
            p=(self.p/'new report #1'/unquote(target)).resolve();self.assertTrue(p.exists(),target)
        self.assertNotRegex(page,r'<(?:script|link)\b')

    def factors(self):
        self.context['bigwigs']={'genome':{'methods':{'rpm':{}},'strands':'both'}}
        d=self.root/'bigwigs/genome/rpm';(d/'tracks').mkdir(parents=True)
        f=d/'normalization_factors.tsv';f.write_text('sample_id\tsize_factor\tmultiplier\taccepted_records\tmethod\tunits\tcohort_id\n'+''.join(n+'\t\t500000\t2\trpm\tRPM\tpool\n' for n in self.names))
        outputs={'normalization_factors.tsv':sha(f)}
        for n in self.names:
            for strand in ('plus','minus'):
                rel='tracks/'+n+'.'+strand+'.bw';(d/rel).write_bytes(b'not-decoded');outputs[rel]='e'*64
        m={'status':'complete','normalization':'rpm','cohort':self.names,'emitted_samples':self.names,'parameters':{'strands':'both'},'units':'RPM','cohort_id':'pool','matrix_sha256':'c'*64,'reference_dictionary_sha256':'d'*64,'outputs_sha256':outputs}
        put(d/'provenance.json',m);return d,m

    def test_normalization_cohort_factors_integrity_and_missing_track(self):
        d,m=self.factors();self.assertEqual(self.run_report()[0],0)
        for key,value in [('cohort',['01','02']),('normalization','raw'),('status','failed')]:
            put(d/'provenance.json',dict(m,**{key:value}));self.assertEqual(self.run_report()[0],1)
        put(d/'provenance.json',m)
        f=d/'normalization_factors.tsv';old=f.read_text();f.write_text(old.replace('500000','600000'));self.assertEqual(self.run_report()[0],1)
        f.write_text(old);(d/'tracks/01.plus.bw').unlink();self.assertEqual(self.run_report()[0],1)

    def test_empty_receipt_missing_artifact_and_invalid_record_fail(self):
        path=self.root/'provenance/run.json';old=path.read_text();put(path,{})
        self.assertEqual(self.run_report()[0],1);path.write_text(old)
        bam=self.root/'genome/02/alignments.bam';bam.unlink();self.assertEqual(self.run_report()[0],1)
        path.write_text('{malformed');self.assertEqual(self.run_report()[0],1)

    def test_explicit_metadata_only(self):
        sheet=self.p/'metadata.tsv';sheet.write_text('sample_id\ttreatment\treplicate\n02\tcontrol\tb\n01\tRpp\ta\n')
        p=self.root/'provenance/run.json';d=json.loads(p.read_text());d['settings']['manifest']=str(sheet);put(p,d)
        code,_,_,rows=self.run_report();self.assertEqual(code,0);self.assertEqual(rows[0]['treatment'],'control');self.assertEqual(rows[1]['replicate'],'a')

    def test_record_only_regeneration_and_cli(self):
        before={str(p):(sha(p),p.stat().st_mtime_ns) for p in self.root.rglob('*') if p.is_file()}
        self.assertEqual(self.run_report()[0],0);self.assertEqual(self.run_report()[0],0)
        after={str(p):(sha(p),p.stat().st_mtime_ns) for p in self.root.rglob('*') if p.is_file()};self.assertEqual(before,after)
        r=subprocess.run([sys.executable,'-m','capzyme_cli','report','--help'],cwd=ROOT,text=True,capture_output=True)
        self.assertEqual(r.returncode,0,r.stderr);self.assertIn('--analysis-dir',r.stdout)


if __name__=='__main__':unittest.main()
