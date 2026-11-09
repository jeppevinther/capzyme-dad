"""Hash contract and real Cutadapt NextSeq-off regression, using synthetic reads."""
import json
import gzip
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'workflow/scripts'))
from tool_settings import parser, resolve, load
from reference_hashes import expected_digest, check_manifest
from unittest.mock import patch
from core import sha256

class ContractTests(unittest.TestCase):
    def test_hash_types_and_case(self):
        self.assertIsNone(expected_digest(None,'pin'))
        self.assertEqual(expected_digest('A'*64,'pin'),'a'*64)
        for v in ('',False,0,[],{},'g'*64,'a'*63):
            with self.subTest(value=v),self.assertRaises(ValueError):expected_digest(v,'pin')

    def test_complete_trusted_manifest_is_binding(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t).resolve();ref=p/'ref.fa';ref.write_text('>rna\nACGT\n');prefix=p/'idx'
            files=[p/('idx'+s) for s in ('.1.bt2','.2.bt2','.3.bt2','.4.bt2','.rev.1.bt2','.rev.2.bt2')]
            for f in files:f.write_text(f.name)
            parts={str(f)[len(str(prefix)):]:sha256(f) for f in files}
            data={'schema_version':1,'aligner':'bowtie2','reference_sha256':sha256(ref),'parts':parts}
            manifest=p/'manifest.json'
            def run(d):
                manifest.write_text(json.dumps(d))
                with patch('genome_common.index_files',return_value=files):return check_manifest(manifest,ref,prefix,'bowtie2')
            self.assertEqual(run(data),data)
            for change in ({'parts':{}},{'parts':{k:v for k,v in list(parts.items())[:1]}},
                           {'parts':dict(parts,**{'.1.bt2':None})},{'reference_sha256':None},
                           {'reference_sha256':'0'*64},{'schema_version':True},{'aligner':'hisat2'},
                           {'parts':dict(parts,**{'.1.bt2':'not-a-hash'})}):
                with self.subTest(change=change),self.assertRaises(ValueError):run(dict(data,**change))
            before={f:sha256(f) for f in [ref,*files]}
            bad=dict(data,parts=dict(parts,**{'.1.bt2':'0'*64}))
            with self.assertRaises(ValueError):run(bad)
            self.assertEqual(before,{f:sha256(f) for f in before})

    def test_standalone_settings_manifest_omitted_or_null(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t).resolve();cfg=p/'settings.json'
            for tool,space in [('rna-align','rna'),('genome-align','genome')]:
                for optional in ({},{'verification_manifest':None}):
                    data={'schema_version':2,'input':{'type':'trimmed_fastq'},
                          'references':{space:{'fasta':'ref.fa','index':{'action':'use_existing','prefix':'idx',**optional}}}}
                    if space=='genome':data['references'][space]['no_decoys']=True
                    cfg.write_text(json.dumps(data))
                    args=[tool,'--settings',str(cfg),'--input-file','reads.fq','--output-dir',str(p/'out')]
                    if space=='genome':args+=['--no-decoys']
                    plan=resolve(parser().parse_args(args))
                    self.assertIn(str(p/'ref.fa'),plan['command'])

    def test_trimming_precedence_and_zero_is_still_enabled(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t).resolve();cfg=p/'s.json'
            data={'schema_version':2,'input':{'type':'raw_fastq'},'trimming':{'adapter':'ACGT','nextseq_trim':None}}
            cfg.write_text(json.dumps(data))
            def plan(*extra):
                return resolve(parser().parse_args(['trim','--settings',str(cfg),'--input-file',str(p/'reads.fq'),'--output-dir',str(p/'out'),*extra]))['command']
            self.assertNotIn('--nextseq-trim',plan())
            self.assertIn('20',plan('--nextseq-trim','20'))
            data['trimming']['nextseq_trim']=20;cfg.write_text(json.dumps(data))
            self.assertNotIn('--nextseq-trim',plan('--no-nextseq-trim'))
            cmd=plan('--nextseq-trim','0');self.assertEqual(cmd[cmd.index('--nextseq-trim')+1],'0')
            for v in (False,'off',-1):
                data['trimming']['nextseq_trim']=v;cfg.write_text(json.dumps(data))
                with self.assertRaises(ValueError):plan()

    @unittest.skipUnless(shutil.which('cutadapt'),'Cutadapt required')
    def test_real_nextseq_off_keeps_high_quality_G_and_quality_trimming_stays_active(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t).resolve();reads=p/'reads.fq';seq='ACCTAACCTAACCTAACCTAA'+'G'*12
            reads.write_text('@high\n'+seq+'\n+\n'+'I'*len(seq)+'\n@low\n'+seq+'\n+\n'+'I'*(len(seq)-4)+'!!!!\n')
            for name,args in [('off',['--no-nextseq-trim']),('on',['--nextseq-trim','20']),('zero',['--nextseq-trim','0'])]:
                r=subprocess.run([sys.executable,str(ROOT/'toolbox.py'),'trim','--input-file',str(reads),'--adapter','TTTTTTTTTTTTTTTT','--quality','20','--minimum-length','1','--output-dir',str(p/name),*args],capture_output=True,text=True)
                self.assertEqual(r.returncode,0,r.stdout+r.stderr)
            def sequences(name):
                with gzip.open(p/name/'trimmed.fastq.gz','rt') as f:return f.read().splitlines()[1::4]
            self.assertEqual(sequences('off'),[seq,seq[:-4]])
            self.assertEqual(sequences('on')[0],seq[:-12])
            self.assertEqual(sequences('zero')[0],seq[:-12])
            methods=json.loads((p/'off/methods.json').read_text())
            self.assertNotIn('--nextseq-trim',methods['command'])
