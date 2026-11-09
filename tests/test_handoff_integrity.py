"""Final accounting must bind bytes across stages, not just read identifiers."""
import copy
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'workflow/scripts'))
from core import sha256
import sequential


@unittest.skipUnless(shutil.which('samtools'), 'SAMtools required for synthetic BAM accounting')
class HandoffIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='handoff integrity ')
        self.addCleanup(self.temp.cleanup)
        self.p = Path(self.temp.name).resolve()
        self.rna, self.genome = self.p/'rna', self.p/'genome'
        for d in (self.rna, self.genome): d.mkdir()
        header = '@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:ref\tLN:100\n'
        accepted = 'rna\t0\tref\t1\t60\t30M\t*\t0\t0\t'+'A'*30+'\t'+'I'*30+'\n'
        unmapped = 'g\t4\t*\t0\t0\t*\t*\t0\t0\t'+'C'*30+'\t'+'I'*30+'\n'
        mapped = 'g\t0\tref\t20\t60\t30M\t*\t0\t0\t'+'C'*30+'\t'+'I'*30+'\n'
        for d, all_records, good in [(self.rna, accepted+unmapped, accepted), (self.genome, mapped, mapped)]:
            for name, body in [('alignments', all_records), ('accepted', good)]:
                sam = d/(name+'.sam'); sam.write_text(header+body)
                subprocess.run(['samtools','view','-b','-o',str(d/(name+'.bam')),str(sam)],check=True,capture_output=True)
        self.handoff = self.rna/'unmapped.fastq'
        self.handoff.write_text('@g\n'+'C'*30+'\n+\n'+'I'*30+'\n')
        (self.genome/'unmapped.fastq').write_text('')
        for d, outcomes in [(self.rna, {'accepted':1,'unmapped':1}), (self.genome, {'accepted':1})]:
            qc = {'primary_outcomes':outcomes,'mapq_minimum':20,'excluded_flags':2308}
            (d/'qc.json').write_text(json.dumps(qc))
            meta = {'status':'complete','inputs_unchanged':True,'input_reads':sum(outcomes.values()),
                    'qc':qc,'parameters':{'output_dir':str(d)},
                    'outputs_sha256':{n:sha256(d/n) for n in ['alignments.bam','accepted.bam','qc.json','unmapped.fastq']}}
            if d == self.genome:
                meta['parameters']['fastq'] = str(self.handoff)
                meta['inputs_sha256'] = {str(self.handoff):sha256(self.handoff)}
            (d/'provenance.json').write_text(json.dumps(meta))

    def meta(self, stage): return json.loads((stage/'provenance.json').read_text())
    def save(self, stage, meta): (stage/'provenance.json').write_text(json.dumps(meta))
    def run_accounting(self, suffix='result'):
        with tempfile.TemporaryFile(mode='w+') as log:
            return sequential.accounting(self.rna, self.genome, self.p/(suffix+'.json'), log)

    def test_unchanged_valid_handoff_and_consumed_hashes(self):
        result = self.run_accounting()
        self.assertTrue(result['handoff_sequence_quality_verified'])
        chain = result['handoff_integrity']
        self.assertEqual({chain[k] for k in ['rna_output_sha256','genome_input_sha256','current_handoff_sha256']}, {sha256(self.handoff)})
        for d in [self.rna,self.genome]:
            for n in ['alignments.bam','accepted.bam','qc.json']:
                self.assertEqual(result['verified_artifacts_sha256'][str(d/n)],sha256(d/n))

    def test_same_id_changed_sequence_or_quality_is_rejected_even_if_genome_pin_updated(self):
        original = self.handoff.read_text()
        gm = self.meta(self.genome)
        for label, replacement in [('sequence',original.replace('C'*30,'G'*30)),
                                   ('quality',original.replace('I'*30,'H'*30))]:
            for update_genome in [False, True]:
                with self.subTest(change=label, genome_updated=update_genome):
                    self.handoff.write_text(replacement)
                    changed = copy.deepcopy(gm)
                    if update_genome: changed['inputs_sha256'][str(self.handoff)] = sha256(self.handoff)
                    self.save(self.genome, changed)
                    with self.assertRaisesRegex(ValueError,'checksum'):
                        self.run_accounting()
                    self.assertFalse((self.p/'result.json').exists())
        # Restoring RNA bytes does not repair a genome input recorded from other bytes.
        self.handoff.write_text(original)
        with self.assertRaisesRegex(ValueError,'checksum'):
            self.run_accounting()

    def test_missing_incomplete_and_inconsistent_provenance_fails_closed(self):
        original = self.meta(self.genome)
        cases = [dict(original, status='failed'), dict(original, inputs_unchanged=False),
                 dict(original, inputs_unchanged='true'), dict(original, outputs_sha256={}),
                 dict(original, inputs_sha256={}), dict(original, qc={}), dict(original, input_reads=2)]
        other = copy.deepcopy(original); other['parameters']['fastq'] = str(self.p/'unrelated/unmapped.fastq'); cases.append(other)
        other = copy.deepcopy(original); other['inputs_sha256'][str(self.handoff)] = 'invalid'; cases.append(other)
        for meta in cases:
            with self.subTest(meta=meta):
                self.save(self.genome,meta)
                with self.assertRaises(ValueError): self.run_accounting()
        (self.genome/'provenance.json').unlink()
        with self.assertRaisesRegex(ValueError,'integrity record'): self.run_accounting()

    def test_missing_or_modified_consumed_bam_and_qc_artifacts(self):
        for d in (self.rna,self.genome):
            for name in ['alignments.bam','accepted.bam','qc.json']:
                p=d/name; original=p.read_bytes()
                with self.subTest(stage=d.name,artifact=name):
                    p.write_bytes(original+b'changed')
                    with self.assertRaisesRegex(ValueError,'artifact checksum'):self.run_accounting()
                    p.unlink()
                    with self.assertRaisesRegex(ValueError,'artifact checksum'):self.run_accounting()
                    p.write_bytes(original)

    def test_published_attempt_paths_and_new_attempt_restart(self):
        for attempt in ['first','restart']:
            original_rna = self.p/'.attempts'/attempt/'rna/result'
            for d in [self.rna,self.genome]:
                meta = self.meta(d)
                meta['parameters']['output_dir'] = str(original_rna if d == self.rna else self.p/'.attempts'/attempt/'genome/result')
                if d == self.genome:
                    source = str(original_rna/'unmapped.fastq')
                    meta['parameters']['fastq'] = source
                    meta['inputs_sha256'] = {source:sha256(self.handoff)}
                self.save(d,meta)
                (d/'publication.json').write_text(json.dumps({
                    'standalone_output_dir':meta['parameters']['output_dir'],
                    'published_output_dir':str(d),'standalone_provenance_sha256':sha256(d/'provenance.json')}))
            result = self.run_accounting(attempt)
            self.assertTrue(result['handoff_sequence_quality_verified'])
        receipt = json.loads((self.rna/'publication.json').read_text())
        receipt['standalone_provenance_sha256'] = '0'*64
        (self.rna/'publication.json').write_text(json.dumps(receipt))
        with self.assertRaisesRegex(ValueError,'Publication'): self.run_accounting()
        (self.rna/'publication.json').unlink()
        with self.assertRaisesRegex(ValueError,'integrity record'): self.run_accounting()

    def test_artifacts_changed_during_accounting_are_rejected(self):
        original = sequential.bam_records
        def changed(path, log):
            yield from original(path,log)
            if Path(path) == self.genome/'accepted.bam':
                with (self.rna/'qc.json').open('a') as handle: handle.write(' ')
        with patch('sequential.bam_records',side_effect=changed):
            with self.assertRaisesRegex(ValueError,'changed during'): self.run_accounting()
        self.assertFalse((self.p/'result.json').exists())


if __name__ == '__main__': unittest.main()
