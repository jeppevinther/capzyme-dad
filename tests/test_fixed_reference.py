"""Fixed-reference selection, hash rejection and immutable-input regressions."""
import copy
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'workflow/scripts'))
from configuration import load_settings
from core import collapse, combine, sha256
from fixed_reference import PARTS, check_fixed_hashes, verify_fixed
from synthetic import generate

AVAILABLE = all(shutil.which(t) for t in ('snakemake', 'bowtie2', 'bowtie2-build',
                                         'bowtie2-inspect', 'samtools', 'bedtools', 'cutadapt'))


@unittest.skipUnless(AVAILABLE, 'Full local workflow and index-inspection tools required')
class FixedReferenceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='fixed reference test ')
        self.addCleanup(self.tmp.cleanup)
        self.p = Path(self.tmp.name).resolve()
        self.config_path = generate(self.p)
        self.config = json.loads(self.config_path.read_text())
        self.inputs = self.p / 'immutable inputs'
        self.inputs.mkdir()
        collapse(self.p / 'ncrna.fa', self.p / 'merged', 30, [], ROOT / 'legacy_scripts/merge_by_prefix_v2.py')
        self.fasta = self.inputs / 'historical.fa'
        combine([self.p / 'viral.fa', self.p / 'merged.fa'], self.fasta)
        self.fasta.write_text(self.fasta.read_text().replace('RNU2__n2', 'FIXED_ONLY'))
        self.prefix = self.inputs / 'historical'
        subprocess.run(['bowtie2-build', str(self.fasta), str(self.prefix)],
                       check=True, capture_output=True)
        self.fixed = {'mode': 'fixed', 'fasta': str(self.fasta),
                      'fasta_sha256': sha256(self.fasta), 'index_prefix': str(self.prefix),
                      'index_format': 'bt2',
                      'index_sha256': {part + '.bt2': sha256(str(self.prefix) + part + '.bt2') for part in PARTS}}
        self.config['reference'] = self.fixed
        self.config['mapping']['keep_diagnostic_sam'] = True
        self.config_path.write_text(json.dumps(self.config))

    def cmd(self, dry=False, force=False):
        command = ['snakemake', '--snakefile', str(ROOT / 'workflow/Snakefile'),
                   '--configfile', str(self.config_path), '--directory', str(self.p), '--cores', '4']
        if dry:
            command += ['--dry-run']
        if force:
            command += ['--forceall']
        return command

    def snapshot(self):
        return {p.name: (sha256(p), p.stat().st_size, p.stat().st_mtime_ns, p.stat().st_mode)
                for p in self.inputs.iterdir()}

    def test_fixed_selection_no_preparation_and_read_only_even_when_forced(self):
        for p in self.inputs.iterdir():
            p.chmod(0o444)
        self.inputs.chmod(0o555)
        def restore():
            self.inputs.chmod(0o755)
            for p in self.inputs.iterdir():
                p.chmod(0o644)
        self.addCleanup(restore)
        before = self.snapshot()
        dry = subprocess.run(self.cmd(dry=True), text=True, capture_output=True)
        self.assertEqual(dry.returncode, 0, dry.stdout + dry.stderr)
        for excluded in ('collapse_ncrna', 'combine_reference', 'bowtie2_index', 'index_reference'):
            self.assertNotIn(excluded, dry.stdout + dry.stderr)
        self.assertIn('verify_fixed_reference', dry.stdout + dry.stderr)
        run = subprocess.run(self.cmd(force=True), text=True, capture_output=True)
        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
        self.assertEqual(self.snapshot(), before)
        self.assertFalse(Path(str(self.fasta) + '.fai').exists())
        observed = (self.p / 'run/counts/5ends.raw.tsv').read_text()
        expected = (self.p / 'expected.tsv').read_text().replace('RNU2__n2', 'FIXED_ONLY')
        self.assertEqual(observed, expected)
        self.assertTrue((self.p / 'run/diagnostics/B.unfiltered.sam').is_file())
        report = json.loads((self.p / 'run/reference/fixed_verified.json').read_text())
        self.assertTrue(report['index_matches_fasta'])
        self.assertEqual(report['sequence_count'], 3)
        self.assertEqual(report['index_format'], 'bt2')
        provenance = json.loads((self.p / 'run/provenance/run.json').read_text())
        self.assertEqual(provenance['settings']['reference_mode'], 'fixed')
        self.assertEqual(provenance['inputs'][str(self.fasta)], self.fixed['fasta_sha256'])
        for part in PARTS:
            self.assertIn(str(self.prefix) + part + '.bt2', provenance['inputs'])

    def test_wrong_fasta_checksum_rejected_before_dag_execution(self):
        self.config['reference']['fasta_sha256'] = '0' * 64
        self.config_path.write_text(json.dumps(self.config))
        result = subprocess.run(self.cmd(dry=True), text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('SHA-256 mismatch', result.stdout + result.stderr)
        self.assertFalse((self.p / 'run').exists())

    def test_wrong_index_checksum_rejected(self):
        bad = copy.deepcopy(self.fixed)
        bad['index_sha256']['.1.bt2'] = '0' * 64
        with self.assertRaisesRegex(ValueError, 'SHA-256 mismatch'):
            check_fixed_hashes(bad)

    def test_malformed_pins_and_mixed_modes_rejected(self):
        for field, value in [('fasta_sha256', 'bad'), ('index_sha256', {'.1.bt2': None}), ('index_format', 'invalid')]:
            cfg = copy.deepcopy(self.config)
            cfg['reference'][field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                load_settings(cfg, ROOT)
        cfg = copy.deepcopy(self.config)
        cfg['reference']['ncrna_fasta'] = str(self.p / 'ncrna.fa')
        with self.assertRaisesRegex(ValueError, 'Unknown reference fields'):
            load_settings(cfg, ROOT)
        missing = copy.deepcopy(self.fixed)
        missing['index_prefix'] = str(self.inputs / 'absent')
        with self.assertRaisesRegex(ValueError, 'Missing fixed-reference input'):
            check_fixed_hashes(missing)

    def test_omitted_null_and_partial_expectations_record_all_inputs(self):
        for mode in ('omitted','null','partial','empty'):
            ref=dict(self.fixed)
            ref.pop('fasta_sha256');ref.pop('index_sha256')
            if mode=='null':ref.update(fasta_sha256=None,index_sha256=None)
            if mode=='partial':ref['index_sha256']={'.1.bt2':self.fixed['index_sha256']['.1.bt2'].upper()}
            if mode=='empty':ref['index_sha256']={}
            before=self.snapshot()
            report=self.p/(mode+'.json');verify_fixed(ref,self.p/(mode+'.sizes'),report)
            got=json.loads(report.read_text())
            self.assertEqual(len(got['inputs_sha256']),7)
            self.assertTrue(got['index_matches_fasta'])
            self.assertEqual(self.snapshot(),before)
        bad=dict(ref,index_sha256={'.1.bt2':'0'*64})
        with self.assertRaisesRegex(ValueError,'mismatch'):check_fixed_hashes(bad)
        for value in ([1],{'.unknown':'0'*64},{'.1.bt2':False}):
            with self.assertRaises(ValueError):check_fixed_hashes(dict(ref,index_sha256=value))

    def test_content_mismatch_blocks_reads_even_with_valid_file_pins(self):
        text = self.fasta.read_text()
        lines = text.splitlines()
        lines[1] = ('A' if lines[1][0] != 'A' else 'C') + lines[1][1:]
        self.fasta.write_text('\n'.join(lines) + '\n')
        self.config['reference']['fasta_sha256'] = sha256(self.fasta)
        self.config_path.write_text(json.dumps(self.config))
        before = self.snapshot()
        run = subprocess.run(self.cmd(), text=True, capture_output=True)
        self.assertNotEqual(run.returncode, 0)
        self.assertIn('does not match the fixed FASTA', run.stdout + run.stderr)
        self.assertFalse((self.p / 'run/trimmed/B.fastq.gz').exists())
        self.assertEqual(self.snapshot(), before)

    def test_large_index_mode_and_ambiguous_prefix(self):
        prefix = self.inputs / 'large'
        subprocess.run(['bowtie2-build', '--large-index', str(self.fasta), str(prefix)],
                       check=True, capture_output=True)
        large = dict(self.fixed, index_prefix=str(prefix), index_format='bt2l',
                     index_sha256={part + '.bt2l': sha256(str(prefix) + part + '.bt2l') for part in PARTS})
        out = self.p / 'verified.json'
        verify_fixed(large, self.p / 'sizes.tsv', out)
        self.assertTrue(json.loads(out.read_text())['index_matches_fasta'])
        # A mixed-family prefix cannot reliably select a specific family in Bowtie2.
        shutil.copyfile(str(self.prefix) + '.1.bt2', str(prefix) + '.1.bt2')
        with self.assertRaisesRegex(ValueError, 'Ambiguous index prefix'):
            check_fixed_hashes(large)


if __name__ == '__main__':
    unittest.main()
