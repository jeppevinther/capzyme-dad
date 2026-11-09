"""Standalone adapter tests; scientific regression remains in the existing suites."""
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'workflow/scripts'))
from core import sha256
from tool_settings import load, parser, resolve


class SettingsTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='toolbox '); self.addCleanup(self.tmp.cleanup)
        self.p = Path(self.tmp.name).resolve()
        self.cfg = self.p/'settings.json'
        self.data = {'schema_version': 2, 'samples': 'samples.tsv', 'output_dir': 'new-count',
                     'input': {'type': 'bam'}, 'references': {'rna': {'fasta': 'ref.fa'}},
                     'counting': {'minimum_mapq': 20, 'exclude_duplicates': True, 'rna_strands': 'plus'}}
        (self.p/'samples.tsv').write_text('sample_id\tinput_file\texperiment\ttreatment\nB\tb.bam\texp1\tMDE\nA\ta.bam\texp1\tRpp\n')

    def plan(self, *extra):
        self.cfg.write_text(json.dumps(self.data))
        return resolve(parser().parse_args(['count', '--settings', str(self.cfg), *extra]))

    def test_precedence_order_paths_boolean_override_and_defaults(self):
        plan = self.plan('--mapq', '0', '--no-exclude-duplicates', '--strand', 'minus')
        cmd = plan['command']
        self.assertEqual(cmd[cmd.index('--mapq')+1], '0')
        self.assertEqual(cmd[cmd.index('--strand')+1], 'minus')
        self.assertNotIn('--exclude-duplicates', cmd)
        self.assertEqual(cmd[-2:], [str(self.p/'b.bam'), str(self.p/'a.bam')])
        self.assertEqual(plan['output_dir'], str(self.p/'new-count'))
        self.assertEqual(plan['policies']['excluded_flags'], 2308)
        self.assertEqual([cmd[i+1] for i, v in enumerate(cmd) if v == '--sample-name'], ['B', 'A'])
        self.data.pop('counting')
        self.assertIn('both', self.plan()['command'])

    def test_unknown_and_unsupported_fields_fail_in_all_sections(self):
        base = copy.deepcopy(self.data)
        for section, key, value in [('execution', 'total_threads', 4), ('counting', 'banana', 1),
                                    ('mapping', 'genome', {'preset': 'fast'}), ('trimming', 'enabled', True),
                                    ('bigwig', 'genome', {'normalizations': ['raw', 'rpm']}),
                                    ('references', 'genome', {'index': {'action': 'build_new'}}),
                                    ('counting', 'minimum_mapq', True), ('counting', 'exclude_duplicates', 'false')]:
            self.data = copy.deepcopy(base); self.data.setdefault(section, {})[key] = value
            with self.subTest(section=section, key=key), self.assertRaises(ValueError): self.plan()
        self.data = base
        with self.assertRaises(ValueError): self.plan('--normalization', 'raw')
        with self.assertRaises(ValueError): self.plan('--mapq', '-1')
        self.data['schema_version'] = 1
        with self.assertRaises(ValueError): self.plan()

    def test_bigwig_resolution_and_invalid_combinations(self):
        self.data['input'] = {'type': 'counts', 'matrix': 'counts.tsv', 'assert_unthresholded': True}
        self.data['bigwig'] = {'rna': {'normalization': 'rpm', 'depths': 'depths.tsv', 'strands': 'plus'}}
        self.cfg.write_text(json.dumps(self.data))
        a = parser().parse_args(['bigwig', '--settings', str(self.cfg), '--normalization', 'raw'])
        with self.assertRaisesRegex(ValueError, 'depths'): resolve(a)
        self.data['bigwig']['rna'].pop('depths'); self.data['bigwig']['rna']['normalization'] = 'raw'
        self.cfg.write_text(json.dumps(self.data))
        plan = resolve(parser().parse_args(['bigwig', '--settings', str(self.cfg)]))
        self.assertEqual(plan['resolved_engine_parameters']['cohort'], ['B', 'A'])
        self.assertEqual(plan['resolved_engine_parameters']['strands'], 'plus')

    def test_json_help_and_preview_without_snakemake_or_yaml(self):
        env = dict(os.environ, PATH='/usr/bin:/bin', PYTHONDONTWRITEBYTECODE='1')
        self.cfg.write_text(json.dumps(self.data))
        for args in (['--help'], ['count', '--settings', str(self.cfg), '--preview']):
            r = subprocess.run([sys.executable, str(ROOT/'toolbox.py'), *args], env=env, cwd='/', text=True, capture_output=True)
            self.assertEqual(r.returncode, 0, r.stderr)
        self.assertFalse((self.p/'new-count').exists())

    def test_explicit_no_decoys_cli_is_a_positive_declaration(self):
        args=parser().parse_args(['genome-align','--input-file','reads.fq','--reference','ref.fa',
            '--index','index','--no-decoys','--output-dir','new-genome'])
        self.assertIs(args.no_decoys,True)
        plan=resolve(args)
        self.assertIn('--no-decoys',plan['command'])

    def test_duplicate_yaml_key_and_bad_sample_columns(self):
        self.cfg.write_text('{"schema_version":2,"counting":{"minimum_mapq":20,"minimum_mapq":0}}')
        with self.assertRaisesRegex(ValueError, 'Duplicate'): load(self.cfg)
        if importlib.util.find_spec('yaml'):
            self.cfg.write_text('schema_version: 2\ncounting:\n  minimum_mapq: 20\n  minimum_mapq: 0\n')
            with self.assertRaisesRegex(ValueError, 'Duplicate'): load(self.cfg)
        (self.p/'samples.tsv').write_text('sample_id\tpassword\nA\tsecret\n')
        with self.assertRaises(ValueError): self.plan()

    @unittest.skipUnless(shutil.which('cutadapt'), 'Cutadapt required')
    def test_trim_is_native_cutadapt_with_effective_settings(self):
        import gzip
        reads = self.p/'raw.fq'; seq = 'ACGTACGTACGTACGTACGTAGATCGGAAGAGCACACGTCT'
        reads.write_text('@x\n'+seq+'\n+\n'+'I'*len(seq)+'\n')
        self.data = {'schema_version': 2, 'input': {'type': 'raw_fastq'}, 'output_dir': 'trim',
                     'trimming': {'adapter': 'AGATCGGAAGAGCACACGTCT', 'quality': 20, 'minimum_length': 15, 'nextseq_trim': 20}}
        self.cfg.write_text(json.dumps(self.data))
        subprocess.run([sys.executable, str(ROOT/'toolbox.py'), 'trim', '--settings', str(self.cfg), '--input-file', str(reads)], check=True, capture_output=True)
        subprocess.run(['cutadapt', '-a', 'AGATCGGAAGAGCACACGTCT', '-q', '20', '-m', '15', '--nextseq-trim', '20', '--cores', '1', '-o', str(self.p/'native.fq.gz'), str(reads)], check=True, capture_output=True)
        self.assertEqual(gzip.open(self.p/'trim/trimmed.fastq.gz', 'rt').read(), gzip.open(self.p/'native.fq.gz', 'rt').read())
        self.assertEqual(json.loads((self.p/'trim/methods.json').read_text())['exit_status'], 0)

    @unittest.skipUnless(all(shutil.which(t) for t in ('hisat2', 'hisat2-build', 'hisat2-inspect', 'hisat2_extract_splice_sites.py', 'hisat2_extract_exons.py', 'bowtie2', 'bowtie2-build', 'samtools', 'bedtools')), 'Genome fixture dependencies required')
    def test_genome_settings_reuses_existing_spliced_engine(self):
        from test_genome_mapping import GenomeMappingTests, write_reads
        f = GenomeMappingTests
        f.setUpClass(); self.addCleanup(f.tearDownClass)
        reads = self.p/'reads.fq'; write_reads(reads, [('unique', f.chrom[140:180]), ('splice', f.chrom[200:220]+f.chrom[420:440])])
        subprocess.run([sys.executable, str(ROOT/'toolbox.py'), 'index-record', '--reference', str(f.fasta), '--index', f.ht, '--aligner', 'hisat2', '--output-dir', str(self.p/'record')], check=True, capture_output=True)
        self.data = {'schema_version': 2, 'input': {'type': 'trimmed_fastq'}, 'output_dir': 'genome',
                     'references': {'genome': {'fasta': str(f.fasta), 'decoys_file': str(f.decoys), 'index': {'action': 'use_existing', 'prefix': f.ht, 'verification_manifest': str(self.p/'record/reference_manifest.json')}}},
                     'mapping': {'genome': {'aligner': 'hisat2'}}}
        self.cfg.write_text(json.dumps(self.data))
        subprocess.run([sys.executable, str(ROOT/'toolbox.py'), 'genome-align', '--settings', str(self.cfg), '--input-file', str(reads)], check=True, capture_output=True)
        subprocess.run([sys.executable, str(ROOT/'workflow/scripts/hisat2_genome.py'), 'map', '--fastq', str(reads), '--reference', str(f.fasta), '--index', f.ht, '--decoys', str(f.decoys), '--output-dir', str(self.p/'native')], check=True, capture_output=True)
        self.assertEqual((self.p/'genome/qc.json').read_bytes(), (self.p/'native/qc.json').read_bytes())
        a = subprocess.check_output(['samtools', 'view', str(self.p/'genome/accepted.bam')])
        b = subprocess.check_output(['samtools', 'view', str(self.p/'native/accepted.bam')])
        self.assertEqual(a, b); self.assertIn(b'N', a)

    def test_collapse_is_unchanged_and_never_overwrites(self):
        ref = self.p/'ref.fa'; ref.write_text('>first\n'+'A'*30+'C\n>long\n'+'A'*30+'CC\n')
        before = sha256(ref)
        cli = [sys.executable, str(ROOT/'toolbox.py'), 'collapse', '--reference', str(ref), '--output-dir', str(self.p/'collapsed')]
        subprocess.run(cli, check=True, capture_output=True)
        subprocess.run([sys.executable, str(ROOT/'legacy_scripts/merge_by_prefix_v2.py'), str(ref), '-k', '30', '-o', str(self.p/'old')], check=True, capture_output=True)
        for suffix in ('fa', 'tsv'):
            self.assertEqual((self.p/'collapsed'/('collapsed.'+suffix)).read_bytes(), (self.p/('old.'+suffix)).read_bytes())
        self.assertNotEqual(subprocess.run(cli, capture_output=True).returncode, 0)
        self.assertEqual(sha256(ref), before)

    @unittest.skipUnless(shutil.which('samtools') and shutil.which('bedtools'), 'BAM tools required')
    def test_count_settings_exact_existing_cli_and_source_integrity(self):
        from test_core import CoordinateTests
        f = CoordinateTests('test_empty_matrix_and_total_threshold'); f.setUp(); self.addCleanup(f.doCleanups)
        bam = f.bam(); before = sha256(bam)
        self.data['references']['rna']['fasta'] = str(f.fasta)
        (self.p/'samples.tsv').write_text(f'sample_id\tinput_file\nB\t{bam}\n')
        self.cfg.write_text(json.dumps(self.data))
        subprocess.run([sys.executable, str(ROOT/'toolbox.py'), 'count', '--settings', str(self.cfg)], check=True, capture_output=True)
        subprocess.run([sys.executable, str(ROOT/'workflow/scripts/count_5ends.py'), '--reference', str(f.fasta),
                        '--sample-name', 'B', '--strand', 'plus', '--exclude-duplicates', '--output', str(self.p/'direct.tsv'), str(bam)], check=True, capture_output=True)
        self.assertEqual((self.p/'new-count/counts.tsv').read_bytes(), (self.p/'direct.tsv').read_bytes())
        self.assertEqual(sha256(bam), before)

    @unittest.skipUnless(shutil.which('bowtie2-build') and shutil.which('bowtie2-inspect') and shutil.which('samtools'), 'Bowtie2 required')
    def test_manifest_native_index_and_rna_mapping_existing_engine(self):
        ref = self.p/'ref.fa'; seq = 'ACGAGTAGCTAGGATCGATCGTACTGATCGTAGCTAGCATGCTAGCTACGATCGTACGATCGTAGCTA'
        ref.write_text('>rna\n'+seq+'\n')
        subprocess.run(['bowtie2-build', str(ref), str(self.p/'idx')], check=True, capture_output=True)
        before = {str(p): sha256(p) for p in self.p.glob('idx*')}
        record = [sys.executable, str(ROOT/'toolbox.py'), 'index-record', '--reference', str(ref), '--index', str(self.p/'idx'), '--aligner', 'bowtie2', '--output-dir', str(self.p/'record')]
        subprocess.run(record, check=True, capture_output=True)
        reads = self.p/'reads.fq'; reads.write_text('@x\n'+seq[:40]+'\n+\n'+'I'*40+'\n')
        self.data = {'schema_version': 2, 'input': {'type': 'trimmed_fastq'}, 'output_dir': 'mapped',
                     'references': {'rna': {'fasta': 'ref.fa', 'index': {'action': 'use_existing', 'prefix': 'idx', 'verification_manifest': 'record/reference_manifest.json'}}}}
        self.cfg.write_text(json.dumps(self.data))
        run = [sys.executable, str(ROOT/'toolbox.py'), 'rna-align', '--settings', str(self.cfg), '--input-file', str(reads)]
        subprocess.run(run, check=True, capture_output=True)
        from rna_handoff import main
        subprocess.run([sys.executable, str(ROOT/'workflow/scripts/rna_handoff.py'), '--fastq', str(reads), '--reference', str(ref), '--index', str(self.p/'idx'), '--norc', '--output-dir', str(self.p/'direct')], check=True, capture_output=True)
        for name in ('qc.json', 'unmapped.fastq'):
            self.assertEqual((self.p/'mapped'/name).read_bytes(), (self.p/'direct'/name).read_bytes())
        for p, h in before.items(): self.assertEqual(sha256(p), h)
        part = self.p/'idx.1.bt2'; part.write_bytes(part.read_bytes()+b'bad')
        result = subprocess.run([*run, '--preview'], text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0); self.assertIn('checksum', result.stderr)

    @unittest.skipUnless(importlib.util.find_spec('pyBigWig'), 'pyBigWig required')
    def test_settings_bigwig_is_byte_identical_to_direct_tool(self):
        (self.p/'sizes.tsv').write_text('rna\t10\n')
        (self.p/'matrix.tsv').write_text('chromosome\tstrand\tstart\tend\tB\tA\nrna\t+\t1\t1\t4\t8\nrna\t-\t10\t10\t2\t0\n')
        self.data = {'schema_version': 2, 'input': {'type': 'counts', 'matrix': 'matrix.tsv', 'assert_unthresholded': True},
                     'references': {'rna': {'chrom_sizes': 'sizes.tsv'}}, 'output_dir': 'tracks',
                     'bigwig': {'rna': {'normalization': 'raw', 'cohort': ['B', 'A'], 'strands': 'both'}}}
        self.cfg.write_text(json.dumps(self.data))
        subprocess.run([sys.executable, str(ROOT/'toolbox.py'), 'bigwig', '--settings', str(self.cfg)], check=True, capture_output=True)
        subprocess.run([sys.executable, str(ROOT/'workflow/scripts/bigwig2.py'), '--matrix', str(self.p/'matrix.tsv'), '--assert-unthresholded', '--chrom-sizes', str(self.p/'sizes.tsv'), '--cohort', 'B', 'A', '--normalization', 'raw', '--output-dir', str(self.p/'direct')], check=True, capture_output=True)
        for p in (self.p/'tracks/tracks').glob('*.bw'):
            self.assertEqual(p.read_bytes(), (self.p/'direct/tracks'/p.name).read_bytes())


if __name__ == '__main__': unittest.main()
