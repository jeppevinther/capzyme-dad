"""Synthetic regression of standalone, BAM-only, FASTQ and historical counters."""
import copy
import csv
import io
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
from bam_count import reference_sizes, prepare_bam
from bam_configuration import load_bam_settings
from compare_matrices import compare
from core import sha256
from synthetic import generate
import test_core as fixtures

CLI = ROOT/'workflow/scripts/count_5ends.py'
SNAKE = ROOT/'workflow/existing_bam.smk'
TOOLS = all(shutil.which(t) for t in ('samtools', 'bedtools', 'snakemake'))


def run(args):
    result = subprocess.run(list(map(str, args)), text=True, capture_output=True)
    if result.returncode:
        raise AssertionError(result.stdout + result.stderr)
    return result.stdout + result.stderr


@unittest.skipUnless(TOOLS, 'samtools, bedtools and Snakemake required')
class ExistingBamTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.CoordinateTests('test_empty_matrix_and_total_threshold')
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.p, self.fasta = self.fixture.p, self.fixture.fasta
        self.b = self.p/'B.bam'
        self.empty = self.p/'empty.bam'
        # Deliberately unsorted and unfiltered BAM; includes every boundary/flag fixture.
        run(['samtools', 'view', '-b', '-o', self.b, self.fixture.sam])
        shutil.copyfile(self.fixture.bam(empty=True), self.empty)
        self.sources = [self.b, self.empty]
        for source in [*self.sources, self.fasta]:
            source.chmod(0o444)
        self.before = {p: (sha256(p), p.stat().st_mtime_ns) for p in [*self.sources, self.fasta]}
        self.config = {'schema_version': 1, 'input_mode': 'existing_bam',
                       'samples': str(self.p/'samples.tsv'), 'output_dir': str(self.p/'run'),
                       'reference': {'fasta': str(self.fasta), 'fasta_sha256': sha256(self.fasta)},
                       'counting': {'mapq': 20, 'strand': 'both', 'exclude_duplicates': False, 'min_total_count': 0}}
        (self.p/'samples.tsv').write_text('sample_id\tbam\n'+''.join(f'{p.stem}\t{p}\n' for p in self.sources))

    def assert_unchanged(self):
        for p, state in self.before.items():
            self.assertEqual((sha256(p), p.stat().st_mtime_ns), state)
        self.assertFalse(Path(str(self.b)+'.bai').exists())

    def workflow(self, config, dry=False):
        cfg = self.p/'config.json'
        cfg.write_text(json.dumps(config))
        return run(['snakemake', '--snakefile', SNAKE, '--configfile', cfg,
                    '--directory', self.p, '--cores', '4', '--resources', 'mem_mb=16384',
                    *(['--dry-run'] if dry else [])])

    def standalone(self, output, mapq=20, strand='both', minimum=0, duplicates=False):
        run([sys.executable, CLI, '-g', self.fasta, '-o', output,
             '--sample-name', 'B', '--sample-name', 'empty', '--mapq', str(mapq),
             '--strand', strand, '--min-total-count', str(minimum),
             *(['--exclude-duplicates'] if duplicates else []), *self.sources])

    def test_three_way_legacy_cli_workflow_and_readonly(self):
        bash = os.environ.get('CAPZYME_TEST_BASH', 'bash')
        if int(subprocess.check_output([bash, '-c', 'echo ${BASH_VERSINFO[0]}'], text=True).strip()) < 4:
            self.skipTest('Legacy counter requires Bash >=4')
        for n, (mapq, strand, minimum) in enumerate([(20, 'both', 0), (0, 'plus', 2), (30, 'minus', 1)]):
            cfg = copy.deepcopy(self.config)
            cfg['output_dir'] = str(self.p/f'run{n}')
            cfg['counting'].update(mapq=mapq, strand=strand, min_total_count=minimum)
            cli, legacy = self.p/f'cli{n}.tsv', self.p/f'legacy{n}.tsv'
            self.standalone(cli, mapq, strand, minimum)
            run([bash, ROOT/'legacy_scripts/count5ends_matrix_v1_3.sh', '-g', self.fasta,
                 '--tmp-dir', self.p, '--mapq', str(mapq), '--strand', strand,
                 '--min-count', str(minimum), '-o', legacy, *self.sources])
            dry = self.workflow(cfg, dry=True)
            self.assertIn('prepare_existing_bam', dry)
            for forbidden in ('map_filter_sort', 'trim_reads', 'collapse_ncrna', 'build_bowtie2'):
                self.assertNotIn(forbidden, dry)
            self.workflow(cfg)
            observed = Path(cfg['output_dir'])/f'counts/5ends.min_total_{minimum}.tsv'
            sizes = Path(cfg['output_dir'])/'reference/sizes.tsv'
            for target in (cli, legacy):
                result = compare(target, observed, sizes, minimum)
                self.assertTrue(result['equal'], result)
                self.assertEqual(result['unequal_cells'], 0)
            self.assert_unchanged()
            self.assertIn('Nothing to be done', self.workflow(cfg, dry=True))

    def test_duplicate_policy_already_filtered_sorted_and_single_bam(self):
        out = self.p/'deduplicated.tsv'
        self.standalone(out, duplicates=True)
        with out.open() as fh:
            rows = list(csv.reader(fh, delimiter='\t'))
        plus1 = [r for r in rows[1:] if r[1:3] == ['+', '1']]
        self.assertEqual(plus1[0][4:], ['1', '0'])
        self.assertTrue(any(r[1:3] == ['+', '6'] for r in rows[1:]))  # QC-fail retained
        single = self.p/'single.tsv'
        accepted = self.fixture.bam(flags=3332)
        run([sys.executable, CLI, '-g', self.fasta, '-o', single,
             '--sample-name', 'B', '--exclude-duplicates', accepted])
        a = list(csv.reader(io.StringIO(out.read_text()), delimiter='\t'))
        b = list(csv.reader(io.StringIO(single.read_text()), delimiter='\t'))
        self.assertEqual([r[:5] for r in a], b)
        cfg = copy.deepcopy(self.config)
        cfg['counting']['exclude_duplicates'] = True
        self.workflow(cfg)
        self.assertEqual(out.read_bytes(), (self.p/'run/counts/5ends.min_total_0.tsv').read_bytes())
        self.assert_unchanged()

    def test_sort_only_unordered_output_copies_and_no_indexes(self):
        sorted_bam = self.fixture.bam()
        sizes = self.p/'sizes.tsv'
        reference_sizes(self.fasta, sizes)
        for name, source, should_sort in [('unordered', self.b, True), ('ordered', sorted_bam, False)]:
            before = (sha256(source), source.stat().st_mtime_ns)
            log_path = self.p/(name+'.jsonl')
            with log_path.open('w') as log:
                result = prepare_bam(source, self.p/(name+'.accepted.bam'), sizes, log=log)
            self.assertEqual(result['sort_performed'], should_sort)
            commands = [json.loads(line)['command'] for line in log_path.read_text().splitlines()
                        if line.startswith('{') and '"command"' in line]
            self.assertEqual(any(c[:2] == ['samtools', 'sort'] for c in commands), should_sort)
            self.assertFalse(any(c[:2] == ['samtools', 'index'] for c in commands))
            self.assertEqual((sha256(source), source.stat().st_mtime_ns), before)
            self.assertFalse(Path(str(source)+'.bai').exists())
        self.assert_unchanged()

    def test_invalid_inputs_and_checksums_rejected(self):
        cfg = copy.deepcopy(self.config)
        cfg['reference']['fasta_sha256'] = '0'*64
        with self.assertRaisesRegex(ValueError, 'SHA-256'):
            load_bam_settings(cfg, ROOT)
        cfg = copy.deepcopy(self.config)
        cfg['counting']['mapq'] = True
        with self.assertRaises(ValueError):
            load_bam_settings(cfg, ROOT)
        badref = self.p/'bad.fa'
        badref.write_text('>other\nAAAA\n')
        result = subprocess.run([sys.executable, str(CLI), '-g', str(badref), '-o', str(self.p/'bad.tsv'), str(self.b)], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('dictionary mismatch', result.stderr)
        # Even explicit output collisions must never truncate an input.
        result = subprocess.run([sys.executable, str(CLI), '-g', str(self.fasta), '-o', str(self.b), str(self.b)], capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assert_unchanged()

    def test_exact_comparator_and_workflow_comparison(self):
        target = self.p/'target.tsv'
        self.standalone(target)
        cfg = copy.deepcopy(self.config)
        cfg['comparison'] = {'expected_matrix': str(target), 'expected_matrix_sha256': sha256(target),
                             'require_column_order': True, 'integer_tolerance': 0, 'expected_features': 8}
        self.workflow(cfg)
        report = json.loads((self.p/'run/comparison/exact.json').read_text())
        self.assertTrue(report['equal'])
        sizes = self.p/'run/reference/sizes.tsv'
        lines = target.read_text().splitlines()
        changed = self.p/'changed.tsv'
        row = lines[1].split('\t'); row[4] = str(int(row[4])+1)
        changed.write_text('\n'.join([lines[0], '\t'.join(row), *lines[2:]])+'\n')
        mismatch = compare(target, changed, sizes)
        self.assertFalse(mismatch['equal'])
        self.assertEqual(mismatch['unequal_cells'], 1)
        self.assertFalse(mismatch['sample_totals_equal'])
        reordered = self.p/'reordered.tsv'
        reordered.write_text('\n'.join([lines[0], *reversed(lines[1:])])+'\n')
        self.assertTrue(compare(target, reordered, sizes)['equal'])
        self.assertFalse(compare(target, reordered, sizes)['byte_equal'])
        changed.write_text('\n'.join([lines[0], *lines[2:]])+'\n')
        self.assertEqual(compare(target, changed, sizes)['missing_features'], 1)
        reordered_columns = self.p/'reordered_columns.tsv'
        rows = [row.split('\t') for row in lines]
        reordered_columns.write_text(''.join('\t'.join([*r[:4], r[5], r[4]])+'\n' for r in rows))
        self.assertFalse(compare(target, reordered_columns, sizes)['sample_order_equal'])
        self.assertFalse(compare(target, reordered_columns, sizes)['equal'])
        changed.write_text(lines[0]+'\n'+lines[1]+'\n'+lines[1]+'\n')
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            compare(target, changed, sizes)
        changed.write_text(lines[0]+'\n'+lines[1].replace('\t1\t1\t', '\t0\t0\t')+'\n')
        with self.assertRaisesRegex(ValueError, 'coordinate'):
            compare(target, changed, sizes)
        self.assert_unchanged()


@unittest.skipUnless(TOOLS and all(shutil.which(t) for t in ('cutadapt', 'bowtie2')), 'Full workflow tools required')
class FastqCompatibilityTests(unittest.TestCase):
    def test_existing_fastq_workflow_matches_standalone(self):
        with tempfile.TemporaryDirectory(prefix='capzyme-fastq-compat-') as tmp:
            cfg = generate(tmp)
            run(['snakemake', '--snakefile', ROOT/'workflow/Snakefile', '--configfile', cfg,
                 '--directory', tmp, '--cores', '4'])
            p = Path(tmp)
            standalone = p/'standalone.tsv'
            args = [sys.executable, CLI, '-g', p/'run/reference/combined.fa', '-o', standalone]
            for name in ['B', 'A', 'empty']:
                args += ['--sample-name', name]
            args += [p/f'run/bam/{name}.accepted.bam' for name in ['B', 'A', 'empty']]
            run(args)
            self.assertEqual(standalone.read_bytes(), (p/'run/counts/5ends.raw.tsv').read_bytes())
