import copy
import csv
import importlib.util
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
sys.path.insert(0, str(ROOT / 'workflow/scripts'))
from core import collapse, combine, matrix, threshold
from configuration import load_settings
from synthetic import generate

LEGACY = ROOT / 'legacy_scripts/merge_by_prefix_v2.py'
TOOLS = all(shutil.which(x) for x in ('samtools', 'bedtools'))


class ReferenceTests(unittest.TestCase):
    def test_k30_longest_and_input_order_ties(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            fasta = p / 'input.fa'
            seq = 'A' * 30
            fasta.write_text(f'>short gene_symbol:RNU2-1\n{seq}C\n'
                             f'>first_long gene_symbol:RNU2-2\n{seq}CC\n'
                             f'>second_long gene_symbol:RNU2-3\n{seq}GG\n'
                             '>short1\nacg\n>short2\nACG\n'
                             '>pseudo description: a pseudogene\nTTTT\n>empty_record\n')
            stats = collapse(fasta, p / 'new', 30, [], LEGACY)
            rows = list(csv.DictReader(io.StringIO((p / 'new.tsv').read_text()), delimiter='\t'))
            self.assertEqual(stats['groups'], 2)
            self.assertEqual(stats['skipped'], {'empty_sequence': 1})
            self.assertEqual(rows[0]['representative_id'], 'first_long')
            self.assertEqual(rows[0]['member_ids'], 'short,first_long,second_long')
            self.assertEqual(rows[0]['merged_id'], 'RNU2__n3')
            self.assertIn(seq + 'CC', (p / 'new.fa').read_text())
            self.assertEqual(rows[1]['member_count'], '2')
            subprocess.run([sys.executable, str(LEGACY), str(fasta), '-k', '30', '-o', str(p / 'old')],
                           check=True, capture_output=True)
            for suffix in ('.fa', '.tsv'):
                self.assertEqual((p / ('new' + suffix)).read_bytes(), (p / ('old' + suffix)).read_bytes())
            self.assertEqual(collapse(fasta, p / 'k70', 70, [], LEGACY)['groups'], 4)

    def test_named_biotype_filter_and_id_collisions(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            fasta = p / 'input.fa'
            fasta.write_text('>a gene_symbol:SAME\nAAAA\n>b gene_symbol:SAME\nCCCC\n'
                             '>drop gene_biotype:lncRNA\nGGGG\n')
            collapse(fasta, p / 'out', 30, ['lncRNA'], LEGACY)
            rows = list(csv.DictReader(io.StringIO((p / 'out.tsv').read_text()), delimiter='\t'))
            self.assertEqual([r['merged_id'] for r in rows], ['SAME__n1', 'SAME__n1__g2'])

    def test_duplicate_and_invalid_combined_references_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            one, two = p / 'one.fa', p / 'two.fa'
            one.write_text('>rna\nACGT\n'); two.write_text('>rna\nTTTT\n')
            with self.assertRaisesRegex(ValueError, 'Duplicate'):
                combine([one, two], p / 'out.fa')
            two.write_text('>other\nACGX\n')
            with self.assertRaisesRegex(ValueError, 'Invalid'):
                combine([two], p / 'out.fa')


class ConfigurationTests(unittest.TestCase):
    def test_manifest_order_and_bad_configurations(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = generate(tmp)
            cfg = json.loads(path.read_text())
            settings = load_settings(cfg, ROOT)
            self.assertEqual(list(settings['samples']), ['B', 'A', 'empty'])
            self.assertEqual(settings['flags'], 2308)
            default_k = copy.deepcopy(cfg)
            del default_k['reference']['prefix_length']
            self.assertEqual(load_settings(default_k, ROOT)['k'], 30)
            cases = []
            for section, field, value in [('counting', 'mapq', -1), ('counting', 'mapq', True),
                                           ('counting', 'strand', 'unknown'), ('reference', 'prefix_length', 0),
                                           ('mapping', 'norc', 'true'), ('preprocessing', 'quality', None)]:
                bad = copy.deepcopy(cfg); bad[section][field] = value; cases.append(bad)
            bad = copy.deepcopy(cfg); bad['output_dir'] = str(ROOT / 'results'); cases.append(bad)
            bad = copy.deepcopy(cfg); bad['output_dir'] = tmp; cases.append(bad)
            bad = copy.deepcopy(cfg); bad['unknown'] = 1; cases.append(bad)
            for case in cases:
                with self.subTest(config=case), self.assertRaises(ValueError):
                    load_settings(case, ROOT)
            manifest = Path(cfg['samples'])
            original_manifest = manifest.read_text()
            manifest.write_text(original_manifest + original_manifest.splitlines()[1] + '\n')
            with self.assertRaisesRegex(ValueError, 'duplicate'):
                load_settings(cfg, ROOT)
            for bad_manifest in (
                original_manifest.replace('B\t', 'start\t', 1),
                original_manifest.replace('B\t', '../bad\t', 1),
                'sample_id\tfastq\tread2\n',
                'sample_id\tfastq\n',
                original_manifest.replace('/A.fastq', '/B.fastq'),
            ):
                manifest.write_text(bad_manifest)
                with self.subTest(manifest=bad_manifest), self.assertRaises(ValueError):
                    load_settings(cfg, ROOT)


@unittest.skipUnless(TOOLS, 'samtools and bedtools required')
class CoordinateTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.p = Path(self.tmp.name)
        self.fasta = self.p / 'reference.fa'
        self.fasta.write_text('>rna\n' + 'A' * 100 + '\n')
        subprocess.run(['samtools', 'faidx', str(self.fasta)], check=True)
        # (name, flag, position1, MAPQ, CIGAR, query length).
        self.alignments = [('first', 0, 1, 20, '10M', 10),
                           ('duplicate_kept', 1024, 1, 60, '10M', 10),
                           ('low', 0, 2, 19, '10M', 10),
                           ('zero_mapq', 0, 3, 0, '10M', 10),
                           ('reverse_boundary', 16, 91, 20, '10M', 10),
                           ('reverse_soft', 16, 50, 60, '2S8M', 10),
                           ('reverse_splice', 16, 20, 60, '5M10N5M', 10),
                           ('reverse_deletion', 16, 60, 60, '5M2D5M', 10),
                           ('insertion', 0, 80, 60, '5M2I5M', 12),
                           ('hard_clip', 0, 40, 60, '2H10M', 10),
                           ('secondary', 256, 4, 60, '10M', 10),
                           ('supplementary', 2048, 5, 60, '10M', 10),
                           ('qcfail_kept', 512, 6, 20, '10M', 10)]
        self.sam = self.p / 'input.sam'
        self.sam.write_text('@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:rna\tLN:100\n' + ''.join(
            f'{name}\t{flag}\trna\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t' + 'A' * length + '\t' + 'I' * length + '\n'
            for name, flag, pos, mapq, cigar, length in self.alignments) + 'unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n')

    def bam(self, mapq=20, flags=2308, empty=False):
        stem = f'{mapq}_{flags}_{empty}'
        filtered, bam = self.p / (stem + '.unsorted.bam'), self.p / (stem + '.bam')
        cmd = ['samtools', 'view', '-b', '-F', str(flags), '-q', str(mapq)]
        if empty:
            cmd.append('-H')
        subprocess.run(cmd + ['-o', str(filtered), str(self.sam)], check=True)
        subprocess.run(['samtools', 'sort', '-o', str(bam), str(filtered)], check=True)
        return bam

    def rows(self, bams, names, strand='both'):
        out = self.p / 'counts.tsv'
        matrix(bams, names, str(self.fasta) + '.fai', out, strand)
        with out.open() as fh:
            return list(csv.reader(fh, delimiter='\t'))

    def test_coordinate_boundaries_filters_zero_fill_and_cigars(self):
        rows = self.rows([self.bam(), self.bam(empty=True)], ['B', 'empty'])
        self.assertEqual(rows[0], ['chromosome', 'strand', 'start', 'end', 'B', 'empty'])
        observed = {(r[1], int(r[2])): int(r[4]) for r in rows[1:]}
        self.assertEqual(observed, {('+', 1): 2, ('+', 6): 1, ('+', 40): 1, ('+', 80): 1,
                                    ('-', 39): 1, ('-', 57): 1, ('-', 71): 1, ('-', 100): 1})
        self.assertTrue(all(r[2] == r[3] and r[5] == '0' for r in rows[1:]))

    def test_mapq_zero_and_optional_duplicate_exclusion(self):
        rows = self.rows([self.bam(mapq=0)], ['A'])
        self.assertEqual(sum(int(r[4]) for r in rows[1:]), 11)
        rows = self.rows([self.bam(flags=3332)], ['A'], 'plus')
        self.assertEqual(sum(int(r[4]) for r in rows[1:]), 4)

    def test_empty_matrix_and_total_threshold(self):
        empty = self.rows([self.bam(empty=True)], ['empty'])
        self.assertEqual(len(empty), 1)
        source = self.p / 'source.tsv'
        source.write_text('chromosome\tstrand\tstart\tend\tB\tA\nrna\t+\t1\t1\t6\t4\nrna\t+\t2\t2\t9\t0\n')
        out = self.p / 'threshold.tsv'
        threshold(source, out, 10)
        self.assertEqual(out.read_text().splitlines()[-1], 'rna\t+\t1\t1\t6\t4')

    def test_unchanged_legacy_counter_when_modern_bash_available(self):
        bash = os.environ.get('CAPZYME_TEST_BASH', 'bash')
        version = subprocess.check_output([bash, '-c', 'echo ${BASH_VERSINFO[0]}'], text=True).strip()
        if int(version) < 4:
            self.skipTest('Unchanged legacy counter requires Bash >=4; set CAPZYME_TEST_BASH')
        bam = self.bam(mapq=0)
        out = self.p / 'legacy.tsv'
        subprocess.run([bash, str(ROOT / 'legacy_scripts/count5ends_matrix_v1_3.sh'), '-g', str(self.fasta),
                        '--tmp-dir', str(self.p), '-o', str(out), str(bam)], check=True, capture_output=True)
        current = self.rows([self.bam()], ['sample'])
        legacy = list(csv.reader(io.StringIO(out.read_text()), delimiter='\t'))
        self.assertEqual(current[0][:4], legacy[0][:4])
        self.assertEqual(sorted(current[1:]), sorted(legacy[1:]))

    def test_legacy_multisample_threshold_and_strands(self):
        bash = os.environ.get('CAPZYME_TEST_BASH', 'bash')
        if int(subprocess.check_output([bash, '-c', 'echo ${BASH_VERSINFO[0]}'], text=True).strip()) < 4:
            self.skipTest('Unchanged legacy counter requires Bash >=4')
        adjacent = self.p / 'adjacent.sam'
        adjacent.write_text('@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:rna\tLN:100\n'
                            'a\t0\trna\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n'
                            'b\t0\trna\t2\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n')
        second = self.p / 'adjacent.bam'
        subprocess.run(['samtools', 'sort', '-o', str(second), str(adjacent)], check=True)
        sources = [second, self.bam(mapq=0), self.bam(empty=True)]
        names = ['B', 'A', 'empty']
        for mapq, strand, minimum in [(20, 'both', 0), (20, 'both', 2), (0, 'plus', 0), (30, 'minus', 0)]:
            with self.subTest(mapq=mapq, strand=strand, minimum=minimum):
                current = self.p / 'new.tsv'
                matrix([second, self.bam(mapq=mapq), sources[2]], names, str(self.fasta) + '.fai', current, strand)
                filtered = self.p / 'new.filtered.tsv'
                threshold(current, filtered, minimum)
                old = self.p / 'old.tsv'
                subprocess.run([bash, str(ROOT / 'legacy_scripts/count5ends_matrix_v1_3.sh'),
                                '-g', str(self.fasta), '--tmp-dir', str(self.p), '--mapq', str(mapq),
                                '--strand', strand, '--min-count', str(minimum), '-o', str(old),
                                *map(str, sources)], check=True, capture_output=True)
                expected = list(csv.reader(io.StringIO(old.read_text()), delimiter='\t'))
                observed = list(csv.reader(io.StringIO(filtered.read_text()), delimiter='\t'))
                self.assertEqual(sorted(observed[1:]), sorted(expected[1:]))


@unittest.skipUnless(all(shutil.which(x) for x in ('snakemake', 'cutadapt', 'bowtie2', 'samtools', 'bedtools')),
                     'Full local workflow dependencies required')
class WorkflowTests(unittest.TestCase):
    def test_synthetic_dry_run_execution_restart_and_threshold(self):
        with tempfile.TemporaryDirectory(prefix='capzyme test ') as tmp:
            cfg = generate(tmp)
            cmd = ['snakemake', '--snakefile', str(ROOT / 'workflow/Snakefile'), '--configfile', str(cfg),
                   '--directory', tmp, '--cores', '4', '--printshellcmds']
            dry = subprocess.run(cmd + ['--dry-run'], text=True, capture_output=True)
            self.assertEqual(dry.returncode, 0, dry.stdout + dry.stderr)
            run = subprocess.run(cmd, text=True, capture_output=True)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            p = Path(tmp)
            self.assertEqual((p / 'run/counts/5ends.raw.tsv').read_text(), (p / 'expected.tsv').read_text())
            filtered = (p / 'run/counts/5ends.min_total_10.tsv').read_text().splitlines()
            self.assertEqual(len(filtered), 2)
            self.assertEqual(filtered[1], 'RNU2__n2\t+\t1\t1\t12\t0\t0')
            provenance = json.loads((p / 'run/provenance/run.json').read_text())
            self.assertEqual(provenance['settings']['mapq'], 20)
            self.assertIn('snakemake', provenance['versions'])
            again = subprocess.run(cmd + ['--dry-run'], text=True, capture_output=True)
            self.assertEqual(again.returncode, 0, again.stdout + again.stderr)
            self.assertIn('Nothing to be done', again.stdout + again.stderr)


if __name__ == '__main__':
    unittest.main()
