import gzip
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/scripts'))
from pilot_validation import subset_fastq


class BoundedSubsetTests(unittest.TestCase):
    def test_exact_prefix_and_no_tail_processing(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            text = '@a\nACG\n+\nIII\n@b\nTTA\n+\nHHH\n'
            source = p / 'reads.fastq.gz'
            with gzip.open(source, 'wt') as fh:
                fh.write(text + 'malformed tail must not be read\n')
            result = subset_fastq(source, p / 'out.fastq', 2)
            self.assertEqual((p / 'out.fastq').read_text(), text)
            self.assertEqual(result['records'], 2)
            self.assertEqual(result['sha256'], hashlib.sha256(text.encode()).hexdigest())

    def test_truncated_malformed_duplicate_and_over_limit_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            for text, limit in [('@a\nAC\n+\nI\n', 1), ('@a\nAC\n+\n', 1),
                                ('@a\nAC\n+\nII\n@a\nTT\n+\nII\n', 2),
                                ('@a\nAC\n+\nII\n', 2)]:
                (p / 'reads.fastq').write_text(text)
                with self.subTest(text=text), self.assertRaises(ValueError):
                    subset_fastq(p / 'reads.fastq', p / 'out.fastq', limit)
            with self.assertRaises(ValueError):
                subset_fastq(p / 'reads.fastq', p / 'out.fastq', 100001)
