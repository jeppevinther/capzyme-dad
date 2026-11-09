"""Independent all-value verification, including deliberate corruptions."""
import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'workflow/scripts'))
from core import sha256
from validate_bigwigs import validate_space
from test_bigwig import HAS_R, RSCRIPT


@unittest.skipUnless(importlib.util.find_spec('pyBigWig'), 'pyBigWig required')
class VerifierTests(unittest.TestCase):
    @unittest.skipUnless(HAS_R, 'actual DESeq2 required')
    def test_actual_deseq_independent_factor_and_all_interval_verification(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d); matrix=p/'matrix.tsv'; sizes=p/'sizes.tsv'
            matrix.write_text('chromosome\tstrand\tstart\tend\tB\tA\nrna\t+\t1\t1\t7\t14\nrna\t-\t10\t10\t3\t6\n')
            sizes.write_text('rna\t10\n')
            out=p/'deseq'
            subprocess.run([sys.executable,str(ROOT/'workflow/scripts/bigwig2.py'),'--matrix',str(matrix),'--assert-unthresholded','--chrom-sizes',str(sizes),'--cohort','B','A','--normalization','deseq2','--deseq2-type','ratio','--rscript',RSCRIPT,'--output-dir',str(out)],check=True,capture_output=True)
            factors=p/'direct_R/factors.tsv'
            subprocess.run([RSCRIPT,'--vanilla',str(ROOT/'tests/independent_deseq2_factors.R'),str(matrix),str(p/'direct_R'),'B','A'],check=True,capture_output=True)
            spec={'matrix':str(matrix),'matrix_sha256':sha256(matrix),'sizes':str(sizes),'samples':['B','A'],'accepted_depths':[10,20],'depth_scope':'artificial accepted records',
                  'arms':[{'method':'deseq2','estimator':'ratio','independent_factors':str(factors),'independent_factors_sha256':sha256(factors),'directory':str(out),'strands':['+','-']}]}
            result=validate_space(spec,p/'verify')
            self.assertEqual(result['arms']['deseq2']['track_count'],4)
            self.assertLessEqual(result['arms']['deseq2']['max_factor_relative_difference'],1e-12)
            with factors.open() as fh: rows=list(csv.DictReader(fh,delimiter='\t'))
            rows[0]['size_factor']=str(float(rows[0]['size_factor'])*1.1)
            with factors.open('w') as f:
                writer=csv.DictWriter(f,fieldnames=['sample_id','size_factor'],delimiter='\t');writer.writeheader();writer.writerows(rows)
            spec['arms'][0]['independent_factors_sha256']=sha256(factors)
            with self.assertRaisesRegex(ValueError,'size factor mismatch'): validate_space(spec,p/'bad-factor')

    def test_all_values_coordinates_dictionaries_empty_and_rejection(self):
        import pyBigWig
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            matrix = p/'matrix.tsv'; matrix.write_text('chromosome\tstrand\tstart\tend\tB\tA\nrna\t+\t1\t1\t7\t1\nrna\t+\t10\t10\t3\t2\n')
            sizes = p/'sizes.tsv'; sizes.write_text('rna\t10\nother\t20\n')
            depths = p/'depths.tsv'; depths.write_text('sample_id\taccepted_records\tscope_id\nB\t10\tRNA_accepted\nA\t3\tRNA_accepted\n')
            arms = []
            for method in ('raw', 'rpm'):
                out = p/method
                cmd = [sys.executable, str(ROOT/'workflow/scripts/bigwig2.py'), '--matrix', str(matrix), '--assert-unthresholded', '--chrom-sizes', str(sizes), '--cohort', 'B', 'A', '--normalization', method, '--output-dir', str(out)]
                if method == 'rpm': cmd += ['--depths', str(depths)]
                subprocess.run(cmd, check=True, capture_output=True)
                arms.append({'method': method, 'directory': str(out), 'strands': ['+', '-']})
            spec = {'matrix': str(matrix), 'matrix_sha256': sha256(matrix), 'sizes': str(sizes), 'samples': ['B', 'A'], 'accepted_depths': [10, 3], 'depth_scope': 'RNA_accepted', 'arms': arms}
            work = p/'good'; work.mkdir()
            report = validate_space(spec, work)
            self.assertEqual(report['status'], 'passed')
            self.assertEqual(sum(t['nonzero_values_checked'] for arm in report['arms'].values() for t in arm['tracks']), 8)
            self.assertEqual(sum(t['empty_zero_sentinel'] for arm in report['arms'].values() for t in arm['tracks']), 4)
            # Even if somebody re-publishes a checksum, wrong decoded values must fail.
            track = p/'raw/tracks/B.plus.bw'
            with pyBigWig.open(str(track), 'w') as bw:
                bw.addHeader([('rna', 10), ('other', 20)])
                bw.addEntries(['rna', 'rna'], [0, 9], ends=[1, 10], values=[8., 3.])
            provenance = p/'raw/provenance.json'; meta = json.loads(provenance.read_text())
            meta['outputs_sha256']['tracks/B.plus.bw'] = sha256(track); provenance.write_text(json.dumps(meta))
            work = p/'bad'; work.mkdir()
            with self.assertRaisesRegex(ValueError, 'interval/value mismatch'): validate_space(spec, work)
            spec['matrix_sha256'] = '0'*64
            work = p/'bad-source'; work.mkdir()
            with self.assertRaisesRegex(ValueError, 'trusted checksum'): validate_space(spec, work)


if __name__ == '__main__': unittest.main()
