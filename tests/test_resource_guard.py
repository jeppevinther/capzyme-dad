"""Tiny Linux-only operational tests; no experimental data or shared installs."""
import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile
import time
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'workflow/scripts'))
from resource_guard import supervise


def stalled_scan(root):
    time.sleep(3)
    return {'apparent_bytes': 0, 'free_bytes': 10**15}


@unittest.skipUnless(sys.platform == 'linux' and importlib.util.find_spec('psutil'), 'Linux + psutil required')
class GuardTests(unittest.TestCase):
    def test_stalled_storage_does_not_block_sampling_or_stop(self):
        with tempfile.TemporaryDirectory() as d:
            result = supervise([sys.executable, '-c', 'import time; time.sleep(10)'], Path(d)/'audit', d,
                               scan=stalled_scan, interval=.05, storage_timeout=.4, rss_mb=16384)
            self.assertEqual(result['stop_reason'], 'storage_observation_stale')
            self.assertNotEqual(result['exit_status'], 0)
            self.assertGreaterEqual(result['samples'], 5)
            self.assertLess(result['max_sample_gap_seconds'], .3)
            self.assertLess(result['elapsed_seconds'], 5)
            self.assertEqual(json.loads((Path(d)/'audit/status.json').read_text())['status'], 'failed')

    def test_normal_completion_and_logs(self):
        with tempfile.TemporaryDirectory() as d:
            result = supervise([sys.executable, '-c', 'import time; print("tiny fixture"); time.sleep(.3)'], Path(d)/'audit', d,
                               interval=.05, rss_mb=16384, free_gib=1)
            self.assertEqual(result['status'], 'passed')
            self.assertIn('tiny fixture', (Path(d)/'audit/command.log').read_text())


if __name__ == '__main__': unittest.main()
