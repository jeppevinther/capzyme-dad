"""Bounded, server-side old/new pilot. No host or experimental paths are embedded."""
import argparse
from collections import Counter
import csv
import datetime
import gzip
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import threading
import time

from core import records, sha256
from fixed_reference import check_fixed_hashes, verify_fixed

LIMIT = 100000
GIB = 1024 ** 3


def subset_fastq(source, output, limit=LIMIT):
    if not 1 <= limit <= LIMIT:
        raise ValueError('Pilot read limit must be between 1 and 100000')
    opener = gzip.open if str(source).endswith('.gz') else open
    ids = set()
    digest = hashlib.sha256()
    with opener(source, 'rt', newline='') as inp, open(output, 'w', newline='') as out:
        for _ in range(limit):
            block = [inp.readline() for _ in range(4)]
            if any(not line for line in block):
                raise ValueError('Source has fewer complete FASTQ records than the pilot limit')
            h, seq, plus, quality = [line.rstrip('\r\n') for line in block]
            if not h.startswith('@') or not plus.startswith('+') or len(seq) != len(quality) or not seq:
                raise ValueError('Malformed FASTQ record within bounded pilot prefix')
            name = h[1:].split()[0]
            if name in ids:
                raise ValueError('Duplicate read ID within pilot prefix')
            ids.add(name)
            # Preserve original line bytes, names, sequence, qualities and ordering.
            text = ''.join(block)
            if not block[-1].endswith('\n'):
                text += '\n'
            out.write(text)
            digest.update(text.encode())
    return {'records': limit, 'sha256': digest.hexdigest(), 'unique_ids': len(ids)}


def sam_records(path):
    if str(path).endswith('.bam'):
        proc = subprocess.Popen(['samtools', 'view', str(path)], stdout=subprocess.PIPE, text=True)
        stream = proc.stdout
    else:
        proc = None
        stream = open(path)
    rows = Counter()
    try:
        for line in stream:
            if line.startswith('@'):
                continue
            fields = line.rstrip('\n').split('\t')
            key = tuple(fields[:11]) + tuple(sorted(fields[11:]))
            rows[key] += 1
    finally:
        stream.close()
    if proc is not None and proc.wait() != 0:
        raise RuntimeError('samtools view failed during comparison')
    return rows


def matrix_rows(path, expected_columns):
    with open(path) as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader)
        if header[:4] != ['chromosome', 'strand', 'start', 'end'] or header[4:] != expected_columns:
            raise ValueError('Unexpected matrix column order/schema: ' + str(path))
        rows = {}
        for row in reader:
            key = (row[0], row[1], int(row[2]), int(row[3]))
            if key in rows:
                raise ValueError('Duplicate matrix feature')
            rows[key] = tuple(map(int, row[4:]))
    return rows


class Pilot:
    def __init__(self, config, root):
        import psutil
        self.psutil = psutil
        self.config, self.root = config, root
        self.stop = threading.Event()
        self.done = threading.Event()
        self.reason = None
        self.current = None
        self.resources = {'peak_tree_rss_bytes': 0, 'max_load1': 0, 'min_available_memory_bytes': None,
                          'nice': os.getpriority(os.PRIO_PROCESS, 0), 'samples': 0}
        self.process = psutil.Process()
        self.affinity = sorted(os.sched_getaffinity(0))[:4]
        os.sched_setaffinity(0, self.affinity)
        self.resources['allowed_cpus'] = self.affinity

    def host_check(self):
        vm = self.psutil.virtual_memory()
        load = os.getloadavg()[0]
        if vm.available < max(64 * GIB, vm.total * .10):
            raise RuntimeError('Host memory pressure: stop rather than changing host')
        if load > self.psutil.cpu_count() * .75:
            raise RuntimeError('Host CPU load too high: stop rather than changing host')
        return {'available_bytes': vm.available, 'total_bytes': vm.total, 'load1': load}

    def monitor(self):
        busy_streak = 0
        with open(self.root / 'resources.jsonl', 'w') as fh:
            while not self.done.is_set():
                vm = self.psutil.virtual_memory()
                cpu_busy = self.psutil.cpu_percent(interval=.5)
                children = self.process.children(recursive=True)
                rss = self.process.memory_info().rss
                for child in children:
                    try:
                        rss += child.memory_info().rss
                    except self.psutil.Error:
                        pass
                load = os.getloadavg()[0]
                record = {'utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                          'tree_rss_bytes': rss, 'host_available_bytes': vm.available,
                          'host_cpu_busy_percent': cpu_busy, 'host_load1': load}
                fh.write(json.dumps(record) + '\n'); fh.flush()
                r = self.resources
                r['samples'] += 1
                r['peak_tree_rss_bytes'] = max(r['peak_tree_rss_bytes'], rss)
                r['max_load1'] = max(r['max_load1'], load)
                r['min_available_memory_bytes'] = min(r['min_available_memory_bytes'] or vm.available, vm.available)
                overloaded = cpu_busy > 85 or load > self.psutil.cpu_count() * .75 or vm.available < max(64 * GIB, vm.total * .10)
                busy_streak = busy_streak + 1 if overloaded else 0
                if rss > 16 * GIB or busy_streak >= 10:
                    self.reason = 'Pilot RSS exceeds 16 GiB' if rss > 16 * GIB else 'Host sustained heavy load/memory pressure'
                    self.stop.set()
                    if self.current is not None and self.current.poll() is None:
                        os.killpg(self.current.pid, signal.SIGTERM)
                    return
                self.done.wait(.5)

    def run(self, name, command, shell=False):
        if self.stop.is_set():
            raise RuntimeError(self.reason)
        self.host_check()
        started = time.time()
        record = {'name': name, 'command': command, 'shell': shell,
                  'utc_start': datetime.datetime.now(datetime.timezone.utc).isoformat()}
        with open(self.root / 'commands.jsonl', 'a') as fh:
            fh.write(json.dumps(record) + '\n')
        cmd = ['/usr/bin/time', '-v', '-o', str(self.root / 'logs' / (name + '.time.txt'))]
        cmd += ['bash', '-o', 'pipefail', '-c', command] if shell else command
        with open(self.root / 'logs' / (name + '.log'), 'w') as log:
            self.current = subprocess.Popen(cmd, cwd=self.root, stdout=log, stderr=subprocess.STDOUT,
                                            start_new_session=True)
            code = self.current.wait()
        record.update(exit_status=code, elapsed_seconds=time.time() - started)
        with open(self.root / 'commands.jsonl', 'a') as fh:
            fh.write(json.dumps(record) + '\n')
        self.current = None
        if code != 0 or self.stop.is_set():
            raise RuntimeError(self.reason or f'Command failed: {name}, exit {code}')

    def main(self):
        import shlex
        q = lambda p: shlex.quote(str(p))
        cfg, root = self.config, self.root
        source = root / 'source'
        self.host_check()
        if self.resources['nice'] < 10:
            raise RuntimeError('Pilot must be launched with nice -n 10')
        required = {'cutadapt': '4.4', 'bowtie2': '2.4.2', 'samtools': '1.21', 'bedtools': 'v2.31.0', 'snakemake': '9.9.0'}
        versions = {}
        for executable, expected in required.items():
            output = subprocess.check_output([executable, '--version'], text=True)
            versions[executable] = output.splitlines()[0]
            if expected not in output.splitlines()[0]:
                raise RuntimeError('Unexpected executable version: ' + executable)
        if len(cfg['samples']) != 4 or [x['sample_id'] for x in cfg['samples']] != ['01', '02', '03', '04']:
            raise ValueError('Only the four authorized libraries, in confirmed order, may run')
        if cfg.get('read_limit') != LIMIT:
            raise ValueError('This pilot is fixed at 100000 records per sample')
        fixed = cfg['reference']
        before_hashes = check_fixed_hashes(fixed)
        inputs = [Path(fixed['fasta']), *(Path(fixed['index_prefix'] + s) for s in fixed['index_sha256'])]
        fai = Path(fixed['fasta'] + '.fai')
        if not fai.is_file():
            raise ValueError('Historical .fai must already exist; do not create one in source storage')
        expected_sizes = [(n, len(s)) for n, s in records(fixed['fasta'])]
        actual_sizes = [(r[0], int(r[1])) for r in csv.reader(fai.open(), delimiter='\t')]
        if expected_sizes != actual_sizes:
            raise ValueError('Historical .fai does not match fixed FASTA')
        inputs.append(fai)
        before_fai_hash = sha256(fai)
        for item in cfg['samples']:
            path = Path(item['raw_fastq'])
            st = path.stat()
            if st.st_size != item['expected_size_bytes'] or not os.access(path, os.R_OK):
                raise ValueError('FASTQ metadata/readability check failed')
            inputs.append(path)
        old_helper = source / 'legacy_scripts/count5ends_matrix_v1_3.sh'
        if sha256(old_helper) != '7bd16ba82d855b1b4bfad32e9dc8795a6a8a4081d9d600d5166337ebb79303b0':
            raise ValueError('Legacy counting script checksum mismatch')
        pre_stats = {str(p): {'size': p.stat().st_size, 'mtime_ns': p.stat().st_mtime_ns} for p in inputs}
        for folder in ['logs', 'subsets', 'baseline', 'preflight']:
            (root / folder).mkdir(exist_ok=True)
        verify_fixed(fixed, root / 'preflight/sizes.tsv', root / 'preflight/reference.json')
        audit = {'versions': versions, 'preflight_host': self.host_check(), 'input_stats_before': pre_stats,
                 'subset_policy': 'first 100000 complete records, shared by both arms',
                 'memory_policy': 'mem_mb=16384 scheduling budget, monitored RSS; nice does not cap memory'}
        (root / 'preflight.json').write_text(json.dumps(audit, indent=2) + '\n')
        monitor = threading.Thread(target=self.monitor, daemon=True); monitor.start()
        results = {'status': 'running', 'samples': {}, 'versions': versions}
        try:
            for item in cfg['samples']:
                sample = item['sample_id']
                dest = root / 'subsets' / (sample + '.fastq')
                subset_start = time.time()
                subset = subset_fastq(item['raw_fastq'], dest)
                results['samples'][sample] = {'subset': subset, 'subset_seconds': time.time() - subset_start}
                if self.stop.is_set():
                    raise RuntimeError(self.reason)
            manifest = root / 'samples.tsv'
            manifest.write_text('sample_id\tfastq\n' + ''.join(s['sample_id'] + '\t' + str(root / 'subsets' / (s['sample_id'] + '.fastq')) + '\n' for s in cfg['samples']))
            new_cfg = {'schema_version': 1, 'samples': str(manifest), 'output_dir': str(root / 'new'), 'reference': fixed,
                       'preprocessing': {'adapter': 'AGATCGGAAGAGCACACGTCT', 'quality': 20, 'minimum_length': 15, 'nextseq_trim': 20},
                       'mapping': {'norc': True, 'preset': 'very-sensitive', 'seed': 0, 'threads': 4, 'keep_diagnostic_sam': True},
                       'counting': {'mapq': 20, 'strand': 'plus', 'exclude_duplicates': False, 'min_total_count': 10}}
            config = root / 'pilot.config.json'; config.write_text(json.dumps(new_cfg, indent=2) + '\n')
            base_cmd = ['snakemake', '--snakefile', str(source / 'workflow/Snakefile'), '--configfile', str(config),
                        '--profile', str(source / 'profiles/pilot'), '--directory', str(root)]
            self.run('new_dry_run', base_cmd + ['--dry-run'])
            # Both arms consume the exact same subsets; all commands inherit affinity and nice.
            for item in cfg['samples']:
                sample = item['sample_id']
                subset = root / 'subsets' / (sample + '.fastq')
                trimmed = root / 'baseline' / (sample + '_trimmed.fastq.gz')
                bam = root / 'baseline' / (sample + '_mapped.bam')
                accepted = root / 'baseline' / (sample + '_accepted.bam')
                self.run('old_trim_' + sample, ['cutadapt', '-a', 'AGATCGGAAGAGCACACGTCT', '-q', '20', '-m', '15', '--nextseq-trim=20', str(subset), '-o', str(trimmed)])
                self.run('old_map_' + sample,
                         f'bowtie2 -p 1 --norc --very-sensitive -x {q(fixed["index_prefix"])} -U {q(trimmed)} '
                         f'2> {q(root / "logs" / ("old_bowtie2_" + sample + ".log"))} | '
                         f'samtools view -b -@ 0 - | samtools sort -@ 0 -m 256M -T {q(root / "baseline" / (sample + ".sorttmp"))} -o {q(bam)} -', shell=True)
                self.run('old_filter_' + sample, ['samtools', 'view', '-b', '-F', '2308', '-q', '20', '-o', str(accepted), str(bam)])
            baseline_bams = [str(root / 'baseline' / (x['sample_id'] + '_mapped.bam')) for x in cfg['samples']]
            for minimum in [0, 10]:
                self.run('old_counts_' + str(minimum), ['bash', str(old_helper), '-g', fixed['fasta'], '--mapq', '20', '--strand', 'plus',
                         '--threads', '1', '--tmp-dir', str(root / 'baseline'), '--min-count', str(minimum),
                         '-o', str(root / 'baseline' / ('counts_' + str(minimum) + '.tsv')), *baseline_bams])
            self.run('new_workflow', base_cmd)
            self.run('new_restart_dry_run', base_cmd + ['--dry-run'])
            for item in cfg['samples']:
                sample = item['sample_id']
                old_trim = root / 'baseline' / (sample + '_trimmed.fastq.gz')
                new_trim = root / 'new/trimmed' / (sample + '.fastq.gz')
                def decompressed_digest(path):
                    h = hashlib.sha256(); lines = 0
                    with gzip.open(path, 'rb') as fh:
                        for line in fh:
                            h.update(line); lines += 1
                    if lines % 4:
                        raise ValueError('Incomplete trimmed FASTQ')
                    return {'sha256': h.hexdigest(), 'records': lines // 4}
                old_d, new_d = decompressed_digest(old_trim), decompressed_digest(new_trim)
                old_aln = sam_records(root / 'baseline' / (sample + '_mapped.bam'))
                new_aln = sam_records(root / 'new/diagnostics' / (sample + '.unfiltered.sam'))
                old_accepted = sam_records(root / 'baseline' / (sample + '_accepted.bam'))
                new_accepted = sam_records(root / 'new/bam' / (sample + '.accepted.bam'))
                record = results['samples'][sample]
                record.update(old_trim=old_d, new_trim=new_d, trimming_equal=old_d == new_d,
                              alignments_equal=old_aln == new_aln, old_alignment_records=sum(old_aln.values()),
                              new_alignment_records=sum(new_aln.values()), old_accepted_records=sum(old_accepted.values()),
                              new_accepted_records=sum(new_accepted.values()), filtering_equal=old_accepted == new_accepted,
                              alignment_records_only_old=sum((old_aln-new_aln).values()),
                              alignment_records_only_new=sum((new_aln-old_aln).values()))
                if not all(record[k] for k in ['trimming_equal', 'alignments_equal', 'filtering_equal']):
                    raise ValueError('Comparison discrepancy in library ' + sample)
            old_cols = [s['sample_id'] + '_mapped' for s in cfg['samples']]
            new_cols = [s['sample_id'] for s in cfg['samples']]
            counts = {}
            for minimum in [0, 10]:
                old_rows = matrix_rows(root / 'baseline' / ('counts_' + str(minimum) + '.tsv'), old_cols)
                new_path = root / 'new/counts' / ('5ends.raw.tsv' if minimum == 0 else '5ends.min_total_10.tsv')
                new_rows = matrix_rows(new_path, new_cols)
                counts[str(minimum)] = {'equal': old_rows == new_rows, 'old_rows': len(old_rows), 'new_rows': len(new_rows),
                                       'sample_sums': [sum(r[i] for r in new_rows.values()) for i in range(4)]}
                if old_rows != new_rows:
                    raise ValueError('Count matrix comparison discrepancy')
            results['counts'] = counts
            results['reference_hashes_unchanged'] = check_fixed_hashes(fixed) == before_hashes
            results['fai_hash_unchanged'] = sha256(fai) == before_fai_hash
            post_stats = {str(p): {'size': p.stat().st_size, 'mtime_ns': p.stat().st_mtime_ns} for p in inputs}
            results['input_metadata_unchanged'] = post_stats == pre_stats
            if not all(results[k] for k in ['reference_hashes_unchanged', 'fai_hash_unchanged', 'input_metadata_unchanged']):
                raise ValueError('Historical input integrity changed')
            results['status'] = 'passed'
        except Exception as exc:
            results['status'] = 'failed'; results['error'] = str(exc)
            raise
        finally:
            self.done.set(); monitor.join(timeout=3)
            results['resources'] = self.resources
            results['stop_reason'] = self.reason
            (root / 'comparison_results.json').write_text(json.dumps(results, indent=2) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True)
    parser.add_argument('--run-directory', required=True)
    args = parser.parse_args()
    root = Path(args.run_directory).resolve()
    if (root / 'comparison_results.json').exists() or (root / 'subsets').exists():
        raise SystemExit('Refusing to overwrite an existing pilot run')
    Pilot(json.loads(Path(args.config).read_text()), root).main()
