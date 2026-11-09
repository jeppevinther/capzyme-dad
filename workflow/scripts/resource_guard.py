#!/usr/bin/env python3
"""Linux process-tree supervision with independent, bounded storage observation.

Requires psutil. Fast sampling writes to node-local temporary storage. A separate
process performs slow storage scans; stale scans stop the workload, never block
RSS/load/stop handling. This is monitoring, not an OS memory/cgroup limit.
"""
import argparse
import json
import multiprocessing as mp
import os
from pathlib import Path
import queue
import shutil
import signal
import subprocess
import tempfile
import time


def storage_scan(root):
    size = 0
    for directory, _, files in os.walk(root):
        for name in files:
            p = Path(directory)/name
            try:
                if not p.is_symlink(): size += p.stat().st_size
            except FileNotFoundError:
                pass  # Transient tool scratch files may disappear during a scan.
    return {'apparent_bytes': size, 'free_bytes': shutil.disk_usage(root).free}


def storage_worker(root, results, interval, scan):
    while True:
        try:
            value = scan(root)
            results.put((time.monotonic(), value, None))
        except Exception as exc: results.put((time.monotonic(), None, str(exc)))
        time.sleep(interval)


def descendants(process, known):
    import psutil
    try: known.update({p.pid: p for p in process.children(recursive=True)})
    except psutil.Error: pass
    live = []
    for pid, p in list(known.items()):
        try:
            if p.is_running() and p.status() != psutil.STATUS_ZOMBIE: live.append(p)
        except psutil.Error: pass
    return live


def terminate(processes):
    import psutil
    # Descendants may start separate sessions; terminate by PID/create-time identity.
    for p in reversed(processes):
        try: p.terminate()
        except psutil.Error: pass
    # Do not reap Popen's child via psutil.wait_procs: that would lose its real
    # exit status before Popen.wait(). Poll /proc instead; Popen owns waitpid.
    deadline = time.monotonic()+3
    alive = list(processes)
    while alive and time.monotonic() < deadline:
        still = []
        for p in alive:
            try:
                if p.is_running() and p.status() != psutil.STATUS_ZOMBIE: still.append(p)
            except psutil.Error: pass
        alive = still
        if alive: time.sleep(.05)
    for p in alive:
        try: p.kill()
        except psutil.Error: pass


def supervise(command, audit, storage_root, *, cores=4, rss_mb=16384,
              storage_gib=20, free_gib=50, interval=1.0, storage_interval=15,
              storage_timeout=60, pressure_seconds=10, scan=storage_scan):
    import psutil
    audit = Path(audit).resolve()
    # All real-time I/O, including STOP, remains local even on a stalled project mount.
    local = Path(tempfile.mkdtemp(prefix='capzyme-monitor-', dir='/tmp'))
    print(json.dumps({'local_monitor_directory': str(local), 'stop_file': str(local/'STOP')}), flush=True)
    started = time.monotonic()
    ctx = mp.get_context('fork')
    results = ctx.Queue()
    worker = ctx.Process(target=storage_worker, args=(str(storage_root), results, storage_interval, scan), daemon=True)
    worker.start()
    flags = {'requested': False}
    def stopped(signum, frame): flags['requested'] = True
    old = {s: signal.signal(s, stopped) for s in (signal.SIGTERM, signal.SIGINT)}
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
    log = (local/'command.log').open('w')
    child = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, env=env, start_new_session=True)
    process = psutil.Process(child.pid)
    if hasattr(process, 'cpu_affinity'):
        process.cpu_affinity(process.cpu_affinity()[:cores])
    known = {process.pid: process}
    psutil.cpu_percent()
    last_storage = started
    storage = None
    peak = 0
    samples = 0
    max_gap = 0
    previous = started
    pressure_start = None
    stop_reason = None
    summary = {'command': command, 'local_monitor_directory': str(local), 'nice': os.getpriority(os.PRIO_PROCESS, 0),
               'limits': {'cores': cores, 'rss_mb': rss_mb, 'storage_gib': storage_gib, 'free_gib': free_gib,
                          'storage_timeout_seconds': storage_timeout, 'pressure_seconds': pressure_seconds}}
    try:
        with (local/'resources.jsonl').open('w', buffering=1) as output:
            while True:
                tick = time.monotonic()
                live = descendants(process, known)
                rss = psutil.Process().memory_info().rss
                if worker.is_alive():
                    try: rss += psutil.Process(worker.pid).memory_info().rss
                    except psutil.Error: pass
                for p in live:
                    try: rss += p.memory_info().rss
                    except psutil.Error: pass
                memory = psutil.virtual_memory()
                cpu = psutil.cpu_percent()
                load = os.getloadavg()[0]
                # Acquisition timestamp precedes storage-result handling and records actual fast sampling.
                acquired = time.monotonic()
                peak = max(peak, rss); samples += 1
                max_gap = max(max_gap, acquired-previous); previous = acquired
                try:
                    while True:
                        finished, storage, error = results.get_nowait()
                        last_storage = finished
                        if error: stop_reason = 'storage_scan_failed: ' + error
                except queue.Empty: pass
                high = cpu > 85 or load > .75 * (psutil.cpu_count() or 1)
                pressure_start = (pressure_start or acquired) if high else None
                reason = stop_reason
                if flags['requested'] or (local/'STOP').exists(): reason = 'stop_requested'
                elif rss > rss_mb * 1024**2: reason = 'rss_limit'
                elif memory.available < max(64*1024**3, memory.total*.1): reason = 'host_memory_reserve'
                elif pressure_start and acquired-pressure_start >= pressure_seconds: reason = 'sustained_host_pressure'
                elif acquired-last_storage > storage_timeout: reason = 'storage_observation_stale'
                elif storage and storage['apparent_bytes'] > storage_gib*1024**3: reason = 'storage_allowance'
                elif storage and storage['free_bytes'] < free_gib*1024**3: reason = 'filesystem_free_floor'
                row = {'elapsed_s': acquired-started, 'acquired_unix_s': time.time(), 'rss_bytes': rss,
                       'available_ram': memory.available, 'host_cpu_percent': cpu, 'load1': load,
                       'storage_age_s': acquired-last_storage, 'storage': storage, 'stop_reason': reason}
                output.write(json.dumps(row)+'\n')
                if reason:
                    stop_reason = reason; terminate(descendants(process, known)); break
                if child.poll() is not None:
                    remaining = [p for p in descendants(process, known) if p.pid != child.pid]
                    if remaining:
                        stop_reason = 'orphan_descendants'; terminate(remaining)
                    break
                time.sleep(max(0, interval-(time.monotonic()-tick)))
        rc = child.wait(timeout=10)
        summary.update(exit_status=rc, stop_reason=stop_reason, elapsed_seconds=time.monotonic()-started,
                       sampled_peak_rss_bytes=peak, samples=samples, max_sample_gap_seconds=max_gap,
                       last_storage=storage, status='passed' if rc == 0 and not stop_reason else 'failed')
    finally:
        terminate(descendants(process, known))
        log.close()
        worker.terminate(); worker.join(timeout=1)
        for s, handler in old.items(): signal.signal(s, handler)
        (local/'status.json').write_text(json.dumps(summary, indent=2)+'\n')
    # Workload is over. Export in a separate bounded process; preserve node-local logs on failure.
    export = subprocess.run([shutil.which('python3') or 'python3', '-c',
                             'import shutil,sys; shutil.copytree(sys.argv[1],sys.argv[2])', str(local), str(audit)], timeout=60)
    if export.returncode: raise RuntimeError('Audit export failed; node-local logs retained: ' + str(local))
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--audit-dir', required=True, help='Fresh external directory; compact logs copied here after execution')
    p.add_argument('--storage-root', required=True)
    p.add_argument('--cores', type=int, default=4)
    p.add_argument('--rss-mb', type=int, default=16384)
    p.add_argument('--storage-gib', type=int, default=20)
    p.add_argument('--free-gib', type=int, default=50)
    p.add_argument('command', nargs=argparse.REMAINDER)
    a = p.parse_args()
    command = a.command[1:] if a.command[:1] == ['--'] else a.command
    if not command or min(a.cores, a.rss_mb, a.storage_gib, a.free_gib) <= 0: p.error('Positive limits and a command required')
    if Path(a.audit_dir).exists(): p.error('audit-dir must be new')
    result = supervise(command, a.audit_dir, a.storage_root, cores=a.cores, rss_mb=a.rss_mb,
                       storage_gib=a.storage_gib, free_gib=a.free_gib)
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result['status'] == 'passed' else 1)


if __name__ == '__main__': main()
