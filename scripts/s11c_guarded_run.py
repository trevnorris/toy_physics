#!/usr/bin/env python3
"""Run one S11c job in a bounded user service, with durable resource telemetry.

Requires a host shell and a systemd user manager. Fails closed if the actual
cgroup limits cannot be verified. No scientific imports or automatic retries.
"""
import argparse
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import uuid

ROOT = Path(__file__).resolve().parents[1]
STORE = ROOT / '_scratch/s11c'
PAUSE = STORE / 'PAUSED_HOST_FREEZE.json'
MEMORY_MAX = 2 * 1024**3
AVAILABLE_MIN = 4 * 1024**3
THREADS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS')


def save(path, value):
    temporary = path.with_suffix(path.suffix + '.new')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def available_memory():
    for line in Path('/proc/meminfo').read_text().splitlines():
        if line.startswith('MemAvailable:'):
            return int(line.split()[1]) * 1024
    raise RuntimeError('MemAvailable unavailable; refusing launch')


def group_path():
    line = next(v for v in Path('/proc/self/cgroup').read_text().splitlines()
                if v.startswith('0::'))
    return Path('/sys/fs/cgroup') / line[3:].lstrip('/')


def sample(group):
    result = {'utc': datetime.now(timezone.utc).isoformat(),
              'hostAvailableBytes': available_memory()}
    for name in ('memory.current', 'memory.peak', 'memory.swap.current',
                 'memory.events', 'memory.pressure', 'io.stat', 'pids.current'):
        path = group / name
        if path.exists():
            result[name] = path.read_text().strip()
    return result


def stop(child):
    if child.poll() is not None:
        return
    os.killpg(child.pid, signal.SIGTERM)
    try:
        child.wait(timeout=5)
    except subprocess.TimeoutExpired:
        os.killpg(child.pid, signal.SIGKILL)
        child.wait()


def child_main(manifest):
    spec = json.loads(manifest.read_text())
    folder = manifest.parent
    os.sched_setaffinity(0, {spec['cpu']})
    group = group_path()
    actual = {name: (group / name).read_text().strip()
              for name in ('memory.max', 'memory.swap.max', 'pids.max')}
    actual.update(group=str(group), nice=os.getpriority(os.PRIO_PROCESS, 0),
                  affinity=sorted(os.sched_getaffinity(0)),
                  threads={name: os.environ.get(name) for name in THREADS})
    save(folder / 'effective-limits.json', actual)
    if not (actual['memory.max'] == str(MEMORY_MAX)
            and actual['memory.swap.max'] == '0'
            and actual['pids.max'] == '32' and actual['nice'] >= 15
            and actual['affinity'] == [spec['cpu']]
            and all(actual['threads'][name] == '1' for name in THREADS)):
        raise RuntimeError('actual resource limits differ; no workload started')
    if available_memory() < AVAILABLE_MIN:
        raise RuntimeError('less than 4 GiB available; no workload started')
    save(folder / 'limit-validation.json', {'verified': True, 'actual': actual})
    started = time.monotonic()
    reason = None
    low_memory = 0
    # The workload is the only scientific child. Its entire descendant tree
    # stays in this service cgroup, even if a descendant creates a new session.
    with (folder / 'resource-samples.jsonl').open('x') as telemetry:
        telemetry.write(json.dumps(sample(group)) + '\n'); telemetry.flush()
        child = subprocess.Popen(spec['command'], cwd=spec['cwd'],
                                 stdin=subprocess.DEVNULL, start_new_session=True)
        try:
            while True:
                try:
                    code = child.wait(timeout=2)
                    break
                except subprocess.TimeoutExpired:
                    observation = sample(group)
                    telemetry.write(json.dumps(observation) + '\n'); telemetry.flush()
                    low_memory = low_memory + 1 if observation['hostAvailableBytes'] < AVAILABLE_MIN else 0
                    if low_memory >= 2:
                        reason = 'host available memory below 4 GiB for two observations'
                    elif time.monotonic() - started >= spec['seconds']:
                        reason = 'wall-time limit'
                    if reason:
                        stop(child); code = 124; break
        finally:
            stop(child)
            telemetry.write(json.dumps(sample(group)) + '\n'); telemetry.flush()
    save(folder / 'child-outcome.json', {'exitCode': code, 'guardReason': reason,
         'wallSeconds': time.monotonic() - started, 'childPid': child.pid})
    return code if code >= 0 else 128 - code


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--log-directory', type=Path)
    parser.add_argument('--seconds', type=int, default=900)
    parser.add_argument('--child-manifest', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.child_manifest:
        return child_main(args.child_manifest)
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command or not args.log_directory or not 1 <= args.seconds <= 900:
        parser.error('fresh --log-directory, --seconds 1..900 and -- command are required')
    if PAUSE.exists():
        raise RuntimeError('S11c work is paused after the host freeze; see ' + str(PAUSE))
    folder = args.log_directory.resolve()
    folder.relative_to(STORE)
    STORE.mkdir(parents=True, exist_ok=True)
    lock = (STORE / 'resource-guard.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    # A killed launcher may leave its bounded service alive. Do not start a
    # second job merely because that launcher no longer holds the file lock.
    active = subprocess.check_output(
        ['systemctl', '--user', 'list-units', '--no-legend', '--plain',
         '--state=active,activating,deactivating', 's11c-guard-*.service'], text=True)
    if active.strip():
        raise RuntimeError('another guarded service remains active; refusing overlap')
    if available_memory() < AVAILABLE_MIN:
        raise RuntimeError('less than 4 GiB available; refusing launch')
    folder.mkdir(parents=True, exist_ok=False)
    unit = 's11c-guard-' + uuid.uuid4().hex[:12]
    spec = {'command': command, 'cwd': str(Path.cwd()), 'unit': unit,
            'seconds': args.seconds, 'cpu': max(os.sched_getaffinity(0)),
            'memoryMax': MEMORY_MAX, 'swapMax': 0, 'tasksMax': 32,
            'minimumHostAvailableBytes': AVAILABLE_MIN,
            'startedUtc': datetime.now(timezone.utc).isoformat()}
    save(folder / 'invocation.json', spec)
    launch = ['systemd-run', '--user', '--quiet', '--wait', '--pipe', '--collect',
              '--unit=' + unit, '--property=MemoryMax=' + str(MEMORY_MAX),
              '--property=MemorySwapMax=0', '--property=TasksMax=32',
              '--property=RuntimeMaxSec=' + str(args.seconds + 10),
              '--property=TimeoutStopSec=5', '--property=KillMode=control-group',
              '--property=OOMPolicy=kill', '--property=Nice=15',
              '--property=IOSchedulingClass=idle']
    launch += ['--setenv=' + name + '=1' for name in THREADS]
    launch += [sys.executable, str(Path(__file__).resolve()),
               '--child-manifest', str(folder / 'invocation.json')]
    started = time.monotonic()
    with (folder / 'stdout').open('xb') as stdout, (folder / 'stderr').open('xb') as stderr:
        try:
            result = subprocess.run(launch, stdout=stdout, stderr=stderr, stdin=subprocess.DEVNULL)
        except BaseException:
            subprocess.run(['systemctl', '--user', 'stop', unit],
                           stdout=stderr, stderr=stderr, check=False)
            raise
    outcome = {'exitCode': result.returncode, 'wallSeconds': time.monotonic()-started,
               'unit': unit, 'stderrBytes': (folder / 'stderr').stat().st_size,
               'limitsVerified': (folder / 'limit-validation.json').exists(),
               'childOutcome': json.loads((folder / 'child-outcome.json').read_text())
               if (folder / 'child-outcome.json').exists() else None}
    # A missing completion receipt never becomes successful scientific evidence.
    code = result.returncode or (0 if outcome['childOutcome'] and
                                 outcome['childOutcome']['exitCode'] == 0 else 1)
    outcome['exitCode'] = code
    save(folder / 'outcome.json', outcome)
    print(json.dumps(outcome, indent=2))
    return code


if __name__ == '__main__':
    sys.exit(main())
