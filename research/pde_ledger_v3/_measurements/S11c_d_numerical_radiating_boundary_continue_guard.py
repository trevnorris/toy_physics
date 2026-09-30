#!/usr/bin/env python3
"""Run one S11c job in a bounded user service, with durable resource telemetry.

Requires a host shell and a systemd user manager. Fails closed if the actual
cgroup limits cannot be verified. No scientific imports or automatic retries.
Task-local no-deadline continuation. Shared guard unchanged.
Defaults remain 2 GiB and 32 tasks. Nondefault limits require explicit user
authorization for the particular job; they do not change other jobs' defaults.
"""
import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import uuid

ROOT = Path('/var/projects/toy_physics')
STORE = ROOT / '_scratch/s11c'
PAUSE = STORE / 'PAUSED_HOST_FREEZE.json'
MEMORY_MAX = 2 * 1024**3
AVAILABLE_MIN = 4 * 1024**3
THREADS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS')


def resource_limits(memory_gib=2, tasks_max=32):
    if memory_gib not in (2, 4) or tasks_max not in (32, 64):
        raise ValueError('unsupported resource limits')
    return {'memoryMax': memory_gib * 1024**3, 'swapMax': 0,
            'tasksMax': tasks_max}


def limits_match(actual, spec):
    # Validate the requested bounds as well as their enforcement. An altered
    # manifest must not turn a missing or unlimited cap into an accepted job.
    return (spec['memoryMax'] in (MEMORY_MAX, 2 * MEMORY_MAX)
            and spec['tasksMax'] in (32, 64) and spec['swapMax'] == 0
            and actual['memory.max'] == str(spec['memoryMax'])
            and actual['memory.swap.max'] == '0'
            and actual['pids.max'] == str(spec['tasksMax'])
            and actual['nice'] >= 15
            and actual['affinity'] == [spec['cpu']]
            and all(actual['threads'][name] == '1' for name in THREADS))


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


def validate_duration(command,folder):
    gate_path=ROOT/'research/pde_ledger_v3/_measurements/S11c_d_numerical_radiating_boundary_continue_gate.json'
    gate=json.loads(gate_path.read_text())
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    if not (gate['status']=='READY_FOR_USER_AUTHORIZED_NO_DEADLINE_BOUNDARY_CONTINUATION'
            and gate['proceedAuthority']=='EXPLICIT_USER_NO_TIME_LIMIT_SAVED_RETURN_CONTINUATION'
            and gate['wallDeadlineSeconds'] is None and gate['nativeDeadlineSeconds'] is None
            and gate['inactivityDeadlineSeconds'] is None
            and command==gate['guardedCommand']
            and str(folder)==str(Path(gate['runRoot'])/'resource-guard')):
        raise RuntimeError('exact no-deadline job authority mismatch')
    for p,h in gate['sourcePins'].items():
        if sha(p)!=h:raise RuntimeError('approved source/helper changed: '+p)
    if sha(gate['inputManifestPath'])!=gate['manifestSha256'] or sha(__file__)!=gate['guardSha256']:
        raise RuntimeError('duration guard/manifest pins differ')
    approval=gate['runtimeApproval'];authority=json.loads(Path(approval['path']).read_text())
    if sha(approval['path'])!=approval['sha256'] or authority['literalUserInstruction']!="Stop it with the time limits. We're losing work because of that crap. Seriously. Yes run a continuation":
        raise RuntimeError('explicit user duration approval unavailable')
    return {'gatePath':str(gate_path),'gateSha256':sha(gate_path),'approval':approval,'wallDeadlineSeconds':None,'nativeDeadlineSeconds':None,'inactivityDeadlineSeconds':None,'sharedGuardUnchanged':True}


def child_main(manifest):
    spec = json.loads(manifest.read_text())
    folder = manifest.parent
    if spec['memoryMax']!=2*1024**3 or spec['tasksMax']!=32 or spec['seconds']!=0 or validate_duration(spec['command'],folder)!=spec['timeAuthorization']:
        raise RuntimeError('parent/child duration authorization differs')
    actual_duration=dict(line.split('=',1) for line in subprocess.check_output(['systemctl','--user','show',spec['unit'],'--property=RuntimeMaxUSec','--property=Restart'],text=True).splitlines() if '=' in line)
    if actual_duration!={'RuntimeMaxUSec':'infinity','Restart':'no'}:
        raise RuntimeError('unlimited service duration/restart mismatch; no workload started')
    save(folder/'duration-validation.json',{'verified':True,'actual':actual_duration,'wallDeadlineSeconds':None,'inactivityDeadlineSeconds':None})
    os.sched_setaffinity(0, {spec['cpu']})
    group = group_path()
    actual = {name: (group / name).read_text().strip()
              for name in ('memory.max', 'memory.swap.max', 'pids.max')}
    actual.update(group=str(group), nice=os.getpriority(os.PRIO_PROCESS, 0),
                  affinity=sorted(os.sched_getaffinity(0)),
                  threads={name: os.environ.get(name) for name in THREADS})
    save(folder / 'effective-limits.json', actual)
    if not limits_match(actual, spec):
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
    parser.add_argument('--seconds', type=int, default=0)
    parser.add_argument('--memory-gib', type=int, choices=(2, 4), default=2,
                        help='whole-job memory cap; 4 requires explicit job authorization')
    parser.add_argument('--tasks-max', type=int, choices=(32, 64), default=32,
                        help='whole-job process/thread cap; 64 requires explicit job authorization')
    parser.add_argument('--child-manifest', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.child_manifest:
        return child_main(args.child_manifest)
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command or not args.log_directory or args.seconds != 0 or args.memory_gib!=2 or args.tasks_max!=32:
        parser.error('this approved job requires --seconds 0 and the exact guarded command')
    if PAUSE.exists():
        raise RuntimeError('S11c work is paused after the host freeze; see ' + str(PAUSE))
    folder = args.log_directory.resolve()
    folder.relative_to(STORE)
    time_authorization = validate_duration(command, folder)
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
            **resource_limits(args.memory_gib, args.tasks_max),
            'minimumHostAvailableBytes': AVAILABLE_MIN,
            'startedUtc': datetime.now(timezone.utc).isoformat(), 'timeAuthorization':time_authorization}
    save(folder / 'invocation.json', spec)
    launch = ['systemd-run', '--user', '--quiet', '--wait', '--pipe', '--collect',
              '--unit=' + unit, '--property=MemoryMax=' + str(spec['memoryMax']),
              '--property=MemorySwapMax=0', '--property=TasksMax=' + str(spec['tasksMax']),
              '--property=RuntimeMaxSec=infinity', '--property=Restart=no',
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
