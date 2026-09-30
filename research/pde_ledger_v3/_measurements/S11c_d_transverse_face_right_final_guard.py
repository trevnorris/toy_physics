#!/usr/bin/env python3
"""Task-local guard for one explicitly approved progress-dependent continuation.

The shared runner is unchanged. Containment, locks, resource sampling and
fail-closed controls are retained; exact duration approval is required.
No total duration cap; stop on one hour without a newly saved result, with
60 seconds for native failure bookkeeping. CPU activity is not progress.

Run one S11c job in a bounded user service, with durable resource telemetry.

Requires a host shell and a systemd user manager. Fails closed if the actual
cgroup limits cannot be verified. No scientific imports or automatic retries.
This job alone uses the expressly approved 16 GiB and desktop-managed priority.
Zero swap, one CPU, 32 tasks, one thread and all lock/host protections remain.
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
MEMORY_MAX = 16 * 1024**3
AVAILABLE_MIN = 4 * 1024**3
THREADS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS')


def resource_limits(memory_gib=16, tasks_max=32):
    if memory_gib != 16 or tasks_max != 32:
        raise ValueError('unsupported resource limits')
    return {'memoryMax': memory_gib * 1024**3, 'swapMax': 0,
            'tasksMax': tasks_max}


def limits_match(actual, spec):
    # Validate the requested bounds as well as their enforcement. An altered
    # manifest must not turn a missing or unlimited cap into an accepted job.
    return (spec['memoryMax'] == MEMORY_MAX
            and spec['tasksMax'] == 32 and spec['swapMax'] == 0
            and actual['memory.max'] == str(spec['memoryMax'])
            and actual['memory.swap.max'] == '0'
            and actual['pids.max'] == str(spec['tasksMax'])
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



def validate_duration(seconds, memory_gib, tasks_max, command, folder):
    """Fail closed on any duration, command or scope outside this approved job."""
    gate_path=ROOT/'research/pde_ledger_v3/_measurements/S11c_d_transverse_face_right_final_gate.json'
    gate=json.loads(gate_path.read_text())
    digest=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    if not (gate['status']=='READY_FOR_ONE_FINAL_RIGHT_END_CONTINUATION'
            and gate['progressDependentRuntimeExplicitlyApproved'] is True and gate['progressStallSeconds']==3600
            and seconds==gate['seconds']==0 and gate['nativeSeconds']==0
            and memory_gib==gate['memoryGiB']==16 and gate['memoryBytes']==16*1024**3 and tasks_max==gate['tasksMax']==32
            and (gate['swapMax'],gate['cpuCount'],gate['nativeThreads'])==(0,1,1)
            and gate['priorityPolicy']=='DESKTOP_MANAGED' and gate['scienceJobOrdinal']==7
            and command==gate['guardedCommand']
            and str(Path(folder).resolve())==str(Path(gate['runRoot'])/'resource-guard')):
        raise RuntimeError('request differs from explicitly approved no-deadline job')
    approval=gate['runtimeApproval']
    if digest(approval['path'])!=approval['sha256']:
        raise RuntimeError('runtime approval receipt changed')
    authorized=json.loads(Path(approval['path']).read_text())
    if not (authorized['status']=='USER_DIRECTED_FINAL_RIGHT_END_CONTINUATION'
            and authorized['userReply']=="So be safe and give it 16 GB. I'll make sure nothing else runs on the box until this is complete."
            and authorized['memoryBytes']==16*1024**3
            and authorized['priorityPolicy']=='DESKTOP_MANAGED'
            and authorized['priorityReply']=='Resume with desktop-managed priority'
            and authorized['scientificExecutionOrdinal']==7
            and authorized['runtimeReply']=="I give permission for running them longer than the normal limit. allow them to run as long as necessary as long as they are making progress"
            and authorized['seconds']==seconds and authorized['nativeSeconds']==0
            and authorized['progressStallSeconds']==3600
            and authorized['maximumNewScientificExecutions']==1
            and authorized['finalAttempt'] is True and authorized['maximumReviewRounds']==1
            and gate['finalAttempt'] is True):
        raise RuntimeError('explicit no-deadline approval unavailable')
    for p,h in ((Path(__file__),gate['guardSha256']),
                (ROOT/'scripts/s11c_guarded_run.py',gate['sharedGuardSha256']),
                (Path(gate['workerPath']),gate['workerSha256']),
                (Path(gate['inputManifestPath']),gate['inputManifestSha256'])):
        if digest(p)!=h:
            raise RuntimeError('approved no-deadline source pin changed: '+str(p))
    return {'gatePath':str(gate_path),'gateSha256':digest(gate_path),
            'approvalPath':approval['path'],'approvalSha256':approval['sha256'],
            'seconds':seconds,'nativeSeconds':0,'sharedGuardUnchanged':True,
            'progressStallSeconds':3600,'guardBookkeepingGraceSeconds':60,
            'runRoot':gate['runRoot']}


def verify_service_duration(spec):
    raw=subprocess.check_output(['systemctl','--user','show',spec['unit'],
        '--property=RuntimeMaxUSec','--property=Restart'],text=True)
    actual=dict(line.split('=',1) for line in raw.splitlines() if '=' in line)
    if actual!={'RuntimeMaxUSec':'infinity','Restart':'no'}:
        raise RuntimeError('service duration or restart differs; no workload started')
    return {'verified':True,'actual':actual,'wallDeadlineSeconds':None,'automaticRetry':False}

def child_main(manifest):
    spec = json.loads(manifest.read_text())
    duration = validate_duration(spec['seconds'], spec['memoryMax']//1024**3,
        spec['tasksMax'], spec['command'], manifest.parent)
    if duration != spec['timeAuthorization']:
        raise RuntimeError('parent/child duration authorization differs')
    folder = manifest.parent
    save(folder / 'duration-validation.json', verify_service_duration(spec))
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
    progress_sequence=0
    progress_seen=set()
    last_progress=started
    progress_root=Path(duration['runRoot'])
    progress_path=progress_root/'progress.json'
    def read_progress():
        if not progress_path.exists(): return None
        event=json.loads(progress_path.read_text())
        if not isinstance(event['sequence'],int) or event['sequence']<1:
            raise RuntimeError('invalid progress sequence')
        receipt=event['receipt']; path=Path(receipt['path']).resolve(strict=True)
        path.relative_to(progress_root/'complete')
        raw=path.read_bytes()
        if len(raw)!=receipt['bytes'] or hashlib.sha256(raw).hexdigest()!=receipt['sha256']:
            raise RuntimeError('progress receipt changed')
        if event['kind'] not in ('PRESERVED_FIXED_POINT_RESULTS','RESTORED_PRIOR_COMPLETE_SUBOPERATION',
                                 'RESTORED_PRIOR_COMPLETE_RETURN','NEW_COMPLETE_SUBOPERATION','NEW_COMPLETE_OPERATION'):
            raise RuntimeError('unsupported progress event')
        return event
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
                    try:
                        progress=read_progress()
                        if progress is not None:
                            if progress['sequence']<progress_sequence:
                                raise RuntimeError('progress sequence decreased')
                            if progress['sequence']>progress_sequence:
                                receipt_path=progress['receipt']['path']
                                if receipt_path in progress_seen:
                                    raise RuntimeError('old result replayed as new progress')
                                progress_seen.add(receipt_path)
                                progress_sequence=progress['sequence'];last_progress=time.monotonic()
                                save(folder/'last-progress.json',dict(event=progress,observedUtc=datetime.now(timezone.utc).isoformat()))
                    except (OSError,ValueError,KeyError,RuntimeError) as error:
                        reason='invalid progress evidence: '+str(error)
                    if not reason and time.monotonic()-last_progress >= 3660:
                        reason='no newly saved result for 3600 seconds plus 60 seconds bookkeeping grace'
                    if low_memory >= 2:
                        reason = 'host available memory below 4 GiB for two observations'
                    elif spec['seconds'] and time.monotonic() - started >= spec['seconds']:
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
    parser.add_argument('--memory-gib', type=int, choices=(16,), default=16,
                        help='whole-job memory cap: explicitly authorized 16 GiB only')
    parser.add_argument('--tasks-max', type=int, choices=(32,), default=32,
                        help='whole-job process/thread cap: 32')
    parser.add_argument('--child-manifest', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.child_manifest:
        return child_main(args.child_manifest)
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command or not args.log_directory or args.seconds != 0:
        parser.error('this task-local runner requires the explicitly approved --seconds 0 (no deadline) job')
    if PAUSE.exists():
        raise RuntimeError('S11c work is paused after the host freeze; see ' + str(PAUSE))
    folder = args.log_directory.resolve()
    folder.relative_to(STORE)
    time_authorization = validate_duration(args.seconds, args.memory_gib, args.tasks_max, command, folder)
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
            'startedUtc': datetime.now(timezone.utc).isoformat(), 'timeAuthorization': time_authorization}
    save(folder / 'invocation.json', spec)
    launch = ['systemd-run', '--user', '--quiet', '--wait', '--pipe', '--collect',
              '--unit=' + unit, '--property=MemoryMax=' + str(spec['memoryMax']),
              '--property=MemorySwapMax=0', '--property=TasksMax=' + str(spec['tasksMax']),
              '--property=RuntimeMaxSec=infinity', '--property=Restart=no',
              '--property=TimeoutStopSec=5', '--property=KillMode=control-group',
              '--property=OOMPolicy=kill',
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
