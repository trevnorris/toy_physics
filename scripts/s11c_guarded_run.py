#!/usr/bin/env python3
"""Run one S11c job with memory containment and no elapsed-time cutoff.

Requires a host shell and a systemd user manager. Fails closed if the actual
cgroup limits cannot be verified. No scientific imports or automatic retries.
Standing user policy: no wall-clock, worker or inactivity deadlines.
The older runner remains frozen for existing source pins.
Defaults remain 2 GiB and 32 tasks. Opt-in --pool jobs share one global 16 GiB
budget across all pool names and at most six distinct physical cores. Legacy
exclusive launches remain unchanged. Native threads remain one. No workload
starts until actual cgroup memory, native address-space, inherited process
affinity and duration enforcement is verified; pooled priority is desktop-managed.

CPU allocation uses systemd CPUAffinity plus inherited process affinity, not a
hard cgroup cpuset. The guard verifies both at startup; it does not authorize
workers to change their allocated affinity. No delegated cpuset controller is
required. The normalization supervisor needs --parallel-prerequisite-read for
pooled overlap (its per-run lock remains exclusive). A live unaccounted service
or stale reservation fails closed. Reservations are never time leases and are
never reclaimed automatically after launcher death. Explicit reconciliation
must first establish terminal systemd state, no queued start and a dead owner;
keep the old registry and run receipts when correcting such stale metadata.
This tool does not authorize sharing mutable scientific output directories.
"""
import argparse
from datetime import datetime, timezone
import fcntl
import json
import os
import re
import resource
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
POOL_MEMORY_MAX = 16 * 1024**3
POOL_MAX_JOBS = 6
POOL_DIRECTORY = STORE / 'resource-pool'
THREADS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS')


def resource_limits(memory_gib=2, tasks_max=32, pool=None):
    allowed = range(1, 17) if pool else (2, 4)
    if memory_gib not in allowed or tasks_max not in (32, 64):
        raise ValueError('unsupported resource limits')
    return {'memoryMax': memory_gib * 1024**3, 'swapMax': 0,
            'tasksMax': tasks_max}


def limits_match(actual, spec):
    # Legacy jobs retain their exact profile. Parallel jobs have an explicit,
    # separately accounted profile and desktop-managed priority.
    pooled = spec.get('pool') is not None
    profile = (isinstance(spec['memoryMax'], int)
               and spec['memoryMax'] % (1024**3) == 0
               and 1024**3 <= spec['memoryMax'] <= POOL_MEMORY_MAX
               and spec.get('poolMemoryMax') == POOL_MEMORY_MAX
               and spec.get('priorityPolicy') == 'desktop-managed'
               and actual.get('nativeAddressSpaceBytes') == [spec['memoryMax']] * 2
               and spec.get('affinityEnforcement') == 'systemd-and-inherited-process'
               and actual.get('systemd.CPUAffinity') == str(spec['cpu'])) if pooled else (
                   spec['memoryMax'] in (MEMORY_MAX, 2 * MEMORY_MAX)
                   and actual['nice'] >= 15)
    return (profile and spec['tasksMax'] in (32, 64) and spec['swapMax'] == 0
            and actual['memory.max'] == str(spec['memoryMax'])
            and actual['memory.swap.max'] == '0'
            and actual['pids.max'] == str(spec['tasksMax'])
            and actual['affinity'] == [spec['cpu']]
            and all(actual['threads'][name] == '1' for name in THREADS))


def acquire_execution_lock(path, pooled=False):
    lock = path.open('a')
    try:
        fcntl.flock(lock, (fcntl.LOCK_SH if pooled else fcntl.LOCK_EX) | fcntl.LOCK_NB)
    except BaseException:
        lock.close()
        raise
    return lock


def active_units():
    text = subprocess.check_output(
        ['systemctl', '--user', 'list-units', '--no-legend', '--plain',
         '--state=active,activating,deactivating', 's11c-guard-*.service'], text=True)
    units = [line.split()[0] for line in text.splitlines() if line.strip()]
    if any(not re.fullmatch(r's11c-guard-[a-zA-Z0-9-]+\.service', unit) for unit in units):
        raise RuntimeError('unrecognized service listing; refusing admission')
    return units


def process_identity(pid):
    try:
        stat = Path('/proc', str(pid), 'stat').read_text()
        return {'pid': pid, 'startTicks': stat.rsplit(')', 1)[1].split()[19],
                'bootId': Path('/proc/sys/kernel/random/boot_id').read_text().strip()}
    except FileNotFoundError:
        return None


def pool_cpus():
    # One logical CPU per physical core, so admitted workers do not compete
    # for sibling threads while idle physical cores are available.
    physical = {}
    for cpu in sorted(os.sched_getaffinity(0), reverse=True):
        topology = Path('/sys/devices/system/cpu', 'cpu' + str(cpu), 'topology')
        key = tuple((topology / name).read_text().strip()
                    for name in ('physical_package_id', 'core_id'))
        physical.setdefault(key, cpu)
    return list(physical.values())[:POOL_MAX_JOBS]


def service_resources(unit):
    raw = subprocess.check_output(['systemctl', '--user', 'show', unit,
        '--property=ControlGroup', '--property=RuntimeMaxUSec', '--property=Restart',
        '--property=CPUAffinity'], text=True)
    info = dict(line.split('=', 1) for line in raw.splitlines() if '=' in line)
    if info.get('RuntimeMaxUSec') != 'infinity' or info.get('Restart') != 'no':
        raise RuntimeError('live service duration/restart mismatch: ' + unit)
    location = info.get('ControlGroup', '')
    if not location.startswith('/') or '..' in Path(location).parts:
        raise RuntimeError('live service has no verifiable cgroup: ' + unit)
    group = Path('/sys/fs/cgroup') / location.lstrip('/')
    actual = {name: (group / name).read_text().strip() for name in
              ('memory.max', 'memory.swap.max', 'memory.current', 'pids.max')}
    actual['systemd.CPUAffinity'] = info.get('CPUAffinity')
    return actual


def read_reservations(directory):
    path = directory / 'reservations.json'
    if not path.exists():
        return {}
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise RuntimeError('malformed resource reservation registry')
    return value


def write_reservations(directory, value):
    path = directory / 'reservations.json'
    temporary = directory / 'reservations.json.new'
    with temporary.open('w') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)
    descriptor = os.open(directory, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def admission_plan(reservations, live, requested, cpus, available,
                   identity=process_identity):
    # A reservation is written before requesting systemd launch. Never silently
    # reclaim a dead launcher's reservation: a queued start could still exist.
    # Unknown active services likewise block, even if a launcher lock vanished.
    if set(live) - set(reservations):
        raise RuntimeError('unaccounted guarded service; refusing pooled overlap')
    occupied = set()
    reserved = 0
    headroom = 0
    for unit, record in reservations.items():
        cap = record['memoryMax']
        cpu = record['cpu']
        if (not isinstance(cap, int) or cap % (1024**3) or not 1024**3 <= cap <= POOL_MEMORY_MAX
                or cpu in occupied or record.get('poolMemoryMax') != POOL_MEMORY_MAX
                or record.get('tasksMax') not in (32, 64)):
            raise RuntimeError('invalid or overlapping saved reservation')
        occupied.add(cpu)
        reserved += cap
        used = 0
        if unit in live:
            actual = live[unit]
            if (actual['memory.max'] != str(cap) or actual['memory.swap.max'] != '0'
                    or actual['pids.max'] != str(record['tasksMax'])
                    or actual['systemd.CPUAffinity'] != str(cpu)):
                raise RuntimeError('live service does not match reservation: ' + unit)
            used = int(actual['memory.current'])
            if used < 0 or used > cap:
                raise RuntimeError('live service memory accounting outside cap')
        elif identity(record['owner']['pid']) != record['owner']:
            raise RuntimeError('unconfirmed stale reservation; explicit reconciliation required: ' + unit)
        headroom += cap - used
    if (requested < 1024**3 or requested > POOL_MEMORY_MAX or requested % (1024**3)
            or reserved + requested > POOL_MEMORY_MAX):
        raise RuntimeError('aggregate 16 GiB reservation budget exceeded')
    if len(reservations) >= POOL_MAX_JOBS:
        raise RuntimeError('aggregate worker budget exhausted')
    free = [cpu for cpu in cpus if cpu not in occupied]
    if not free:
        raise RuntimeError('no unreserved physical CPU available')
    if available < AVAILABLE_MIN + headroom + requested:
        raise RuntimeError('host reserve plus outstanding memory commitments unavailable')
    return {'cpu': free[0], 'reservedBeforeBytes': reserved,
            'outstandingHeadroomBytes': headroom,
            'hostAvailableBytes': available, 'poolMemoryMax': POOL_MEMORY_MAX,
            'poolMaxJobs': POOL_MAX_JOBS, 'activeUnits': sorted(live)}


def reserve_pool(spec, directory=POOL_DIRECTORY):
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / 'admission.lock').open('a') as admission:
        fcntl.flock(admission, fcntl.LOCK_EX)
        records = read_reservations(directory)
        live = {unit: service_resources(unit) for unit in active_units()}
        plan = admission_plan(records, live, spec['memoryMax'], pool_cpus(), available_memory())
        spec.update(plan)
        spec['owner'] = process_identity(os.getpid())
        if not spec['owner']:
            raise RuntimeError('cannot establish launcher identity')
        records[spec['unit'] + '.service'] = {key: spec[key] for key in
            ('pool', 'poolMemoryMax', 'memoryMax', 'tasksMax', 'cpu', 'owner', 'logDirectory')}
        write_reservations(directory, records)
        return plan


def release_pool(spec, directory=POOL_DIRECTORY):
    with (directory / 'admission.lock').open('a') as admission:
        fcntl.flock(admission, fcntl.LOCK_EX)
        unit = spec['unit'] + '.service'
        if unit in active_units():
            raise RuntimeError('service still active; reservation retained: ' + unit)
        records = read_reservations(directory)
        record = records.get(unit)
        if record is None or record['owner'] != spec['owner']:
            raise RuntimeError('reservation ownership mismatch; refusing release')
        del records[unit]
        write_reservations(directory, records)


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
    properties = ['--property=RuntimeMaxUSec', '--property=Restart']
    if spec.get('pool') is not None:
        properties.append('--property=CPUAffinity')
    raw = subprocess.check_output(['systemctl', '--user', 'show', spec['unit'],
                                   *properties], text=True)
    unit_properties = dict(line.split('=', 1) for line in raw.splitlines() if '=' in line)
    duration = {key: unit_properties.get(key) for key in ('RuntimeMaxUSec', 'Restart')}
    if duration != {'RuntimeMaxUSec': 'infinity', 'Restart': 'no'}:
        raise RuntimeError('unlimited runtime/restart verification failed; no workload started')
    save(folder / 'duration-validation.json', {'verified': True, 'actual': duration,
        'wallDeadlineSeconds': None, 'inactivityDeadlineSeconds': None})
    if spec.get('pool') is not None:
        resource.setrlimit(resource.RLIMIT_AS, (spec['memoryMax'], spec['memoryMax']))
    os.sched_setaffinity(0, {spec['cpu']})
    group = group_path()
    actual = {name: (group / name).read_text().strip()
              for name in ('memory.max', 'memory.swap.max', 'pids.max')}
    actual.update(group=str(group), nice=os.getpriority(os.PRIO_PROCESS, 0),
                  affinity=sorted(os.sched_getaffinity(0)),
                  threads={name: os.environ.get(name) for name in THREADS})
    if spec.get('pool') is not None:
        actual['nativeAddressSpaceBytes'] = list(resource.getrlimit(resource.RLIMIT_AS))
        actual['systemd.CPUAffinity'] = unit_properties.get('CPUAffinity')
        actual['affinityEnforcement'] = 'systemd-and-inherited-process'
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
        environment = os.environ.copy()
        environment.pop('S11C_POOLED_GUARD_MANIFEST', None)
        if spec.get('pool') is not None:
            environment['S11C_POOLED_GUARD_MANIFEST'] = str(manifest.resolve())
        child = subprocess.Popen(spec['command'], cwd=spec['cwd'], env=environment,
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
    parser.add_argument('--seconds', type=int, default=0,
                        help='legacy compatibility argument; recorded but never imposes a deadline')
    parser.add_argument('--memory-gib', type=int, default=2,
                        help='whole-job cap: legacy 2/4 GiB; pooled 1 through 16 GiB')
    parser.add_argument('--pool', help='opt in to shared 16 GiB / six-worker admission; named purpose')
    parser.add_argument('--tasks-max', type=int, choices=(32, 64), default=32,
                        help='whole-job process/thread cap; 64 requires explicit job authorization')
    parser.add_argument('--child-manifest', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.child_manifest:
        return child_main(args.child_manifest)
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command or not args.log_directory:
        parser.error('fresh --log-directory and -- command are required')
    if PAUSE.exists():
        raise RuntimeError('S11c work is paused after the host freeze; see ' + str(PAUSE))
    if args.pool and not re.fullmatch(r'[a-zA-Z0-9][a-zA-Z0-9_-]{0,63}', args.pool):
        parser.error('--pool must be a short alphanumeric purpose name')
    profile = resource_limits(args.memory_gib, args.tasks_max, args.pool)
    folder = args.log_directory.resolve()
    folder.relative_to(STORE)
    STORE.mkdir(parents=True, exist_ok=True)
    lock = acquire_execution_lock(STORE / 'resource-guard.lock', pooled=bool(args.pool))
    # Legacy and pooled launchers cannot overlap. A surviving service/reservation
    # also blocks a legacy launcher whose predecessor lost its lock.
    if not args.pool:
        if active_units() or read_reservations(POOL_DIRECTORY):
            raise RuntimeError('another guarded service or reservation remains; refusing overlap')
    if available_memory() < AVAILABLE_MIN:
        raise RuntimeError('less than 4 GiB available; refusing launch')
    folder.mkdir(parents=True, exist_ok=False)
    unit = 's11c-guard-' + uuid.uuid4().hex[:12]
    spec = {'command': command, 'cwd': str(Path.cwd()), 'unit': unit,
            'seconds': 0, 'requestedSecondsArgumentIgnored': args.seconds,
            'wallDeadlineSeconds': None, 'inactivityDeadlineSeconds': None,
            'cpu': max(os.sched_getaffinity(0)),
            **profile,
            'minimumHostAvailableBytes': AVAILABLE_MIN,
            'startedUtc': datetime.now(timezone.utc).isoformat()}
    if args.pool:
        spec.update(pool=args.pool, priorityPolicy='desktop-managed', logDirectory=str(folder),
                    affinityEnforcement='systemd-and-inherited-process')
        reserve_pool(spec)
    save(folder / 'invocation.json', spec)
    launch = ['systemd-run', '--user', '--quiet', '--wait', '--pipe', '--collect',
              '--unit=' + unit, '--property=MemoryMax=' + str(spec['memoryMax']),
              '--property=MemorySwapMax=0', '--property=TasksMax=' + str(spec['tasksMax']),
              '--property=RuntimeMaxSec=infinity', '--property=Restart=no',
              '--property=TimeoutStopSec=5', '--property=KillMode=control-group',
              '--property=OOMPolicy=kill', '--property=IOSchedulingClass=idle']
    if args.pool:
        launch += ['--property=CPUAffinity=' + str(spec['cpu'])]
    else:
        launch += ['--property=Nice=15']
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
            if args.pool:
                release_pool(spec)
            raise
    if args.pool:
        release_pool(spec)
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
