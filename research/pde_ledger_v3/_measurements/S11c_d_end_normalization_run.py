#!/usr/bin/env python3
"""Run one normalization constructor or validator, keeping durable live logs."""
import argparse
from datetime import datetime, timezone
import fcntl
import json
import os
import resource
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
STORE = ROOT.parents[1]/'_scratch/s11c'


def prerequisite_lock_mode(parallel):
    if not parallel:
        return fcntl.LOCK_EX
    # Only the resource-verified pooled guard grants this opt-in. The historical
    # producer still takes EX, while complete-prerequisite readers can share SH.
    manifest = Path(os.environ['S11C_POOLED_GUARD_MANIFEST']).resolve()
    manifest.relative_to(STORE)
    spec = json.loads(manifest.read_text())
    validation = json.loads((manifest.parent / 'limit-validation.json').read_text())
    cgroup = next(line[3:] for line in Path('/proc/self/cgroup').read_text().splitlines()
                  if line.startswith('0::'))
    group = str(Path('/sys/fs/cgroup') / cgroup.lstrip('/'))
    if (not spec.get('pool') or not validation.get('verified')
            or spec.get('priorityPolicy') != 'desktop-managed'
            or validation['actual']['group'] != group
            or sorted(os.sched_getaffinity(0)) != [spec['cpu']]
            or tuple(resource.getrlimit(resource.RLIMIT_AS)) != (spec['memoryMax'],) * 2):
        raise ValueError('shared prerequisite reads require verified pooled containment')
    return fcntl.LOCK_SH


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-root',type=Path,required=True)
    parser.add_argument('--stage',required=True)
    parser.add_argument('--parallel-prerequisite-read', action='store_true',
                        help='verified pooled guard only; share completed historical prerequisite read lock')
    parser.add_argument('command',nargs=argparse.REMAINDER)
    args = parser.parse_args()
    base = args.run_root.resolve();base.relative_to(STORE)
    base.mkdir(parents=True,exist_ok=True)
    command = args.command[1:] if args.command[:1]==['--'] else args.command
    if not command or Path(command[0]).name not in ('python','python3'):
        raise ValueError('Python stage command required')
    lock = (base/'normalization.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    # The preceding queue owns a separate lock at its historical durable root.
    previous = STORE/'s11c-thickness-coordinate-20260914'
    prior_lock = (previous/'endpoint_controller.lock').open('a')
    fcntl.flock(prior_lock,prerequisite_lock_mode(args.parallel_prerequisite_read)|fcntl.LOCK_NB)
    progress = json.loads((previous/'endpoint_continuation/progress.json').read_text())
    if progress[-1]['event'] != 'endpoint_plan_complete':
        raise ValueError('endpoint prerequisite queue has not completed')
    stdout_path = base/(args.stage+'.stdout')
    stderr_path = base/(args.stage+'.stderr')
    record = {'stage':args.stage,'command':command,
        'prerequisiteLock': 'shared-read' if args.parallel_prerequisite_read else 'exclusive',
        'stdout':str(stdout_path),
        'stderr':str(stderr_path),'startedUtc':datetime.now(timezone.utc).isoformat()}
    pointer = base/'active.json'
    def save():
        temporary = pointer.with_suffix('.new')
        temporary.write_text(json.dumps(record,indent=2)+'\n');temporary.replace(pointer)
    started = time.monotonic()
    with stdout_path.open('xb') as stdout, stderr_path.open('xb') as stderr:
        child = subprocess.Popen(command,cwd=ROOT,stdout=stdout,stderr=stderr)
        record.update({'childPid':child.pid,'status':'running'});save()
        code = child.wait()
    record.update({'status':'completed' if code==0 else 'failed','exitCode':code,
        'wallSeconds':time.monotonic()-started,'stderrBytes':stderr_path.stat().st_size,
        'finishedUtc':datetime.now(timezone.utc).isoformat()})
    # Persist the child outcome even if subsequent log bookkeeping fails.
    save()
    if '--run-directory' in command and any(Path(word).name == 'S11c_d_end_normalization_check.py' for word in command[1:]):
        destination = Path(command[command.index('--run-directory')+1]).resolve()
        destination.relative_to(base)
        if destination.exists():
            # Keep the original live logs as well; no symlink into /tmp.
            for source, name in ((stdout_path,'full.out'),(stderr_path,'stderr.txt')):
                target = destination/name
                if target.exists():raise ValueError('constructor log target already exists')
                target.hardlink_to(source)
    save()
    (base/(args.stage+'.invocation.json')).write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))
    return code


if __name__ == '__main__':
    sys.exit(main())
