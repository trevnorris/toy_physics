#!/usr/bin/env python3
"""Silently wait for an owned local job; queue only completion/error events.

Linux pidfds bind the wait to the original process, not a reusable PID. No
model call or scheduled Codex task runs while waiting. Notifications use the
installed CLI's documented `codex queue --thread ... --message ...` command.
"""
import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import select
import shutil
import subprocess
import uuid


def save(path, value):
    temporary = path.with_suffix('.new')
    with temporary.open('w') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n'); stream.flush(); os.fsync(stream.fileno())
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--pid', type=int, required=True, help='Owned job supervisor PID')
    parser.add_argument('--thread', type=uuid.UUID, required=True)
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--message-file', type=Path, required=True)
    parser.add_argument('--error-log', type=Path, help='Only for a log where nonempty means an issue')
    parser.add_argument('--codex-bin', default=shutil.which('codex'))
    args = parser.parse_args()
    base = args.directory.resolve(); base.mkdir(parents=True, exist_ok=True)
    lock = (base/'watcher.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    if (base/'state.json').exists():
        raise ValueError('watcher already recorded; inspect its state before creating another')
    message = args.message_file.read_text()
    if not args.codex_bin:
        raise ValueError('codex CLI unavailable')
    state = {'status':'arming', 'watcherPid':os.getpid(), 'jobPid':args.pid,
        'thread':str(args.thread), 'startedUtc':datetime.now(timezone.utc).isoformat(),
        'sourceSha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'messageSha256':hashlib.sha256(message.encode()).hexdigest(), 'notifications':{}}
    def notify(event):
        # Record the attempt before contacting the daemon. An ambiguous result
        # is never automatically retried and cannot produce duplicate wake-ups.
        state['notifications'][event] = {'status':'sending'}; save(base/'state.json',state)
        prompt = ('User-authorized local script watcher: '+event+'.\n'
            'Watcher record: '+str(base/'state.json')+'\n'+message)
        try:
            result = subprocess.run([args.codex_bin,'queue','--thread',str(args.thread),'--message',prompt],
                capture_output=True, text=True, timeout=60)
            (base/(event+'.stdout')).write_text(result.stdout)
            (base/(event+'.stderr')).write_text(result.stderr)
            state['notifications'][event] = {'status':'queued' if result.returncode==0 else 'failed',
                'exitCode':result.returncode, 'finishedUtc':datetime.now(timezone.utc).isoformat()}
        except subprocess.TimeoutExpired:
            state['notifications'][event] = {'status':'delivery_unknown_timeout'}
        save(base/'state.json',state)

    try:
        descriptor = os.pidfd_open(args.pid)
    except ProcessLookupError:
        state['status']='job_already_exited'; save(base/'state.json',state)
        notify('completion'); return
    state['status']='waiting'; save(base/'state.json',state)
    poller = select.poll(); poller.register(descriptor,select.POLLIN)
    try:
        while True:
            # This is an OS wait, not a model/assistant polling loop. Optional
            # error-file checks are plain local Python and emit no progress.
            ended = bool(poller.poll(15000 if args.error_log else None))
            if ended:
                state['status']='job_exited'; save(base/'state.json',state)
                notify('completion'); break
            if args.error_log and args.error_log.exists() and args.error_log.stat().st_size:
                if 'issue' not in state['notifications']:
                    state['errorLog']=str(args.error_log.resolve())
                    notify('issue')
    finally:
        os.close(descriptor)
    state['status']='finished'; save(base/'state.json',state)


if __name__ == '__main__':
    main()
