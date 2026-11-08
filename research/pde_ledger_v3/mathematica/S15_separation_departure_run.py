#!/usr/bin/env python3
"""Sequential, fresh-kernel launcher for the self-contained Wolfram engine.

Run this launcher through s11c_guarded_run.py. It never imports or transforms
mathematical expressions: the base pass emits an object-name manifest, and each
fresh kernel reconstructs one named item-4 object from the Wolfram source.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import threading


def sha256(path: Path) -> str:
    with path.open('rb') as source:
        return hashlib.file_digest(source, 'sha256').hexdigest()


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


DECISION_FAILURES = {
    f'region/component_{i}/{kind}_labels_{{drain1, drain2}}':
        'User decision 2026-11-06: complete item-4 conditions were not returned by the recorded prior attempts. '
        'This not-obtained label is set by that decision; the exponent-sensitive attempt was not made in this final run. '
        'The retained sufficient conditions are printed above. The run record lists the attempt history and the larger-memory follow-up.'
    for i in range(1,4) for kind in ('momentum','rate')
}

# User-authorized item-4 failure policies. They do not limit base construction.
ATTEMPT_SECONDS = 6 * 60 * 60


def memory_peak() -> int:
    for line in Path('/proc/self/cgroup').read_text().splitlines():
        if line.startswith('0::'):
            return int((Path('/sys/fs/cgroup') / line[3:].lstrip('/') / 'memory.peak').read_text())
    raise RuntimeError('Cgroup peak memory accounting is unavailable.')


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--transcript', type=Path, required=True)
    parser.add_argument('--work-directory', type=Path, required=True)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--exclude-object', action='append', default=[])
    parser.add_argument('--only-object', action='append', default=[])
    args = parser.parse_args()
    source = args.source.resolve(strict=True)
    # Do not resolve the transcript: it must not follow an annex link.
    transcript = args.transcript.absolute()
    work = args.work_directory.absolute()
    if transcript.is_symlink():
        raise RuntimeError('Transcript must be a plain file; remove its link itself before a fresh run.')
    if args.resume:
        if not transcript.is_file() or not (work/'execution.json').is_file():
            raise RuntimeError('Resume needs the existing plain transcript and execution record.')
    else:
        if transcript.exists():raise RuntimeError('Transcript must be absent before a fresh run.')
        work.mkdir(parents=True, exist_ok=False)
    transcript.parent.mkdir(parents=True, exist_ok=True)
    pinned_hash = sha256(source)
    command = ['math', '-script', str(source)]
    state = {'source': str(source), 'source_sha256': pinned_hash,
             'transcript': str(transcript), 'command': command,
             'started': utc(), 'kernels': []}
    state_path = work / 'execution.json'
    if args.resume:
        state=json.loads(state_path.read_text())
        if state['source_sha256']!=pinned_hash:raise RuntimeError('Resume source differs from the original construction.')
    state.setdefault('segments',[]).append({'started':utc(),'guard_run':os.environ.get('S15_WL_RUN_NAME',''), 'resume':args.resume,'exclude':args.exclude_object,'only':args.only_object})

    def save() -> None:
        temporary = work / 'execution.next.json'
        temporary.write_text(json.dumps(state, indent=2)+'\n')
        temporary.replace(state_path)

    def run_kernel(mode: str, object_name: str = '', attempt_failure: str = '') -> None:
        if sha256(source) != pinned_hash:
            raise RuntimeError('Source changed during this run.')
        environment = os.environ.copy()
        environment.update(S15_WL_MODE=mode, S15_WL_OBJECT=object_name,
                           S15_WL_TRANSCRIPT=str(transcript), S15_WL_ATTEMPT_FAILURE=attempt_failure)
        record = {'mode': mode, 'object': object_name, 'started': utc(),
                  'command': command, 'source_sha256': pinned_hash, 'attempt_failure_label': attempt_failure,
                  'guard_run':os.environ.get('S15_WL_RUN_NAME','')}
        state['kernels'].append(record)
        print('WORKING_ON: '+(object_name or 'base pass')+' / start fresh Wolfram kernel', flush=True)
        start = time.monotonic()
        # Exactly one CAS child is alive. The user authorized stopping an
        # item-4 attempt at six hours then labelling it.
        child = subprocess.Popen(command, env=environment, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT, text=True, bufsize=1)
        record['pid'] = child.pid
        save()
        activity = {'last_step': '', 'last_output': start, 'attempt_start': None,
                    'memory_exhausted': False}
        stopped_kind = None
        child_log = work / ('kernel_' + str(len(state['kernels'])) + '.stdout')
        def read_output():
            with child_log.open('w') as kernel_output:
                for line in child.stdout:
                    print(line, end='', flush=True)
                    kernel_output.write(line)
                    kernel_output.flush()
                    activity['last_output'] = time.monotonic()
                    if line.startswith('WORKING_ON:'):
                        activity['last_step'] = line.strip()
                        if 'item4 exponent-sensitive /' in line and activity['attempt_start'] is None:
                            activity['attempt_start'] = time.monotonic()
                            record['exact_attempt_started'] = utc()
                    if 'No more memory available.' in line:
                        activity['memory_exhausted'] = True
        reader = threading.Thread(target=read_output)
        reader.start()
        while child.poll() is None:
            now = time.monotonic()
            attempt_start = activity['attempt_start']
            if mode == 'item4' and not attempt_failure and attempt_start is not None:
                if now - attempt_start >= ATTEMPT_SECONDS:
                    stopped_kind = 'six_hour_stop'
                if stopped_kind:
                    record['author_stop_kind'] = stopped_kind
                    record['attempt_seconds_at_stop'] = now - attempt_start
                    record['kernel_seconds_at_stop'] = now - start
                    record['last_step_at_stop'] = activity['last_step']
                    record['guarded_run_peak_at_stop'] = memory_peak()
                    save()
                    print('STOPPING_ITEM4_ATTEMPT: '+json.dumps(record), flush=True)
                    # Stop only this failed CAS child; the guarded launcher stays
                    # alive to print its label and continue the other objects.
                    child.kill()
                    break
            time.sleep(1)
        reader.join()
        last_step = activity['last_step']
        memory_exhausted = activity['memory_exhausted']
        code = child.wait()
        record.update(last_step=last_step, memory_exhausted=memory_exhausted,
                      guarded_run_peak_bytes=memory_peak(), stdout=str(child_log))
        record.update(returncode=code, ended=utc(), wall_seconds=time.monotonic()-start)
        save()
        print('FINISHED_KERNEL: '+json.dumps(record), flush=True)
        if code:
            if mode == 'item4' and (memory_exhausted or stopped_kind) and not attempt_failure:
                kind = stopped_kind or 'memory_exhausted'
                if kind=='memory_exhausted' and os.environ.get('S15_ALONE_MAXIMUM_CAP')!='1':
                    state.setdefault('pending_alone_memory_retries',{})[object_name]=record
                    save()
                    print('ITEM4_PENDING_ALONE_MEMORY_RETRY: '+object_name,flush=True)
                    return
                cause = ('ran out of memory' if kind == 'memory_exhausted' else
                         'was stopped after 6 hours or more without returning')
                label = ('The exponent-sensitive native Euler radial, differential-matrix, '
                         'series and parameter-stratum construction '+cause+'. Run work directory: '
                         +str(work)+'. Last reported step: '+last_step
                         +'; kernel elapsed seconds: '+str(record['wall_seconds'])
                         +'; exact-attempt elapsed seconds: '+str(record.get('attempt_seconds_at_stop', time.monotonic()-(activity['attempt_start'] or start)))
                         +'; guarded run peak memory at failure: '+str(record['guarded_run_peak_bytes'])
                         +' bytes. Partial attempt: '+str(child_log)+'. No complete classification was returned.')
                state.setdefault('failed_objects', {})[object_name] = {'kind':kind,'label':label}
                save()
                # This pass only reconstructs the retained sufficient conditions
                # and prints the failed-attempt label. It does not retry the attempt.
                run_kernel('item4', object_name, label)
                return
            raise RuntimeError(f'Wolfram {mode} kernel exited {code}: {object_name}')
        marker = ('END_OF_BASE_PASS: "S15 separation departure Wolfram objects"'
                  if mode == 'base' else
                  'END_OF_ITEM4_PASS: '+json.dumps(object_name))
        with transcript.open('rb') as output:
            output.seek(max(0, transcript.stat().st_size - 65536))
            tail = output.read().decode('utf-8', errors='replace')
        if marker not in tail.splitlines():
            raise RuntimeError(f'Kernel omitted its completion marker: {mode} {object_name}')

    save()
    if not args.resume:run_kernel('base')
    manifest = []
    prefix = 'ITEM4_WORKER_MANIFEST_ENTRY: '
    with transcript.open() as output:
        for line in output:
            if line.startswith(prefix):
                name = json.loads(line[len(prefix):])
                if not isinstance(name, str) or name in manifest:
                    raise RuntimeError('Invalid or duplicate object-name manifest entry.')
                manifest.append(name)
    state['item4_manifest'] = manifest
    save()
    # Reconstruct the region objects first; every manifest entry still runs once.
    worker_order = sorted(manifest, key=lambda name: (not name.startswith('region/'), not ('region/' in name and '{drain1, drain2}' in name)))
    state['item4_worker_order'] = worker_order
    save()
    completed={k['object'] for k in state['kernels'] if k.get('mode')=='item4' and k.get('returncode')==0}
    for name in worker_order:
        if name in completed or name in args.exclude_object or (args.only_object and name not in args.only_object):continue
        failure = DECISION_FAILURES.get(name, '')
        if failure:
            state.setdefault('failed_objects', {})[name] = {'kind': 'user_decision_2026_11_06', 'label': failure}
            save()
        run_kernel('item4', name, failure)
    completed={k['object'] for k in state['kernels'] if k.get('mode')=='item4' and k.get('returncode')==0}
    remaining=[name for name in worker_order if name not in completed]
    state['remaining_objects']=remaining
    state['segments'][-1].update(ended=utc(),remaining=remaining)
    save()
    if remaining:
        print('WOLFRAM_EXECUTION_PENDING_OBJECTS: '+json.dumps(remaining),flush=True)
        return 0
    with transcript.open('a') as output:
        output.write('END_OF_TRANSCRIPT: "S15 separation departure Wolfram objects"\n')
        output.flush()
    state['ended'] = utc()
    state['transcript_sha256'] = sha256(transcript)
    save()
    print('FINISHED_WOLFRAM_EXECUTION: '+str(state_path), flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
