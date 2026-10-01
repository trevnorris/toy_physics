#!/usr/bin/env python3
"""Launch the approved uniform check once, after its completion hook is armed."""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import select
import shutil
import subprocess
import sys
import time
import traceback

ROOT = Path('/var/projects/toy_physics')
M = ROOT / 'research/pde_ledger_v3/_measurements'
PREFIX = 'S11c_d_near_unity_uniform'
GATE = M / (PREFIX + '_gate.json')
MANIFEST = M / (PREFIX + '_inputs.json')
RUN = ROOT / '_scratch/s11c/s11c-parallel-near-unity-20261001/uniform-01'
HOOK = ROOT / 'scripts/codex_job_watch.py'
THREAD = '01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE = '''The user-authorized selected uniform near-unity check has completed or reported an error.
Root: /var/projects/toy_physics/_scratch/s11c/s11c-parallel-near-unity-20261001/uniform-01.
Inspect coordinator/launch, actual guard invocation/duration/effective limits/resources/outcomes/logs,
supervisor near_unity_uniform invocation/active, strict scientific stderr, stdout/checks identity,
all checks/failure/incomplete operations, complete operands/returns and source/input posthashes.
Process exit is not scientific acceptance. Preserve every incomplete or failed result; no automatic retry.
Both fresh method reviewers literally CLEAR FOR THIS SELECTED UNIFORM METHOD; the frozen28-file
packet bf55b921706cc1d2c3df225fe87ec8a20405df782de771f6e50f8d60506f3752 and raw reports remain
in sibling uniform-method-review-r2. Source review is not worker/result acceptance. Read canonical
review_r2_record, implementation readiness, exact manifest/gate and launch hashes.
Scope is strict-rest-bulk LAB_HELD/RHO4_CONSTANT, omega3, saved tangents, positive effective bulk-speed
bindings across exact selected modal matching speeds. Inspect actual original-source lift/strong-pencil,
two-leg physical-depth/native radical, exact restriction/group direction/current, fresh native face and
chemical/memory/port joins, grazing limits, units and responsive controls. Unknown remains unresolved.
No defect integration/sweep, primitive calibration, drain response, complete mode census or loss claim.
Upstream scoped evidence and any missing constant-end applicability remain explicit dependencies.
Run uses no wall/native/CPU/inactivity deadline,8GiB native/cgroup in the16GiB pooled budget,zero swap,
one assigned CPU,32tasks,one thread,4GiB host reserve,desktop priority and unchanged shared guard around
existing normalization supervisor. Verify actual enforcement, RuntimeMaxUSec=infinity and Restart=no.
Scientific restoration only inside containment; JSON/source/hash/opaque-byte completion inspection is
lightweight. No new scientific validator or external packet is authorized by completion alone.
Stop after uniform evidence and applicable upstream findings, before a defect sweep. Keep all history,
Lean/S11_lean and protected builder suffix f01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2
unchanged. Canonical records outside scratch; scratch never committed. No model polling or recurring task.
'''


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n'); stream.flush(); os.fsync(stream.fileno())


def verify():
    gate, manifest = read(GATE), read(MANIFEST)
    assert gate['status'] == 'READY_FOR_SELECTED_NEAR_UNITY_UNIFORM'
    assert gate['methodClearance'] is True and gate['pooledExecution'] is True
    assert gate['launcherSha256'] == sha(__file__)
    assert gate['manifestSha256'] == sha(MANIFEST)
    assert gate['completionMessageSha256'] == hashlib.sha256(MESSAGE.encode()).hexdigest()
    assert gate['outputDirectory'] == str(RUN / 'complete')
    assert gate['sourcePins'] == manifest['sourcePins']
    for path, digest in gate['sourcePins'].items():
        assert sha(path) == digest, path
    for item in [*manifest['packets'].values(), manifest['physicalInput'],
                 manifest['plan'], manifest['reviewRecord']]:
        path = Path(item['path'])
        assert path.stat().st_size == item['bytes'] and sha(path) == item['sha256'], str(path)
    assert gate['command'] == [sys.executable, str(ROOT / 'scripts/s11c_guarded_run.py'),
        '--pool', 's11c-near-unity', '--memory-gib', '8', '--tasks-max', '32',
        '--log-directory', str(RUN / 'resource-guard'), '--', sys.executable,
        str(M / 'S11c_d_end_normalization_run.py'), '--parallel-prerequisite-read',
        '--run-root', str(RUN), '--stage', 'near_unity_uniform', '--', sys.executable,
        '-u', str(M / (PREFIX + '.py')), '--out', str(RUN / 'complete'),
        '--inputs', str(MANIFEST), '--gate', str(GATE)]
    return gate, manifest


def coordinate(descriptor):
    try:
        # Startup handshake only; this is not a computation deadline.
        ready, _, _ = select.select([descriptor], [], [], 30)
        if not ready or os.read(descriptor, 1) != b'1':
            raise RuntimeError('Completion hook not armed; refusing scientific launch')
        os.close(descriptor)
        gate, _ = verify()
        start = time.monotonic()
        save(RUN / 'stage-start.json', {'utc': datetime.now(timezone.utc).isoformat(),
             'gateSha256': sha(GATE), 'command': gate['command']})
        with (RUN / 'guard-launch.stdout').open('x') as out, (RUN / 'guard-launch.stderr').open('x') as err:
            result = subprocess.run(gate['command'], cwd=ROOT, stdin=subprocess.DEVNULL,
                                    stdout=out, stderr=err)
        save(RUN / 'coordinator-outcome.json', {'exitCode': result.returncode,
             'wallSeconds': time.monotonic() - start, 'scientificAcceptance': 'PENDING_ACTUAL_INSPECTION'})
        return result.returncode
    except BaseException:
        save(RUN / 'coordinator-failure.json', {'traceback': traceback.format_exc(), 'automaticRetry': False})
        return 1


def launch():
    gate, _ = verify()
    assert shutil.which('codex'), 'completion queue unavailable'
    RUN.mkdir(exist_ok=False)
    copies = {}
    for source, digest in {**gate['sourcePins'], str(GATE): sha(GATE), str(MANIFEST): sha(MANIFEST)}.items():
        path = Path(source)
        target = RUN / 'source' / path.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        assert sha(target) == digest
        copies[source] = {'path': str(target), 'sha256': digest, 'bytes': target.stat().st_size}
    save(RUN / 'source-snapshots.json', copies)
    (RUN / 'completion-message.txt').write_text(MESSAGE)
    reader, writer = os.pipe()
    with (RUN / 'coordinator.stdout').open('x') as out, (RUN / 'coordinator.stderr').open('x') as err:
        job = subprocess.Popen([sys.executable, str(Path(__file__).resolve()), '--coordinator', str(reader)],
              cwd=ROOT, stdin=subprocess.DEVNULL, stdout=out, stderr=err, pass_fds=(reader,), start_new_session=True)
    os.close(reader)
    try:
        with (RUN / 'watcher.stdout').open('x') as out, (RUN / 'watcher.stderr').open('x') as err:
            watcher = subprocess.Popen([sys.executable, str(HOOK), '--pid', str(job.pid),
                '--thread', THREAD, '--directory', str(RUN / 'completion-watcher'),
                '--message-file', str(RUN / 'completion-message.txt')], cwd=ROOT,
                stdin=subprocess.DEVNULL, stdout=out, stderr=err, start_new_session=True)
        save(RUN / 'launch.json', {'jobPid': job.pid, 'watcherPid': watcher.pid,
             'gateSha256': sha(GATE), 'command': gate['command'], 'thread': THREAD, 'automaticRetry': False})
        for _ in range(50):
            path = RUN / 'completion-watcher/state.json'
            if path.exists() and read(path)['status'] == 'waiting':
                os.write(writer, b'1')
                print(json.dumps({'runDirectory': str(RUN), 'jobPid': job.pid,
                      'watcherPid': watcher.pid, 'hookStatus': 'waiting'}))
                return
            if job.poll() is not None or watcher.poll() is not None:
                raise RuntimeError('Coordinator/watcher exited before arming')
            time.sleep(.1)
        raise RuntimeError('Hook did not arm')
    finally:
        os.close(writer)


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--coordinator':
        sys.exit(coordinate(int(sys.argv[2])))
    if len(sys.argv) != 1:
        raise ValueError('Unexpected launch arguments')
    launch()
