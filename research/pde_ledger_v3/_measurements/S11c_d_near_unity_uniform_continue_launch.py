#!/usr/bin/env python3
"""Launch the approved saved-return uniform continuation once, after its completion hook is armed."""
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
PREFIX = 'S11c_d_near_unity_uniform_continue'
GATE = M / (PREFIX + '_gate.json')
MANIFEST = M / (PREFIX + '_inputs.json')
RUN = ROOT / '_scratch/s11c/s11c-parallel-near-unity-20261001/uniform-continuation-01'
HOOK = ROOT / 'scripts/codex_job_watch.py'
THREAD = '01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE = '''The user-authorized saved-return near-unity denominator continuation has completed or reported an error.
Root: /var/projects/toy_physics/_scratch/s11c/s11c-parallel-near-unity-20261001/uniform-continuation-01.
Inspect coordinator/launch/stage-start, actual shared-guard invocation/duration/enforced limits/resources/outcomes/logs,
supervisor near_unity_uniform_continue invocation/active, strict scientific stderr, stdout/checks identity,
checks/failure/incomplete operations, journal operands/returns, prior argument joins and every posthash.
Process exit is not acceptance. Preserve all failures; no automatic retry or validator.
User said: k. let's do that next step. This approves one exact saved-denominator certification and continuation
of the unfinished points/controls. Original uniform-01 remains UNRESOLVED:27 completed returns,44 domain refusals,
32 files/46689296 bytes,5383 blobs; own exact LEFT/RIGHT matches and four two-sided limits supported.
Verify all27 complete returns restored without functions, byte-identical32-file prior copy, actual source/input
argument joins, completed point reuse and preserved original refusals. No original source/native/limit replay.
Only exact symbol/EndBinding attribute context and unsaved small physical-point substitution are reconstituted.
The old domain values are reused without rebinding. Read every certificate: original flags or exact real/imaginary
components, literal zero reconstruction and signed nonzero component. No numerical smallness or unknown-to-true
promotion. Certificate self-checks include zero/nonfinite/unbound rejection. All other point/current/face/sheet
formulas and tolerances stay in the pinned original worker; the nongrazing native sheet controls must actually run.
Both original method reports literally CLEAR FOR THIS SELECTED UNIFORM METHOD, packet
bf55b921706cc1d2c3df225fe87ec8a20405df782de771f6e50f8d60506f3752; no new independent worker/result review.
The repair is an exact domain-decision implementation, not new equations, a relaxed tolerance or new method.
Strict-rest-bulk LAB_HELD/RHO4_CONSTANT,omega3,saved tangents,effective modal-speed sweep only; nonuniform/direct
mixed-grade applicability stays conditional. No defect sweep, loss, primitive calibration or draining inference.
No wall/native/CPU/inactivity deadline.8GiB native/cgroup in16GiB pooled budget,zero swap,one assigned CPU,
32tasks,one thread,4GiB host reserve,desktop priority,unchanged shared guard around normalization supervisor.
Verify RuntimeMaxUSec=infinity,Restart=no and actual resources. No scheduler changes or unguarded fallback.
Scientific payload restoration only within containment; completion inspection JSON/source/hash/opaque bytes only.
Stop after uniform evidence and applicable upstream findings. Preserve every prior result/failure/review/incident,
Lean/S11_lean,shared guard and protected suffixf01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2.
Canonical records outside scratch; scratch never committed. Hook targets01a0e01b-ef84-7192-817f-584cda5d339b.
No model polling or recurring tasks.
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
    assert gate['status'] == 'READY_FOR_SAVED_UNIFORM_DOMAIN_CONTINUATION'
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
        '--run-root', str(RUN), '--stage', 'near_unity_uniform_continue', '--', sys.executable,
        '-u', str(M / (PREFIX + '.py')), '--out', str(RUN / 'complete'),
        '--inputs', str(MANIFEST), '--gate', str(GATE)]
    assert gate['authorityRecord'] == manifest['authorityRecord']
    for item in (manifest['authorityRecord'], manifest['priorManifest']):
        assert sha(item['path']) == item['sha256']
    for relative, item in manifest['priorComplete']['files'].items():
        path = Path(manifest['priorComplete']['path']) / relative
        assert path.stat().st_size == item['bytes'] and sha(path) == item['sha256']
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
