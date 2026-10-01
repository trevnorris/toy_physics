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
PREFIX = 'S11c_upstream_mixed_tanh'
GATE = M / (PREFIX + '_gate.json')
MANIFEST = M / (PREFIX + '_inputs.json')
RUN = Path('/var/projects/toy_physics/_scratch/s11c/s11c-mixed-tanh-20261001/diagnostic-01')
HOOK = ROOT / 'scripts/codex_job_watch.py'
THREAD = '01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE = 'The user-authorized focused mixed-tanh diagnostic has completed or reported an error.\nRoot: /var/projects/toy_physics/_scratch/s11c/s11c-mixed-tanh-20261001/diagnostic-01.\nInspect coordinator/launch/stage-start, shared guard duration/enforced limits/resources/outcomes/logs, supervisor upstream_mixed_tanh invocation/active, strict scientific stderr, stdout/checks identity, failure/incomplete operation, all exact input/return/identity/control JSON, source snapshots and every posthash. Process exit is not physics acceptance. No scientific restoration outside containment; completion inspection uses JSON/source/hash/opaque bytes only.\nThis is the one authorized upper-face bare-impedance eta*sigma_W diagnostic on actual tanh and common physical outgoing dispersion, omega3 cs10 saved tangents, profile kin0 kout1/10. No light excitation, closed slab response, defect sweep, integral or loss computation. Inspect actual native source/flat/linear/normal/profile/jet joins, original boundary residuals and coefficient, delta/PV transfer cancellation, both transfer assignments, sign/phase ingredients, branch endpoints/tail, units and responsive native slope/sheet controls. Unknown remains unresolved. Do not infer an incomplete closed operator solely from a bare direct slot.\nBoth reports finished before local repairs. Literal Claude NEEDS REVISION and Grok CLEAR FOR THIS FOCUSED MIXED-GRADE INSTRUMENT remain preserved at6d645f57 and sibling build-review, packetb7b6b8dd5648e95714703d41a0db7cf0f7c6b07caad9dcfece6f137f8efec361. Both supported the selected mathematics. Local N1 persistence and N4 routing fixes save operands before checks and use already-computed native factors with exact equality guards. N2 simplifier failure was not established; N3 duplicate fragments was refuted by actual AST. No new independent CLEAR or reviewer rerun; gate records independentBuildClearance=false and standing user tooling-fix authority. Eight stdlib tests pass, not a scientific result. Preserve all reviewed bytes and actual repair record.\nNo deadline of any kind. Shared guard around existing normalization supervisor,4GiB native/cgroup within16GiB pool,zero swap,one assignedCPU/thread,32tasks,4GiB host reserve,desktop-managed priority,RuntimeMaxUSec=infinity,Restart=no. No scheduler change, fallback or automatic retry. Preserve every failure and partial return. STOP after diagnostic evidence and disposition; no additional validator, producer regeneration, closed-operator repair or defect sweep is authorized by this completion. Keep selected uniform results, all history/incident/review debt,Lean/S11_lean and protected builder suffixf01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2. Canonical records outside scratch; scratch never committed. Hook targets THIS session01a0e01b-ef84-7192-817f-584cda5d339b. No model polling or recurring task.\n'


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
    assert gate['status'] == 'READY_FOR_FOCUSED_MIXED_TANH_DIAGNOSTIC'
    assert gate['independentBuildClearance'] is False and gate['pooledExecution'] is True
    assert gate['launcherSha256'] == sha(__file__)
    assert gate['workerSha256'] == sha(M / (PREFIX + '.py'))
    assert gate['manifestSha256'] == sha(MANIFEST)
    assert gate['completionMessageSha256'] == hashlib.sha256(MESSAGE.encode()).hexdigest()
    assert gate['outputDirectory'] == str(RUN / 'complete')
    assert gate['sourcePins'] == manifest['sourcePins']
    for path, digest in gate['sourcePins'].items():
        assert sha(path) == digest, path
    review = read(M / (PREFIX + '_review_record.json'))
    repair = read(M / (PREFIX + '_repair_record.json'))
    assert review['allMetadataChecksPassed'] is True and review['noScienceRun'] is True
    assert review['reports']['claude']['literalVerdict'] == 'NEEDS REVISION'
    assert review['reports']['grok']['literalVerdict'] == 'CLEAR FOR THIS FOCUSED MIXED-GRADE INSTRUMENT'
    assert repair['status'] == 'LOCAL_TOOLING_REPAIRS_TESTED_NO_FRESH_INDEPENDENT_CLEAR'
    assert repair['workerSha256'] == gate['workerSha256'] and repair['tests']['exitCode'] == 0
    assert repair['independentBuildClearance'] is False
    assert gate['reviewRecordSha256'] == sha(M / (PREFIX + '_review_record.json'))
    assert gate['repairRecordSha256'] == sha(M / (PREFIX + '_repair_record.json'))
    assert manifest['resources']['memoryBytes'] == 4*1024**3
    assert manifest['resources']['durationLimits'] is None
    assert gate['command'] == [sys.executable, str(ROOT / 'scripts/s11c_guarded_run.py'),
        '--pool', 's11c-near-unity', '--memory-gib', '4', '--tasks-max', '32',
        '--log-directory', str(RUN / 'resource-guard'), '--', sys.executable,
        str(M / 'S11c_d_end_normalization_run.py'), '--parallel-prerequisite-read',
        '--run-root', str(RUN), '--stage', 'upstream_mixed_tanh', '--', sys.executable,
        '-u', str(M / (PREFIX + '.py')), '--out', str(RUN / 'complete'),
        '--inputs', str(MANIFEST)]
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
