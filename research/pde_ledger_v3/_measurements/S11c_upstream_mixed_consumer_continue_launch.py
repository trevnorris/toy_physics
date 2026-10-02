#!/usr/bin/env python3
"""Launch the authorized selected source/consumer diagnostic once, after its completion hook is armed."""
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
PREFIX = 'S11c_upstream_mixed_consumer_continue'
GATE = M / (PREFIX + '_gate.json')
MANIFEST = M / (PREFIX + '_inputs.json')
RUN = Path('/var/projects/toy_physics/_scratch/s11c/s11c-mixed-consumer-20261001/continuation-01')
HOOK = ROOT / 'scripts/codex_job_watch.py'
THREAD = '01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE = "The renewed user-authorized saved-evidence source/consumer continuation has completed or reported an error.\nRoot: /var/projects/toy_physics/_scratch/s11c/s11c-mixed-consumer-20261001/continuation-01.\nInspect actual coordinator/launch/stage-start, shared guard invocation/duration/enforced limits/resources/outcomes/logs, supervisor upstream_mixed_consumer_continue invocation/active, strict stderr, stdout/checks identity, failure/incomplete operation, all saved argument joins/control inputs/returns and raw/cancelled reconstruction residuals, prior-copy index, source snapshots and all posthashes. Exit status is not acceptance. JSON/source/hash inspection only outside containment; no scientific restoration or extra validator. Preserve any failure; no automatic retry.\nUser renewed work after failed result747fd8cc. Original130files/2914095bytes and all previous artifacts stay pinned unchanged. Reuse published partial observations, not completed top-level returns (old index empty). Restore selected actual operands and symbol assumptions; no source/grade/census/restriction/kernel/trace/integral replay. Only unsaved small callable/binding/formal-label context and unpublished control calculations are reconstructed. The original unfinished control body is source-extracted with identical mathematical assignments; per-record persistence added. The reconstruction helper saves jets/Dummy map/coefficients/expanded raw residual and exact cancel(together(raw)) before requiring literal zero and independent-jet-free coefficients. If raw is nonzero-looking and cancellation is exactly zero, record that actual representation finding; never infer it before inspecting. A remaining nonzero/unknown is fatal. No numeric tolerance or equation change.\nLiteral source reviews remain Claude NEEDS REVISION/Grok CLEAR FOR THIS SELECTED SOURCE/CONSUMER INSTRUMENT at3f1193e7, earlier pair0b30905b; no fresh independent CLEAR or external export. Local tests are tooling/source-equivalence only. Inspect actual scalar/divergence omissions through THETA/E_W, native pressure-slot ablation, form routing and unit joins. Scope only direct kernel(1,1)source(0,0)consumer(0,0),upper face,strict-rest-bulk LAB_HELD/RHO4_CONSTANT,omega3,cs10,saved tangents,off-shell kin0/kout0.1. Controls show sensitivity,not independent correctness. No whole-operator decoupling,loss,lower-face,calibration,drain,production repair or defect sweep.\nShared guard around existing normalization supervisor,4GiB native/cgroup in16GiB pool,zero swap,one CPU/thread,32tasks,4GiBhost reserve,desktop priority,no deadlines,RuntimeMaxUSec=infinity,Restart=no. Hook armed first to session01a0e01b-ef84-7192-817f-584cda5d339b. No scheduler change,fallback,model polling or recurring task. Preserve Lean/S11_lean,shared guard,protected builder suffixf01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2 and every old result/failure/review/incident. Canonical records outside scratch; scratch never committed. Stop with this diagnostic's evidence before production/defect work.\n"


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
    assert gate['status'] == 'READY_FOR_SELECTED_MIXED_CONSUMER_CONTINUATION'
    assert gate['independentBuildClearance'] is False and gate['pooledExecution'] is True
    assert gate['launcherSha256'] == sha(__file__)
    assert gate['workerSha256'] == sha(M / (PREFIX + '.py'))
    assert gate['manifestSha256'] == sha(MANIFEST)
    assert gate['completionMessageSha256'] == hashlib.sha256(MESSAGE.encode()).hexdigest()
    assert gate['outputDirectory'] == str(RUN / 'complete')
    assert gate['sourcePins'] == manifest['sourcePins']
    for path, digest in gate['sourcePins'].items():
        assert sha(path) == digest, path
    review = read(M / ('S11c_upstream_mixed_consumer_review_r2_record.json'))
    repair = read(M / (PREFIX + '_repair_record.json'))
    assert all(review['checks'].values()) and review['scientificJobLaunched'] is False
    assert review['reviewers']['claude']['literalVerdict'] == 'NEEDS REVISION'
    assert review['reviewers']['grok']['literalVerdict'] == 'CLEAR FOR THIS SELECTED SOURCE/CONSUMER INSTRUMENT'
    assert repair['status'] == 'LOCAL_CONTINUATION_TOOLING_TESTED_NO_FRESH_INDEPENDENT_CLEAR'
    assert repair['workerSha256'] == gate['workerSha256'] and repair['tests']['exitCode'] == 0
    assert repair['independentBuildClearance'] is False
    assert gate['reviewRecordSha256'] == sha(M / ('S11c_upstream_mixed_consumer_review_r2_record.json'))
    assert gate['repairRecordSha256'] == sha(M / (PREFIX + '_repair_record.json'))
    assert gate['authoritySha256'] == sha(M / (PREFIX + '_execution_authority.json'))
    assert read(M / (PREFIX + '_execution_authority.json'))['scienceExecutionsAuthorized'] == 1
    assert gate['reviewedWorkerSha256'] == review['sourceSnapshotWorkerSha256']
    assert manifest['resources']['memoryBytes'] == 4*1024**3
    assert manifest['resources']['durationLimits'] is None
    assert gate['command'] == [sys.executable, str(ROOT / 'scripts/s11c_guarded_run.py'),
        '--pool', 's11c-near-unity', '--memory-gib', '4', '--tasks-max', '32',
        '--log-directory', str(RUN / 'resource-guard'), '--', sys.executable,
        str(M / 'S11c_d_end_normalization_run.py'), '--parallel-prerequisite-read',
        '--run-root', str(RUN), '--stage', 'upstream_mixed_consumer_continue', '--', sys.executable,
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
