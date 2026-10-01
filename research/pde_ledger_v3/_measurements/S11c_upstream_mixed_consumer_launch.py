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
PREFIX = 'S11c_upstream_mixed_consumer'
GATE = M / (PREFIX + '_gate.json')
MANIFEST = M / (PREFIX + '_inputs.json')
RUN = Path('/var/projects/toy_physics/_scratch/s11c/s11c-mixed-consumer-20261001/diagnostic-01')
HOOK = ROOT / 'scripts/codex_job_watch.py'
THREAD = '01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE = 'The user-authorized selected direct mixed-term source/consumer diagnostic has completed or reported an error.\nRoot: /var/projects/toy_physics/_scratch/s11c/s11c-mixed-consumer-20261001/diagnostic-01.\nInspect coordinator/launch/stage-start, shared guard invocation/duration/enforced resources/outcomes/logs, supervisor upstream_mixed_consumer invocation/active, strict scientific stderr, stdout/checks byte identity, failure/incomplete operation, all journal input/return and native source/census/grade/restriction/consumer/control/unit evidence, source snapshots and every posthash. Exit status is not scientific acceptance. Completion inspection is JSON/source/hash/opaque bytes only; no scientific restoration outside containment. Preserve failure and STOP; no automatic retry or new job.\nBoth source-only reviews finished before local repairs. Literal Claude NEEDS REVISION/Grok CLEAR FOR THIS SELECTED SOURCE/CONSUMER INSTRUMENT remain preserved at3f1193e7 and sibling build-review-r2; packet3c10f7462d48c5ee0b56f77e754e0c3878c5d1e892391a216c6aff2bf5ade731. Earlier pair at0b30905b remains. No fresh independent CLEAR or new review. Local repairs enforce existing profile/zero-grade contracts, evaluate the same exact constant nonzero predicate via real/imaginary components only when needed, and save all control/unit evidence before final response guards. Sixteen stdlib tests pass; scientific operands unchanged (only nonzero-state assignment differs among153 prior assignments). No higher-grade calculation added; gate records independentBuildClearance=false. Reviewer CLI stderr is not scientific stderr.\nScope is only upper-face direct kernel(1,1) times source(0,0) times consumer(0,0), strict-rest-bulk LAB_HELD/RHO4_CONSTANT omega3 cs10, saved tangents and off-shell profile kin0/kout0.1. Inspect complete native chemical/density/velocity and pressure/jet census, exact source joins and normalization, regularity/extracted grade, arbitrary curl/scalar/longitudinal restriction, actual THETA/E_W contraction, epsilon/units, scalar/divergence source omissions through actual rows and pressure-slot ablation. Controls establish sensitivity, not independent correctness; reverse U absence is census-based and wrong routing form-only. Unknown is unresolved. Whole inherited convolution remains symbolic and enters once, never another middle integral. No lower-face correction, full mode, current, loss, primitive calibration, drain response, production repair or defect sweep.\nNo deadlines of any kind. Shared guard around existing normalization supervisor,4GiB native/cgroup within16GiB pooled budget,zero swap,one assignedCPU/thread,32tasks,4GiBhost reserve,desktop-managed priority,RuntimeMaxUSec=infinity,Restart=no. No scheduler change, unguarded fallback or retry. Stop with diagnostic evidence and applicable disposition. Preserve every prior accepted/failure/incident/review artifact,Lean/S11_lean,shared guard and protected builder suffixf01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2. Canonical records outside scratch;scratch never committed. Hook targets session01a0e01b-ef84-7192-817f-584cda5d339b;no model polling or recurring task.\n'


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
    assert gate['status'] == 'READY_FOR_SELECTED_MIXED_CONSUMER_DIAGNOSTIC'
    assert gate['independentBuildClearance'] is False and gate['pooledExecution'] is True
    assert gate['launcherSha256'] == sha(__file__)
    assert gate['workerSha256'] == sha(M / (PREFIX + '.py'))
    assert gate['manifestSha256'] == sha(MANIFEST)
    assert gate['completionMessageSha256'] == hashlib.sha256(MESSAGE.encode()).hexdigest()
    assert gate['outputDirectory'] == str(RUN / 'complete')
    assert gate['sourcePins'] == manifest['sourcePins']
    for path, digest in gate['sourcePins'].items():
        assert sha(path) == digest, path
    review = read(M / (PREFIX + '_review_r2_record.json'))
    repair = read(M / (PREFIX + '_local_repair_r2_record.json'))
    assert all(review['checks'].values()) and review['scientificJobLaunched'] is False
    assert review['reviewers']['claude']['literalVerdict'] == 'NEEDS REVISION'
    assert review['reviewers']['grok']['literalVerdict'] == 'CLEAR FOR THIS SELECTED SOURCE/CONSUMER INSTRUMENT'
    assert repair['status'] == 'LOCAL_TOOLING_REPAIRS_TESTED_NO_FRESH_INDEPENDENT_CLEAR'
    assert repair['workerSha256'] == gate['workerSha256'] and repair['tests']['exitCode'] == 0
    assert repair['independentBuildClearance'] is False
    assert gate['reviewRecordSha256'] == sha(M / (PREFIX + '_review_r2_record.json'))
    assert gate['repairRecordSha256'] == sha(M / (PREFIX + '_local_repair_r2_record.json'))
    assert gate['authoritySha256'] == sha(M / (PREFIX + '_execution_authority.json'))
    assert read(M / (PREFIX + '_execution_authority.json'))['scienceExecutionsAuthorized'] == 1
    assert gate['reviewedWorkerSha256'] == review['sourceSnapshotWorkerSha256']
    assert manifest['resources']['memoryBytes'] == 4*1024**3
    assert manifest['resources']['durationLimits'] is None
    assert gate['command'] == [sys.executable, str(ROOT / 'scripts/s11c_guarded_run.py'),
        '--pool', 's11c-near-unity', '--memory-gib', '4', '--tasks-max', '32',
        '--log-directory', str(RUN / 'resource-guard'), '--', sys.executable,
        str(M / 'S11c_d_end_normalization_run.py'), '--parallel-prerequisite-read',
        '--run-root', str(RUN), '--stage', 'upstream_mixed_consumer', '--', sys.executable,
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
