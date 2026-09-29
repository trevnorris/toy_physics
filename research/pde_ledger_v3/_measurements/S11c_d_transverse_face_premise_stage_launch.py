#!/usr/bin/env python3
"""Launch the exact cleared premise worker once, with its completion hook first."""
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
M = ROOT/'research/pde_ledger_v3/_measurements'
BASE = ROOT/'_scratch/s11c/s11c-d-transverse-face-20260928'
RUN = BASE/'production'
PREFIX = 'S11c_d_transverse_face_premise'
WORKER = M/(PREFIX+'.py')
INPUT = M/(PREFIX+'_inputs.json')
GATE = M/(PREFIX+'_gate.json')
APPROVAL = M/(PREFIX+'_authorization.json')
REVIEW = M/(PREFIX+'_review_r5_record.json')
PREPARATION = M/(PREFIX+'_review_r5_preparation.json')
GUARD = ROOT/'scripts/s11c_guarded_run.py'
SUPERVISOR = M/'S11c_d_end_normalization_run.py'
HOOK = ROOT/'scripts/codex_job_watch.py'
THREAD = '01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE = '''The explicitly approved job 1 transverse/face premise instrument has completed or reported an error.
Root: /var/projects/toy_physics/_scratch/s11c/s11c-d-transverse-face-20260928/production.
Inspect coordinator/launch, guard invocation/effective limits/limit validation/resource samples/child outcome/outcome and logs; supervisor transverse_face_premise invocation/active; strict worker stderr; stdout/checks identity; checks/failure/incomplete operation; complete operation operands/returns/index, artifact index, all posthashes and actual saved premise/face/current/coverage/control evidence. Exit code is not scientific acceptance. Preserve all failures and partial output, no retry.

User explicitly said "Approve jobs 1-2 only" after the source-only frequency comparison committed at d2f12b76. Job 1 is the exact fifth-round independently cleared worker, SHA256 ebf2739a2f8dca9ed4130b89606a36a27be29ef5d5c0915756b7118fcbd6fae5. Both literal reviews CLEAR FOR THIS BOUNDED TRANSVERSE/FACE BUILD; this is build clearance, not premise or loss acceptance. Job 2 is ONE separately guarded, independently reviewed saved-output reader, prepared against these actual complete or partial returns and selected old end/forcing/contrast payloads needed for the omega=1/omega=3 comparison. No constructor/root/mode/transform/integral replay or repair disguised as validation. Do not ask repetitive science-stage approval already covered; exact external packet consent requirements still apply. Review job 2 before running it. Never restore scientific payloads during lightweight JSON/hash inspection.

Both jobs have 900s outer/840s native, 2GiB/zero swap/one CPU/nice15/32tasks/one native thread via unchanged shared scripts/s11c_guarded_run.py around S11c_d_end_normalization_run.py. No unlimited inheritance, overlap, fallback or automatic retry. One worker, hook first to this session 01a0e01b-ef84-7192-817f-584cda5d339b. This is job 1 of the four-job ceiling; only jobs 1-2 are authorized. After job 2 return actual results/costs/coverage and omega=1 versus omega=3 comparison, including which consumed candidates need clearance. STOP before the one-day method-route decision, five-day authoring phase or jobs 3-4. No first-order outgoing-field/power construction or loss magnitude is supplied by these two jobs.

Keep availability/absence/unresolved distinctions, partial coverage and limited reference face evidence. c2 is velocity-only, orientation-blind/index-literal. Zero reference drive does not establish all-grade matched lossless ends or zero total loss. Retain physical damping. Thickness classification remains deferred; calibration of the analog-light frequency band is explicitly OPEN. No full Green/FORM/A11/A12, radiation fraction, calibrated light claim or historical review-debt erasure.

Preserve all accepted science/failures/incident/review literals, Lean/S11_lean/shared guard and protected suffix f01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2. Canonical records outside scratch; runtime scratch remains ignored, never committed. Explicit diff exclusions for generated tracked data above 1MiB. No model polling or recurring tasks.
'''


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()


def route(path):
    path = Path(path)
    return dict(path=str(path), canonicalPath=str(path.resolve(strict=True)),
                bytes=path.stat().st_size, sha256=sha(path))


def read(path):
    return json.loads(path.read_text())


def save(path, value):
    with path.open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n'); f.flush(); os.fsync(f.fileno())


def preflight():
    gate = read(GATE)
    assert gate['status'] == 'READY_FOR_ONE_GUARDED_TRANSVERSE_FACE_PREMISE_JOB'
    assert gate['scienceAuthorizationAfterCostedPlanGo'] is True
    assert gate['userScienceApprovalRecord'] == route(APPROVAL)
    approval = read(APPROVAL)
    assert approval['userReply'] == 'Approve jobs 1-2 only'
    assert approval['authorizedScienceJobOrdinals'] == [1, 2]
    assert approval['stopAfterJob2AndFrequencyComparison'] is True
    assert approval['workerSha256'] == gate['workerSha256'] == sha(WORKER)
    assert approval['inputManifestSha256'] == gate['inputManifestSha256'] == sha(INPUT)
    assert gate['independentBuildClearance'] is True and gate['scienceJobOrdinal'] == 1
    assert gate['reviewRecord'] == route(REVIEW)
    review = read(REVIEW)
    assert review['independentBuildClearance'] is True
    assert all(v['literalVerdict'] == 'CLEAR FOR THIS BOUNDED TRANSVERSE/FACE BUILD'
               for v in review['reviews'].values()) and len(review['reviews']) == 2
    assert gate['launcherSha256'] == sha(Path(__file__))
    assert gate['completionMessageSha256'] == hashlib.sha256(MESSAGE.encode()).hexdigest()
    assert gate['guardSha256'] == sha(GUARD) and gate['supervisorSha256'] == sha(SUPERVISOR)
    assert gate['completionHookSha256'] == sha(HOOK)
    assert (gate['seconds'], gate['nativeSeconds'], gate['memoryGiB'], gate['tasksMax']) == (900, 840, 2, 32)
    assert (gate['swapMax'], gate['cpuCount'], gate['nice'], gate['nativeThreads']) == (0, 1, 15, 1)
    assert gate['automaticRetry'] is False
    for item in read(PREPARATION)['sourceRecords']:
        assert sha(ROOT/item['path']) == item['sourceSha256'], item['path']
    spec = read(INPUT)
    for name, record in spec['inputs'].items():
        assert route(record['path']) == record, name
    for join in spec['checkpointJoins']:
        assert read(Path(join['checkpoint']))['artifacts'][join['artifact']]['sha256'] == spec['inputs'][join['inputKey']]['sha256']
    assert read(Path(spec['inputs']['physicalInputFile']['path'])) == spec['physicalInput']
    builder = (M/'S11c_d_sympy_builder_report.md').read_bytes()
    suffix = builder[builder.index(b'## Retained user-approved solver/export contract'):]
    assert hashlib.sha256(suffix).hexdigest() == gate['protectedBuilderSuffixSha256']
    return dict(workerSha256=sha(WORKER), inputManifestSha256=sha(INPUT),
                gateSha256=sha(GATE), launcherSha256=sha(Path(__file__)))


def command():
    return ['/usr/bin/python3', str(GUARD), '--log-directory', str(RUN/'resource-guard'),
            '--seconds', '900', '--memory-gib', '2', '--tasks-max', '32', '--',
            '/usr/bin/python3', str(SUPERVISOR), '--run-root', str(RUN),
            '--stage', 'transverse_face_premise', '--', '/usr/bin/python3', '-u', str(WORKER),
            '--input-manifest', str(INPUT), '--gate-receipt', str(GATE),
            '--run-directory', str(RUN/'complete')]


def coordinate(descriptor):
    try:
        ready, _, _ = select.select([descriptor], [], [], 30)
        if not ready or os.read(descriptor, 1) != b'1':
            raise RuntimeError('Completion hook did not arm; no workload launched')
        os.close(descriptor)
        pins = preflight()
        assert pins == read(RUN/'launch.json')['pins']
        started = time.monotonic()
        save(RUN/'stage-start.json', dict(startedUtc=datetime.now(timezone.utc).isoformat(),
             command=command(), pins=pins, automaticRetry=False))
        with (RUN/'guard-launch.stdout').open('x') as out, (RUN/'guard-launch.stderr').open('x') as err:
            result = subprocess.run(command(), cwd=ROOT, stdin=subprocess.DEVNULL, stdout=out, stderr=err)
        save(RUN/'coordinator-outcome.json', dict(exitCode=result.returncode,
             wallSeconds=time.monotonic()-started, scientificAcceptance='PENDING',
             finishedUtc=datetime.now(timezone.utc).isoformat()))
        if result.returncode:
            (RUN/'issues.log').write_text('Guarded stage failed; inspect preserved records. No retry.\n')
        return result.returncode
    except BaseException:
        failure = traceback.format_exc()
        if not (RUN/'coordinator-outcome.json').exists():
            save(RUN/'coordinator-outcome.json', dict(status='COORDINATOR_ERROR_NO_RETRY',
                 traceback=failure, scientificAcceptance='UNACCEPTED', finishedUtc=datetime.now(timezone.utc).isoformat()))
        (RUN/'issues.log').write_text(failure)
        raise


def launch():
    pins = preflight()
    assert shutil.which('codex'), 'Completion queue unavailable'
    RUN.mkdir(exist_ok=False)
    snapshot = RUN/'source'; snapshot.mkdir()
    paths = {ROOT/item['path'] for item in read(PREPARATION)['sourceRecords']}
    paths.update((GATE, APPROVAL, Path(__file__).resolve(), REVIEW, PREPARATION, HOOK,
                  M/'S11c_d_pilot_frequency_comparison.md', M/'S11c_d_outgoing_field_power_costed_plan.md'))
    for path in sorted(paths):
        target = snapshot/path.relative_to(ROOT); target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target); assert sha(target) == sha(path)
    with (RUN/'completion-message.txt').open('x') as f:
        f.write(MESSAGE); f.flush(); os.fsync(f.fileno())
    read_fd, write_fd = os.pipe()
    with (RUN/'coordinator.stdout').open('x') as out, (RUN/'coordinator.stderr').open('x') as err:
        job = subprocess.Popen(['/usr/bin/python3', str(Path(__file__)), '--coordinator', str(read_fd)],
              cwd=ROOT, stdin=subprocess.DEVNULL, stdout=out, stderr=err,
              pass_fds=(read_fd,), start_new_session=True)
    os.close(read_fd)
    try:
        with (RUN/'watcher.stdout').open('x') as out, (RUN/'watcher.stderr').open('x') as err:
            watcher = subprocess.Popen(['/usr/bin/python3', str(HOOK), '--pid', str(job.pid),
                 '--thread', THREAD, '--directory', str(RUN/'completion-watcher'),
                 '--message-file', str(RUN/'completion-message.txt'), '--error-log', str(RUN/'issues.log')],
                 cwd=ROOT, stdin=subprocess.DEVNULL, stdout=out, stderr=err, start_new_session=True)
        save(RUN/'launch.json', dict(jobPid=job.pid, watcherPid=watcher.pid, pins=pins,
             command=command(), thread=THREAD, maximumScientificWorkers=1,
             completionMessageSha256=hashlib.sha256(MESSAGE.encode()).hexdigest()))
        for _ in range(50):
            state = RUN/'completion-watcher/state.json'
            if state.exists() and read(state)['status'] == 'waiting':
                os.write(write_fd, b'1')
                print(json.dumps(dict(jobPid=job.pid, watcherPid=watcher.pid,
                      hookStatus='waiting', runDirectory=str(RUN))))
                return
            if watcher.poll() is not None or job.poll() is not None:
                raise RuntimeError('Coordinator or hook exited before arming')
            time.sleep(.1)
        raise RuntimeError('Completion hook did not arm')
    finally:
        os.close(write_fd)


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--coordinator':
        sys.exit(coordinate(int(sys.argv[2])))
    if len(sys.argv) != 1:
        raise ValueError('Unexpected launch arguments')
    launch()
