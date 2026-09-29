#!/usr/bin/env python3
"""Launch the exact cleared job-2 saved reader once, with its completion hook first."""
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
RUN = BASE/'saved-reader'
PREFIX = 'S11c_d_transverse_face_saved_reader'
WORKER = M/(PREFIX+'.py')
INPUT = M/(PREFIX+'_inputs.json')
GATE = M/(PREFIX+'_gate.json')
APPROVAL = M/'S11c_d_transverse_face_premise_authorization.json'
REVIEW = M/(PREFIX+'_review_r2_record.json')
PREPARATION = M/(PREFIX+'_review_r2_preparation.json')
GUARD = ROOT/'scripts/s11c_guarded_run.py'
SUPERVISOR = M/'S11c_d_end_normalization_run.py'
HOOK = ROOT/'scripts/codex_job_watch.py'
THREAD = '01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE = 'The already-approved job 2 saved-output reader has completed or reported an error.\nRoot: /var/projects/toy_physics/_scratch/s11c/s11c-d-transverse-face-20260928/saved-reader.\nInspect coordinator/launch, guard invocation/effective limits/limit validation/resource samples/child outcome/outcome and all logs; supervisor transverse_face_saved_reader invocation/active; strict worker stderr; stdout/checks identity; checks/failure/incompleteOperation; complete operation inputs/returns/index, artifact index, original source-object copies, launch-source routes/posthashes and all 168 input posthashes. Exit code is not scientific acceptance. Preserve every partial result or failure; no automatic retry.\n\nBoth fresh second-round reviewers literally CLEAR FOR THIS BOUNDED SAVED-OUTPUT READER. Exact worker SHA256 54ea52fd5aaea52eaca23314a227b6891be2856f7a9246697efe827042db1eaf and manifest d41990bb06196d830471445d8430ec9030be1641cffa51b99b843e25c547d0b5 remain unchanged from reviewed commit f9f5553e. Packet SHA256 2ef5f7e88cc2527c68098bb6613ad1982ba6b59cc208fbd6742ff768af4491ad. See reader review_r2_record/disposition and gate. This is build clearance, not physical acceptance. The first reader and both NEEDS REVISION reports remain at fcfbc9a2; no review rerun or external submission authorized.\n\nUser approved jobs 1-2 only, with necessary longer duration permitted. This reader uses ordinary 900s outer/840s native, 2GiB/zero swap/one CPU/nice15/32tasks/one native thread via unchanged scripts/s11c_guarded_run.py around S11c_d_end_normalization_run.py. No longer exception selected. No overlap, fallback or retry. Job 1 remains preserved: three 35s grade-extraction timeouts, 11 complete returns, 3 unresolved calls, zero map points and no face premise; 115.198955s, peak196616192bytes, zero swap/events, all33 inputs intact,98files/6040877bytes. Do not interpret unknowns as absent channels or no loss.\n\nThe reader restores only saved job1 inputs/returns and selected historical response/remainder/coefficient-system objects, plus prior readable end/forcing JSON. No producer import, binding/grade retry, root/mode/solve/transform/integral/current or new field/power work. Inspect actual branch/unit/operation/coverage joins, six selected frequency/end rows, historical contrast arrays and bookkeeping labels, supplied grade/incidence/unit conventions, saved open-flux scope, physical-input joins and retained candidate verdicts. Saved source identity and direct-retained-difference are not independent physical evidence; old halving vectors are not automatically a Born/full or physical-power anchor. No scientific payload restoration outside containment. Lightweight JSON/source/hash inspection is allowed.\n\nAfter inspection return actual job2 results/costs/coverage and omega1 versus omega3 pilot comparison, including consumed-candidate clearance debts and whether a longer unfinished-premise continuation has a concrete route. STOP before the one-day method-route decision, five-day authoring phase or jobs3-4; do not launch a continuation merely because longer runtime was permitted. Only jobs1-2 were approved. No outgoing field or loss magnitude supplied. Calibration of analog-light frequencies remains explicitly OPEN. No full Green/FORM/A11/A12, lossless-end or radiating-witness claim.\n\nPreserve all prior accepted results/failures/incident/review literals, Lean/S11_lean/shared guard and protected suffix f01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2. Canonical concise reports outside scratch; scratch never committed. Explicit diff exclusions for new generated tracked data above1MiB. Local hook targets this session01a0e01b-ef84-7192-817f-584cda5d339b; no model polling or recurring tasks.\n'


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
    assert gate['status'] == 'READY_FOR_ONE_GUARDED_TRANSVERSE_FACE_SAVED_READER_JOB'
    assert gate['scienceAuthorizationAfterCostedPlanGo'] is True
    assert gate['userScienceApprovalRecord'] == route(APPROVAL)
    approval = read(APPROVAL)
    assert approval['userReply'] == 'Approve jobs 1-2 only'
    assert approval['authorizedScienceJobOrdinals'] == [1, 2]
    assert approval['stopAfterJob2AndFrequencyComparison'] is True
    assert approval['methodRouteDecisionAuthorized'] is False
    assert approval['fieldPowerConstructionAuthorized'] is False
    assert gate['workerSha256'] == sha(WORKER)
    assert gate['inputManifestSha256'] == sha(INPUT)
    assert gate['independentBuildClearance'] is True and gate['scienceJobOrdinal'] == 2
    assert gate['reviewRecord'] == route(REVIEW)
    review = read(REVIEW)
    assert review['independentBuildClearance'] is True
    assert all(v['literalVerdict'] == 'CLEAR FOR THIS BOUNDED SAVED-OUTPUT READER'
               for v in review['reviews'].values()) and len(review['reviews']) == 2
    for path in (WORKER, INPUT):
        assert review['packetState']['fileHashes'][str(path.relative_to(ROOT))] == sha(path)
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
    assert spec['scienceJobOrdinal'] == 2 and spec['automaticRetry'] is False
    for join in spec['checkpointJoins']:
        value = read(Path(spec['inputs'][join['checkpointKey']]['path']))
        for component in join['jsonPath']:
            value = value[component]
        target = spec['inputs'][join['inputKey']]
        assert value['sha256'] == target['sha256'] and value['bytes'] == target['bytes']
    builder = (M/'S11c_d_sympy_builder_report.md').read_bytes()
    suffix = builder[builder.index(b'## Retained user-approved solver/export contract'):]
    assert hashlib.sha256(suffix).hexdigest() == gate['protectedBuilderSuffixSha256']
    return dict(workerSha256=sha(WORKER), inputManifestSha256=sha(INPUT),
                gateSha256=sha(GATE), launcherSha256=sha(Path(__file__)))


def command():
    return ['/usr/bin/python3', str(GUARD), '--log-directory', str(RUN/'resource-guard'),
            '--seconds', '900', '--memory-gib', '2', '--tasks-max', '32', '--',
            '/usr/bin/python3', str(SUPERVISOR), '--run-root', str(RUN),
            '--stage', 'transverse_face_saved_reader', '--', '/usr/bin/python3', '-u', str(WORKER),
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
