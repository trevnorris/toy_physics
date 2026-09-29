#!/usr/bin/env python3
"""Launch the exact cleared saved-reader continuation once, with its completion hook first."""
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
RUN = BASE/'saved-reader-continuation'
PREFIX = 'S11c_d_transverse_face_saved_reader_continue'
WORKER = M/(PREFIX+'.py')
INPUT = M/(PREFIX+'_inputs.json')
GATE = M/(PREFIX+'_gate.json')
APPROVAL = M/'S11c_d_transverse_face_premise_authorization.json'
REVIEW = M/(PREFIX+'_review_record.json')
PREPARATION = M/(PREFIX+'_review_preparation.json')
GUARD = ROOT/'scripts/s11c_guarded_run.py'
SUPERVISOR = M/'S11c_d_end_normalization_run.py'
HOOK = ROOT/'scripts/codex_job_watch.py'
CONTINUATION = M/(PREFIX+'_authorization.json')
CENSUS = M/(PREFIX+'_codec_census.json')
SCANNER = M/(PREFIX+'_pickle_metadata.py')
THREAD = '01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE = "The explicitly approved saved-reader continuation has completed or reported an error.\nRoot: /var/projects/toy_physics/_scratch/s11c/s11c-d-transverse-face-20260928/saved-reader-continuation.\nInspect actual coordinator/launch/guard invocation/effective limits/limit validation/resource samples/outcomes/logs, supervisor transverse_face_saved_reader_continue invocation/active, scientific stderr, stdout/checks identity, checks/failure/incompleteOperation, all journal operands/returns/indices, prior-argument/incomplete identity records, complete-summary-reuse and every source/launch posthash. Never accept science from exit code. Preserve all incomplete work and any failure; no automatic retry.\n\nBoth fresh reviewers literally CLEAR FOR THIS BOUNDED SAVED-READER CONTINUATION; exact worker ab3ff8d17ee94b6738d2716685ed32e1de96ba508e77355f971d53a9d0e489cb and manifest ab9ee0cdb5613411453afddf6ef8b5fe6903891b513082c7539943c5d2e0142b unchanged from8e18ddc1. Raw reports/packet68ad4e2e56a6163ac0506fb9fa5be8ec2809328592c2fcc60162783ade6780a0 remain in sibling saved-reader-continue-build-review. See continuation review_record/disposition and gate. Local static pickletools metadata of all20 prior blobs found no globals beyond the earlier census; no prior payload was restored outside containment. Metadata is not proof of unpickling. No optional edits or reviewer rerun.\n\nUser said Continue after job2's formatter failure, authorizing one narrow corrected continuation after independent clearance. This is job2 role but scientific execution3 under the four-execution ceiling. All9 complete prior returns must be RESTORED_PRIOR_COMPLETE_RETURN without calling their functions; source-join-restore-branch-context is the first resumed incomplete call. Verify the original20 operand/return blobs and3 complete summaries copied byte-for-byte, their actual argument/provenance matches, and the failed input bundle. Comparisons use reader-journal restorations, not necessarily byte-identical original source representations. Bare UndefinedFunction metadata must be inventoried without invoking it. No producer import, failed premise binding/grade/current/root/mode/solve/transform/integral or new field/power work.\n\nOne ordinary900s outer/840s native job under unchanged scripts/s11c_guarded_run.py around S11c_d_end_normalization_run.py;2GiB,zero swap,one CPU,nice15,32tasks,one native thread. No duration exception,overlap,fallback or retry. Scientific restoration only under containment; source/JSON/hash inspection is lightweight. All244 input pins, original failed71files/6020056bytes, earlier job1's98files and protected builder suffix f01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2 must remain unchanged.\n\nAfter actual saved evidence inspection, return results/costs/coverage and omega1/omega3 comparison, consumed-candidate clearance debts and concrete missing pieces, then STOP. Do not launch a premise-duration continuation, one-day method-route phase, five-day field/power build or another job. Job1 remains three35s source-binding timeouts,11complete returns,zero map points,no face premise. Old halving vectors are coefficient-polynomial bookkeeping, not automatically Born/full or physical-power anchors. Calibration of analog-light frequency is OPEN; no lossless-end,total loss,full Green/FORM/A11/A12 or radiating-witness claim.\n\nPreserve every accepted/failed result and historical review/incident debt. Leave Lean/S11_lean/shared guard untouched. Canonical status/reports outside scratch; scratch never committed. Add explicit diff exclusions for generated tracked data above1MiB. Hook targets session01a0e01b-ef84-7192-817f-584cda5d339b; no model polling or recurring task.\n"


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
    assert gate['status'] == 'READY_FOR_ONE_GUARDED_TRANSVERSE_FACE_SAVED_READER_CONTINUATION'
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
    assert gate['scientificExecutionOrdinal'] == 3
    assert gate['continuationAuthorizationRecord'] == route(CONTINUATION)
    continuation = read(CONTINUATION)
    assert continuation['approvedSavedReaderContinuation'] is True
    assert continuation['originalJobsApproval'] == route(APPROVAL)
    assert continuation['scientificExecutionOrdinal'] == 3 and continuation['maximumNewExecutions'] == 1
    assert continuation['scientificExecutionsUsedBefore'] == 2 and continuation['scientificExecutionCeiling'] == 4
    assert continuation['completedOperationReplayAuthorized'] is False
    assert continuation['unfinishedPremiseBindingAuthorized'] is False
    assert continuation['fieldPowerConstructionAuthorized'] is False
    assert continuation['stopAfterReaderAndFrequencyComparison'] is True
    assert gate['priorCodecCensus'] == route(CENSUS) and gate['priorCodecScanner'] == route(SCANNER)
    census = read(CENSUS)
    assert census['manifestSha256'] == sha(INPUT) and census['uniquePickles'] == 20
    assert census['scientificRestorations'] == census['reducersExecuted'] == 0
    assert RUN == Path(gate['runRoot']) and RUN/'complete' == Path(gate['resultDirectory'])
    assert gate['reviewRecord'] == route(REVIEW)
    review = read(REVIEW)
    assert review['independentBuildClearance'] is True
    assert all(v['literalVerdict'] == 'CLEAR FOR THIS BOUNDED SAVED-READER CONTINUATION'
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
    assert spec['scientificExecutionOrdinal'] == 3
    assert spec['inputs']['continuationAuthorization'] == route(CONTINUATION)
    assert not RUN.is_relative_to(Path(spec['priorRoot'])) and not RUN.is_relative_to(Path(spec['job1Root']))
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
            '--stage', 'transverse_face_saved_reader_continue', '--', '/usr/bin/python3', '-u', str(WORKER),
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
    paths.update((CONTINUATION, CENSUS, SCANNER, M/(PREFIX+'_review_disposition.md'), GATE, APPROVAL, Path(__file__).resolve(), REVIEW, PREPARATION, HOOK,
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
