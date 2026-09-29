#!/usr/bin/env python3
"""Launch the exact cleared execution-4 fixed point once, with its completion hook first."""
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
BASE = ROOT/'_scratch/s11c/s11c-d-transverse-face-fixed-point-20260929'
RUN = BASE/'production'
PREFIX = 'S11c_d_transverse_face_fixed_point'
WORKER = M/(PREFIX+'.py')
INPUT = M/(PREFIX+'_inputs.json')
GATE = M/(PREFIX+'_gate.json')
APPROVAL = M/(PREFIX+'_authorization.json')
SCOPE = M/(PREFIX+'_scope.md')
CENSUS = M/(PREFIX+'_codec_census.json')
REVIEW = M/(PREFIX+'_review_r2_record.json')
PREPARATION = M/(PREFIX+'_review_r2_preparation.json')
GUARD = ROOT/'scripts/s11c_guarded_run.py'
SUPERVISOR = M/'S11c_d_end_normalization_run.py'
HOOK = ROOT/'scripts/codex_job_watch.py'
THREAD = '01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE = 'The explicitly approved execution 4 fixed-point transverse/current and limited reference-face job has completed or reported an error.\nRoot: /var/projects/toy_physics/_scratch/s11c/s11c-d-transverse-face-fixed-point-20260929/production.\nInspect actual coordinator/launch/stage-start, guard invocation/effective limits/limit validation/resource samples/child-outcome/outcome/logs, supervisor transverse_face_fixed_point invocation/active, strict scientific stderr, stdout/checks byte identity, checks/failure/incompleteOperation, journal and suboperation operands/returns/indices, all142 input posthashes and source snapshots. Exit status is not scientific acceptance. Preserve failures and incomplete work; no automatic retry.\n\nBoth fresh second-round reviewers literally CLEAR FOR THIS BOUNDED FIXED-POINT TRANSVERSE/FACE BUILD. Exact worker64cd6ea8ec7907edacd5fce9687c2072d3a235bd270c0ea483e8b76f18af9888 and manifest01dc1264f1ac59490880b5d186e6b52cf65f2afb288dd9a1054741e57a0fc41c remain unchanged from b8a1dd2c. Raw reports and packet58f12be41da26a4e28a564f97be8018286467cc42f41c5fa2508539f9d850624 remain in sibling build-review-r2; original NEEDS REVISION reports/source remain in build-review and d163c1fe/32e3e9fc. No new review or external submission. Review clearance is prospective, not a runtime result.\n\nUser approved execution4 only at omega1 and tangents(1/5,1/10), LAB_HELD/RHO4_CONSTANT: bind fixed frequency/tangents/materials before epsilon_shape-squared extraction; carry epsilon/eta/sigma where present through extraction; then declared end specialization. Ordinary shared scripts/s11c_guarded_run.py around existing S11c_d_end_normalization_run.py, 900s outer/840s native,180s per end binding,2GiB/zero swap/one CPU/nice15/32tasks/one native thread. No unlimited inheritance, overlap, fallback or retry. This uses the fourth and final scientific execution slot; after inspection STOP, no additional scientific job or validator is authorized.\n\nVerify all11 prior COMPLETE returns restored without calling their functions and all98 earlier premise files copied byte-for-byte, including3 incomplete input bundles. Original three35s timeouts/zero-map/no-face result and both saved-reader runs remain intact. New work is only three unfinished fixed-point bindings and limited face/current joins against seeds16/17. Inspect specialization receipts, carrier/source/branch/unit residuals, source-bind receipt-only JSON, real saved point/depth/current operands, per-drive applicability and responsive counts, native E-velocity omission movement, both selected seed statuses and shared LEFT/REFERENCE hash provenance. Both selected joins are required for SELECTED_PREMISE_SUPPORTED; partial, nonfinite, unavailable or unvisited remains unresolved. All-unloaded U controls do not demonstrate sensitivity. c2 is orientation-blind velocity-index evidence only. Unknown is not absent.\n\nThe static opcode census covered31 routes/29 unique files without restoration. Its older manifest hash was reconciled against identical current selected routes/bytes; duplicate LEFT/REFERENCE seeds explain representative-route differences. Runtime codec/structural equality and source-symbol compatibility remain actual checks, not assumed from metadata. Reviewers flagged remaining raw denominator/face/source-term costs; keep timeouts and saved suboperations instead of repairing/retrying.\n\nReport actual results/costs/coverage and stop. Also qualify the saved omega1 finite current balance: four quoted baseline ratios are four incident columns, not four physical cases; no loss resolved at declared1e-6 current resolution is not a certified physical bound. Pure-second-order interference omissions belong to the retained series and must not be automatically attributed to the full finite solve. Finite regulator sensitivity does not establish its absolute effect or a subtractable absorption floor. See S11c_d_transverse_face_fixed_point_saved_current_note.md. No total-loss magnitude, omega3 premise, all-grade lossless-end, calibrated analog-light, full Green/FORM/A11/A12 or radiating-witness claim. Calibration stays OPEN.\n\nScientific payload restoration remains inside containment; completion inspection uses JSON/source/hash metadata only. Preserve every accepted/failure artifact, incident/review history, Lean/S11_lean/shared guard and protected suffixf01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2. Canonical progress outside scratch; scratch never committed. Explicit diff exclusion for generated tracked data over1MiB. Hook targets THIS session01a0e01b-ef84-7192-817f-584cda5d339b. No model polling or recurring task.\n'


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
    assert gate['status'] == 'READY_FOR_ONE_GUARDED_TRANSVERSE_FACE_FIXED_POINT_JOB'
    assert gate['userScienceApprovalRecord'] == route(APPROVAL)
    approval = read(APPROVAL)
    assert approval['scientificExecutionOrdinal'] == 4 and approval['scientificExecutionsUsedBefore'] == 3
    assert approval['scientificExecutionCeiling'] == 4 and approval['maximumNewScientificExecutions'] == 1
    assert approval['stopAfterResults'] is True and approval['completedOperationReplayAuthorized'] is False
    assert approval['fullMapAuthorized'] is False and approval['firstOrderOutgoingFieldOrPowerAuthorized'] is False
    assert approval['automaticRetry'] is False and approval['durationException'] is False
    assert approval['frequency'] == '1' and approval['tangentialComponents'] == ['1/5','1/10']
    assert approval['localBindingSeconds'] == 180
    pins = dict(workerSha256=sha(WORKER), inputManifestSha256=sha(INPUT),
                scopeSha256=sha(SCOPE), authorizationSha256=sha(APPROVAL))
    assert all(gate[k] == v for k,v in pins.items())
    assert gate['independentBuildClearance'] is True and gate['scienceJobOrdinal'] == 4
    assert gate['scientificExecutionOrdinal'] == 4 and gate['scientificExecutionCeiling'] == 4
    assert gate['reviewRecord'] == route(REVIEW)
    review = read(REVIEW)
    assert review['independentBuildClearance'] is True and all(review[k] == v for k,v in pins.items())
    assert len(review['reviews']) == 2 and {v['engine'] for v in review['reviews']} == {'claude','grok'}
    for leg in review['reviews']:
        assert leg['literalVerdict'] == 'CLEAR FOR THIS BOUNDED FIXED-POINT TRANSVERSE/FACE BUILD'
        assert leg['verdict'] == 'CLEAR' and route(leg['output']['path']) == leg['output']
    for path in (WORKER, INPUT):
        assert review['packetState']['fileHashes'][str(path.relative_to(ROOT))] == sha(path)
    assert gate['launcherSha256'] == sha(Path(__file__))
    assert gate['completionMessageSha256'] == hashlib.sha256(MESSAGE.encode()).hexdigest()
    assert gate['guardSha256'] == sha(GUARD) and gate['supervisorSha256'] == sha(SUPERVISOR)
    assert gate['completionHookSha256'] == sha(HOOK)
    assert (gate['seconds'],gate['nativeSeconds'],gate['localBindingSeconds'],gate['maximumEndBindings']) == (900,840,180,3)
    assert (gate['memoryGiB'],gate['tasksMax'],gate['swapMax'],gate['cpuCount'],gate['nice'],gate['nativeThreads']) == (2,32,0,1,15,1)
    assert gate['automaticRetry'] is False and gate['durationException'] is False
    assert gate['codecCensus'] == route(CENSUS)
    assert review['codecCensusJoin']['allSelectedRoutesIdenticalBetweenManifests'] is True
    assert review['codecCensusJoin']['allCensusedBytesMatch'] is True
    for item in read(PREPARATION)['sourceRecords']:
        assert sha(ROOT/item['path']) == item['sourceSha256'], item['path']
    spec = read(INPUT)
    assert spec['scientificExecutionOrdinal'] == spec['scientificExecutionCeiling'] == 4
    assert spec['case'] == 'LAB_HELD_RHO4_CONSTANT' and spec['fixedPoint'] == dict(omega='1',tangents=['1/5','1/10'])
    assert spec['inputs']['fixedPointScope'] == route(SCOPE) and spec['inputs']['fixedPointAuthorization'] == route(APPROVAL)
    assert RUN == Path(spec['runRoot']) == Path(gate['runRoot'])
    assert RUN/'complete' == Path(gate['resultDirectory']) and not RUN.is_relative_to(Path(spec['priorRoot']))
    for name, record in spec['inputs'].items():
        assert route(record['path']) == record, name
    for join in spec['checkpointJoins']:
        value = read(Path(join['checkpoint']))['artifacts'][join['artifact']]
        assert value['sha256'] == spec['inputs'][join['inputKey']]['sha256']
    builder = (M/'S11c_d_sympy_builder_report.md').read_bytes()
    suffix = builder[builder.index(b'## Retained user-approved solver/export contract'):]
    assert hashlib.sha256(suffix).hexdigest() == gate['protectedBuilderSuffixSha256']
    return dict(**pins, gateSha256=sha(GATE), launcherSha256=sha(Path(__file__)))


def command():
    return ['/usr/bin/python3', str(GUARD), '--log-directory', str(RUN/'resource-guard'),
            '--seconds', '900', '--memory-gib', '2', '--tasks-max', '32', '--',
            '/usr/bin/python3', str(SUPERVISOR), '--run-root', str(RUN),
            '--stage', 'transverse_face_fixed_point', '--', '/usr/bin/python3', '-u', str(WORKER),
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
    paths.update((GATE, APPROVAL, SCOPE, CENSUS, Path(__file__).resolve(), REVIEW, PREPARATION, HOOK,
                  M/(PREFIX+'_review_r2_disposition.md'), M/(PREFIX+'_saved_current_note.md')))
    for path in sorted(paths):
        target = snapshot/path.relative_to(ROOT); target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target); assert sha(target) == sha(path)
    save(RUN/'source-snapshots.json', {str(p.relative_to(ROOT)): dict(source=route(p), snapshot=route(snapshot/p.relative_to(ROOT))) for p in sorted(paths)})
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
