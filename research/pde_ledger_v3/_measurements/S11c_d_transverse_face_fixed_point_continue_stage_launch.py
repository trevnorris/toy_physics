#!/usr/bin/env python3
"""Launch the exact cleared progress-dependent fixed-point continuation once, with its completion hook first."""
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
RUN = BASE/'continuation'
PREFIX = 'S11c_d_transverse_face_fixed_point_continue'
WORKER = M/(PREFIX+'.py')
INPUT = M/(PREFIX+'_inputs.json')
GATE = M/(PREFIX+'_gate.json')
APPROVAL = M/(PREFIX+'_authorization.json')
SCOPE = M/(PREFIX+'_scope.md')
CENSUS = M/(PREFIX+'_codec_census.json')
REVIEW = M/(PREFIX+'_review_record.json')
PREPARATION = M/(PREFIX+'_review_preparation.json')
GUARD = M/(PREFIX+'_guard.py')
SHARED_GUARD = ROOT/'scripts/s11c_guarded_run.py'
SUPERVISOR = M/'S11c_d_end_normalization_run.py'
HOOK = ROOT/'scripts/codex_job_watch.py'
THREAD = '01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE = "The user-approved progress-dependent fixed-point saved-return continuation has completed or reported an error.\nRoot: /var/projects/toy_physics/_scratch/s11c/s11c-d-transverse-face-fixed-point-20260929/continuation.\nInspect actual coordinator/launch/stage-start; guard invocation/timeAuthorization/duration-validation/effective-limits/limit-validation/resource-samples/last-progress/child-outcome/outcome/logs; supervisor transverse_face_fixed_point_continue invocation/active; strict scientific stderr; stdout/checks identity; checks/failure/incompleteOperation; all journal/suboperation operands/returns and progress events, restored-suboperations, copy indices, source snapshots and all432 posthashes. Exit status is not scientific acceptance. Preserve incomplete work and failures, no automatic retry.\n\nBoth fresh reviewers literally CLEAR FOR THIS PROGRESS-DEPENDENT FIXED-POINT CONTINUATION. Reviewed worker5d38c329e8e15cff091c09584c1d33806fd8318f66a779fc10f2ae52d00e3135, manifest e78fc287fec25006d555ffb107d1fdcc6ef222554b3221d546c11331d05bb63b and task-local guard1c2bef4f691555e98e20eec1ccb76c8861cfd92ce827ddc36caa19afc12f5ecd are unchanged from3b831062. Raw reports/packet c88b13119288be87d01ce40eed1c9e3ada483d30015fdebad3548feb555197f5 remain in sibling continuation-build-review. Source-only clearance, not a runtime premise. Optional hardening was recorded but no edits/review loop made. Claude's literal source-only disclaimer is preserved. Grok's8967-byte stderr is reviewer diagnostics, not scientific stderr.\n\nUser explicitly allowed longer running while progress is made. This is one same-scope continuation, scientific execution5 after the original four, not a reset of counts. No total wall/native cap; native3600s without a new distinct saved result stops, outer guard3660s includes60s bookkeeping grace. Verify actual systemd RuntimeMaxUSec=infinity/Restart=no. CPU activity/started work and *-input checkpoints do not count. All other controls stay2GiB/zero swap/one CPU/nice15/32tasks/one native thread, host-memory checks, global lock and overlap refusal, native address-space limit and existing normalization supervisor. Shared scripts/s11c_guarded_run.py is unchanged; task-local continue_guard.py supplies only this approved duration/progress exception. No fallback, restart or retry.\n\nVerify all11 original COMPLETE returns restored without calling functions, all98 earlier premise files and all286 execution4 files/8104854bytes copied byte-for-byte, and33 saved suboperation bundles restored as reached. Actual original/fixed/current argument joins and source/branch/unit joins are required. REFERENCE/LEFT resumed interface-power carrier extraction; RIGHT resumed bulk-normal current density. Current scalar operation remains Poly(entry,epsilon_shape).nth(2), with per-entry attempted-input and coefficient returns. Some explicitly listed unsaved small context is reconstructed, not falsely labelled restored. Preserve the source-specialization/transverse/grade/residual/raw/control evidence and actual selected seeds16/17. Both selected joins required for selected support; partial/nonfinite/unvisited remains unresolved. Reference/c2 evidence is limited, not all-grade matched-end or power acceptance. A saved carrier is provenance-linked, not independently revalidated; pin/copy failures inside restoration may be local exceptions, which must be treated as integrity failures during inspection. Malformed progress can leave no child-outcome; missing outcome never implies success.\n\nThe previous execution4 remains preserved:553.558worker seconds,peak509546496bytes,zero swap/events,142 inputs and86 snapshots intact,three180s carrier timeouts,zero point/seed/face results. Earlier failed/successful reader and original premise remain intact. Report actual continued result/cost/coverage and STOP; no additional scientific validator, retry, method-route campaign, omega3 or field/power build is authorized. Saved omega1 finite current has no loss resolved at declared1e-6, not a physical bound. Four baseline ratios are incident columns, not four cases. Pure-second-order omissions belong to retained series, not automatically full finite solve; regulator sensitivity is not its absolute effect or a subtractable floor. Total loss/analog-light calibration/full Green/FORM/A11/A12/radiating witness remain OPEN.\n\nScientific payload restoration only inside containment; completion source/JSON/hash metadata inspection is lightweight. Preserve all accepted/failure/incident/review history, Lean/S11_lean/shared guard and protected suffix f01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2. Canonical records outside scratch; scratch never committed. Explicit diff exclusions for new generated tracked data above1MiB. Hook to THIS session01a0e01b-ef84-7192-817f-584cda5d339b. No model polling or recurring task.\n"


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
    gate=read(GATE);approval=read(APPROVAL);spec=read(INPUT);review=read(REVIEW)
    assert gate['status']=='READY_FOR_ONE_PROGRESS_GUARDED_FIXED_POINT_CONTINUATION'
    assert gate['runtimeApproval']==gate['userScienceApprovalRecord']==route(APPROVAL)
    assert approval['status']=='USER_APPROVED_PROGRESS_DEPENDENT_FIXED_POINT_CONTINUATION'
    assert approval['scientificExecutionOrdinal']==5 and approval['originalExecutionsUsed']==4
    assert approval['maximumNewScientificExecutions']==1 and approval['stopAfterResults'] is True
    assert approval['userReply']=='I give permission for running them longer than the normal limit. allow them to run as long as necessary as long as they are making progress'
    assert (approval['seconds'],approval['nativeSeconds'],approval['progressStallSeconds'],approval['guardBookkeepingGraceSeconds'])==(0,0,3600,60)
    assert approval['automaticRetry'] is False
    pins=dict(workerSha256=sha(WORKER),inputManifestSha256=sha(INPUT),
              scopeSha256=sha(SCOPE),authorizationSha256=sha(APPROVAL),guardSha256=sha(GUARD))
    assert all(gate[k]==v and review[k]==v for k,v in pins.items())
    assert gate['independentBuildClearance'] is True and review['independentBuildClearance'] is True
    assert gate['scienceJobOrdinal']==gate['scientificExecutionOrdinal']==5
    assert gate['maximumNewScientificExecutions']==1 and gate['reviewRecord']==route(REVIEW)
    assert len(review['reviews'])==2 and {leg['engine'] for leg in review['reviews']}=={'claude','grok'}
    for leg in review['reviews']:
        assert leg['literalVerdict']=='CLEAR FOR THIS PROGRESS-DEPENDENT FIXED-POINT CONTINUATION'
        assert leg['verdict']=='CLEAR' and route(leg['output']['path'])==leg['output']
    for path in (WORKER,INPUT,GUARD,SCOPE,APPROVAL):
        assert review['packetState']['fileHashes'][str(path.relative_to(ROOT))]==sha(path)
    assert gate['launcherSha256']==sha(Path(__file__))
    assert gate['completionMessageSha256']==hashlib.sha256(MESSAGE.encode()).hexdigest()
    assert gate['sharedGuardSha256']==sha(SHARED_GUARD) and gate['supervisorSha256']==sha(SUPERVISOR)
    assert gate['completionHookSha256']==sha(HOOK)
    assert (gate['seconds'],gate['nativeSeconds'],gate['progressStallSeconds'],gate['guardBookkeepingGraceSeconds'],gate['maximumEndBindings'])==(0,0,3600,60,3)
    assert (gate['memoryGiB'],gate['tasksMax'],gate['swapMax'],gate['cpuCount'],gate['nice'],gate['nativeThreads'])==(2,32,0,1,15,1)
    assert gate['automaticRetry'] is False and gate['progressDependentRuntimeExplicitlyApproved'] is True
    assert gate['codecCensus']==route(CENSUS)
    assert review['codecCensusJoin']['currentManifestHashMatches'] is True
    assert review['codecCensusJoin']['allCensusedBytesMatch'] is True
    for item in read(PREPARATION)['sourceRecords']:
        assert sha(ROOT/item['path'])==item['sourceSha256'],item['path']
    assert spec['scientificExecutionOrdinal']==spec['scientificExecutionCeiling']==5
    assert spec['case']=='LAB_HELD_RHO4_CONSTANT' and spec['fixedPoint']==dict(omega='1',tangents=['1/5','1/10'])
    assert spec['inputs']['continuationScope']==route(SCOPE) and spec['inputs']['continuationAuthorization']==route(APPROVAL)
    assert RUN==Path(spec['runRoot'])==Path(gate['runRoot'])
    assert RUN/'complete'==Path(gate['resultDirectory']) and not RUN.is_relative_to(Path(spec['resumeRoot']))
    for name,record in spec['inputs'].items():assert route(record['path'])==record,name
    for join in spec['checkpointJoins']:
        value=read(Path(join['checkpoint']))['artifacts'][join['artifact']]
        assert value['sha256']==spec['inputs'][join['inputKey']]['sha256']
    expected_command=command();split=expected_command.index('--')
    assert gate['guardedCommand']==expected_command[split+1:]
    builder=(M/'S11c_d_sympy_builder_report.md').read_bytes()
    suffix=builder[builder.index(b'## Retained user-approved solver/export contract'):]
    assert hashlib.sha256(suffix).hexdigest()==gate['protectedBuilderSuffixSha256']
    return dict(**pins,gateSha256=sha(GATE),launcherSha256=sha(Path(__file__)))


def command():
    return ['/usr/bin/python3', str(GUARD), '--log-directory', str(RUN/'resource-guard'),
            '--seconds', '0', '--memory-gib', '2', '--tasks-max', '32', '--',
            '/usr/bin/python3', str(SUPERVISOR), '--run-root', str(RUN),
            '--stage', 'transverse_face_fixed_point_continue', '--', '/usr/bin/python3', '-u', str(WORKER),
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
                  GUARD, SHARED_GUARD, SUPERVISOR, M/(PREFIX+'_review_disposition.md'),
                  M/'S11c_d_transverse_face_fixed_point_saved_current_note.md'))
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
