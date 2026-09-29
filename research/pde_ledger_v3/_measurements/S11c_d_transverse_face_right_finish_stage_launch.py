#!/usr/bin/env python3
"""Launch the exact cleared 16 GiB RIGHT continuation once, with its completion hook first."""
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
RUN = BASE/'right-finish'
PREFIX = 'S11c_d_transverse_face_right_finish'
WORKER = M/(PREFIX+'.py')
INPUT = M/(PREFIX+'_inputs.json')
GATE = M/(PREFIX+'_gate.json')
APPROVAL = M/(PREFIX+'_authorization.json')
SCOPE = M/(PREFIX+'_scope.md')
CENSUS = M/(PREFIX+'_codec_census.json')
REVIEW = M/(PREFIX+'_review_r2_record.json')
PREPARATION = M/(PREFIX+'_review_r2_preparation.json')
GUARD = M/(PREFIX+'_guard.py')
SHARED_GUARD = ROOT/'scripts/s11c_guarded_run.py'
SUPERVISOR = M/'S11c_d_end_normalization_run.py'
HOOK = ROOT/'scripts/codex_job_watch.py'
THREAD = '01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE = 'The user-approved 16 GiB RIGHT-end fixed-point saved-return continuation has completed or reported an error.\nRoot: /var/projects/toy_physics/_scratch/s11c/s11c-d-transverse-face-fixed-point-20260929/right-finish.\nInspect actual coordinator/launch/stage-start; guard invocation/timeAuthorization/duration-validation/effective-limits/limit-validation/resource-samples/last-progress/child-outcome/outcome/logs; supervisor transverse_face_right_finish invocation/active; strict scientific stderr; stdout/checks byte identity; checks/failure/incompleteOperation; all journal/suboperation operands/returns/progress events, restored-suboperations, copy index, source snapshots and all1364 posthashes. Exit code is not scientific acceptance. Preserve incomplete work/failure, no automatic retry or new job.\n\nBoth corrected fresh reviews literally CLEAR FOR THIS 16 GIB RIGHT-END CONTINUATION. Exact worker f41ca1c84ce0d49ed35d413a36afe6fadab17049e3db3d6878839bec259fa708, manifest1ef4f620efb38de7235f8efec65c377aec36844a38f9e056c80b3d477e3425e3, task-local guard38278b15fbd27c7f5d36fb8f42b25c17ce60446c675dd326e0c642d500adad8f unchanged from3aba519b. Raw reports/packeta3f4b309a4af9761476c08ccd123c440e0e0a00df0ec555d092e542a1c716293 remain in sibling right-finish-build-review-r2. Earlier literal Claude NEEDS REVISION/Grok CLEAR and reviewed bytes remain at2d720c7d/right-finish-build-review. Accepted correction joins all3 helper hashes from the independent review to actual-file-checked gate. No optional code edits or review loops. Grok1431-byte stderr was reviewer CLI diagnostics, not scientific stderr. Source-only clearance does not accept runtime premise.\n\nUser expressly approved16GB, longer work while saved progress continues, and desktop-managed priority. This is execution6, one RIGHT-only saved-return continuation after execution5, not a reset. Both cgroup and native address-space cap16GiB; zero swap,one CPU,32tasks,one native thread,4GiB host available-memory reserve,global lock/overlap refusal and existing normalization supervisor retained. No nice requirement or scheduler exception/host configuration change. No overall wall/native deadline;3600s without distinct saved result stops natively, outer guard adds60s grace. Verify systemd RuntimeMaxUSec=infinity/Restart=no. CPU activity/input checkpoints are not progress. Shared scripts/s11c_guarded_run.py untouched. No fallback/restart/retry.\n\nVerify all23 prior COMPLETE operations restored without calling their functions, including original11 and complete REFERENCE/LEFT source/face/selected-seed results. All1311files/31272777bytes of execution5 (including earlier98/286-file archives) copied byte-for-byte. Only RIGHT-source-bind resumed, with15 preserved suboperation bundles. Actual original/fixed/continued/current arguments, source and summary ancestry join. Scalar52-54 input/returns must restore without Poly; input55 must restore and compare actual entry/epsilon before first new Poly(entry,epsilon_shape).nth(2). Entire recovered context and slab current restore. DomainMatrix/IntegerRing are the two explicit storage-codec additions. Check runtime decoding/summary equality and applyfunc order from actual joins, never static census alone. Scientific constructors/producers/completed calculations are not replayed.\n\nInspect actual source/carrier/grade/branch/unit/current/face/domain/control evidence and both selected RIGHT seeds16/17. Both required for selected support; unvisited/partial/nonfinite/unavailable stays unresolved. c2 is limited orientation-blind velocity/index evidence. Restore integrity/argument failures may be caught as LOCAL_OPERATION_EXCEPTION; treat those as integrity failures, not physics. Cgroup OOM/host reserve stop or late persistence error can omit usual completion records; missing receipt never means success. More memory/prospective review does not guarantee completion. Saved-return provenance is not independent revalidation.\n\nExecution5 stays preserved: selected REFERENCE/LEFT supported, RIGHT MemoryError in scalar55 before its face/seed checks;4830.204worker seconds including789.212priority suspension,peak2069409792bytes,zero swap/events,432 unchanged input pins. Its temporary scheduler exception was removed. No need for any further scheduler work. Keep all older failures/accepted outputs/incident/review debt, Lean/S11_lean and builder suffix f01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2 untouched.\n\nAfter completion inspect and return actual right-end result/cost/coverage and STOP. No further validator, retry, new method, omega3, full field or physical-power campaign is authorized. Saved omega1 finite current has no loss resolved at declared1e-6, not a physical bound; four baseline ratios are incident columns, not cases; pure-second-order omissions belong to retained series, not automatically full finite solve; regulator sensitivity is not its absolute effect or a subtractable floor. Total loss/analog-light calibration/full Green/FORM/A11/A12/radiating witness remain OPEN. Scientific restoration only inside containment; completion JSON/source/hash reads are lightweight. Canonical records outside scratch, runtime scratch never committed; add explicit diff exclusions for new generated tracked data above1MiB. Hook targets this session01a0e01b-ef84-7192-817f-584cda5d339b. No model polling/recurring tasks.\n'


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
    assert gate['status']=='READY_FOR_ONE_16GIB_RIGHT_END_CONTINUATION'
    assert gate['runtimeApproval']==gate['userScienceApprovalRecord']==route(APPROVAL)
    assert approval['status']=='USER_APPROVED_16GIB_RIGHT_END_CONTINUATION'
    assert approval['scientificExecutionOrdinal']==6 and approval['originalExecutionsUsed']==5
    assert approval['maximumNewScientificExecutions']==1 and approval['stopAfterResults'] is True
    assert approval['userReply']=="So be safe and give it 16 GB. I'll make sure nothing else runs on the box until this is complete."
    assert approval['runtimeReply']=='I give permission for running them longer than the normal limit. allow them to run as long as necessary as long as they are making progress'
    assert approval['priorityReply']=='Resume with desktop-managed priority'
    assert (approval['seconds'],approval['nativeSeconds'],approval['progressStallSeconds'],approval['guardBookkeepingGraceSeconds'])==(0,0,3600,60)
    assert approval['memoryBytes']==16*1024**3 and approval['priorityPolicy']=='DESKTOP_MANAGED'
    assert approval['automaticRetry'] is False
    pins=dict(workerSha256=sha(WORKER),inputManifestSha256=sha(INPUT),
              scopeSha256=sha(SCOPE),authorizationSha256=sha(APPROVAL),guardSha256=sha(GUARD),
              sharedGuardSha256=sha(SHARED_GUARD),supervisorSha256=sha(SUPERVISOR))
    assert all(gate[k]==v and review[k]==v for k,v in pins.items())
    assert gate['independentBuildClearance'] is True and review['independentBuildClearance'] is True
    assert gate['scienceJobOrdinal']==gate['scientificExecutionOrdinal']==6
    assert gate['maximumNewScientificExecutions']==1 and gate['reviewRecord']==route(REVIEW)
    assert len(review['reviews'])==2 and {leg['engine'] for leg in review['reviews']}=={'claude','grok'}
    for leg in review['reviews']:
        assert leg['literalVerdict']=='CLEAR FOR THIS 16 GIB RIGHT-END CONTINUATION'
        assert leg['verdict']=='CLEAR' and route(leg['output']['path'])==leg['output']
    for path in (WORKER,INPUT,GUARD,SHARED_GUARD,SUPERVISOR,SCOPE,APPROVAL):
        assert review['packetState']['fileHashes'][str(path.relative_to(ROOT))]==sha(path)
    assert gate['launcherSha256']==sha(Path(__file__))
    assert gate['completionMessageSha256']==hashlib.sha256(MESSAGE.encode()).hexdigest()
    assert gate['completionHookSha256']==sha(HOOK)
    assert (gate['seconds'],gate['nativeSeconds'],gate['progressStallSeconds'],gate['guardBookkeepingGraceSeconds'],gate['maximumEndBindings'])==(0,0,3600,60,1)
    assert (gate['memoryGiB'],gate['memoryBytes'],gate['tasksMax'],gate['swapMax'],gate['cpuCount'],gate['nativeThreads'])==(16,16*1024**3,32,0,1,1)
    assert gate['priorityPolicy']=='DESKTOP_MANAGED'
    assert gate['automaticRetry'] is False and gate['progressDependentRuntimeExplicitlyApproved'] is True
    assert gate['codecCensus']==route(CENSUS)
    assert review['codecCensusJoin']['currentManifestHashMatches'] is True
    assert review['codecCensusJoin']['allCensusedBytesMatch'] is True
    for item in read(PREPARATION)['sourceRecords']:
        assert sha(ROOT/item['path'])==item['sourceSha256'],item['path']
    assert spec['scientificExecutionOrdinal']==spec['scientificExecutionCeiling']==6
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
            '--seconds', '0', '--memory-gib', '16', '--tasks-max', '32', '--',
            '/usr/bin/python3', str(SUPERVISOR), '--run-root', str(RUN),
            '--stage', 'transverse_face_right_finish', '--', '/usr/bin/python3', '-u', str(WORKER),
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
                  GUARD, SHARED_GUARD, SUPERVISOR, M/(PREFIX+'_review_r2_disposition.md'),
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
