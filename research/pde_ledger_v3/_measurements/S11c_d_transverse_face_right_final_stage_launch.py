#!/usr/bin/env python3
"""Launch the final RIGHT continuation with its explicit minor-fix authority once, with its completion hook first."""
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
RUN = BASE/'right-final'
PREFIX = 'S11c_d_transverse_face_right_final'
WORKER = M/(PREFIX+'.py')
INPUT = M/(PREFIX+'_inputs.json')
GATE = M/(PREFIX+'_gate.json')
APPROVAL = M/(PREFIX+'_authorization.json')
SCOPE = M/(PREFIX+'_scope.md')
CENSUS = M/'S11c_d_transverse_face_right_finish_codec_census.json'
MINOR = M/(PREFIX+'_minor_fix_authorization.json')
REPAIR = M/(PREFIX+'_minor_repair_record.json')
REVIEW = M/(PREFIX+'_review_record.json')
PREPARATION = M/(PREFIX+'_review_preparation.json')
GUARD = M/(PREFIX+'_guard.py')
SHARED_GUARD = ROOT/'scripts/s11c_guarded_run.py'
SUPERVISOR = M/'S11c_d_end_normalization_run.py'
HOOK = ROOT/'scripts/codex_job_watch.py'
THREAD = '01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE = 'The user-authorized final RIGHT continuation has completed or reported an error.\nRoot: /var/projects/toy_physics/_scratch/s11c/s11c-d-transverse-face-fixed-point-20260929/right-final.\nInspect coordinator/launch/stage-start, guard invocation/timeAuthorization/duration-validation/effective-limits/limit-validation/resource samples/last-progress/outcomes and logs, supervisor transverse_face_right_final invocation/active, scientific stderr, stdout/checks byte identity, checks/failure/incomplete operation, all complete operands/returns/indices/suboperations/progress/copy/source snapshots and1370 input posthashes. Process exit is not acceptance. Preserve failure and STOP; no retry, standalone diagnostic, validator or further tooling job.\n\nExecution7 is the one final RIGHT continuation after six historical executions. User subsequently said: If it only needs minor fixes then fix it then run it. Both reviewed reports finished before the narrow correction; literal Claude NEEDS REVISION and Grok CLEAR FOR THIS FINAL RIGHT-END CONTINUATION are preserved, not relabelled. Reviewed baselinef4a4f029/packet6a5329f3bc998142e96fc863cb1c0b0012962b70155c985e5fd25419e57041ae remains in sibling right-final-build-review. No fresh independent CLEAR or second review round. Gate honestly records user-authorized minor serialization repair, with independentBuildClearance=false.\n\nClaude identified the unchanged summary equality risk. Metadata/source inspection found saved frequencyDomain strings 1 > 0 and StrictGreaterThan/NEWOBJ without evaluate=False; installed SymPy default construction evaluates. The minor correction adds an exact StrictGreaterThan loader adapter that returns the original class with the actual saved arguments and evaluate=False. No replacement truth value, expression normalization or nonfatal summary fallback. Stdlib synthetic NEWOBJ tests passed; scientific payloads were not restored during repair. Actual runtime joins remain required; this is not proof every summary difference is cured. See minor_repair_record/checks/evidence/authority.\n\nRestore23 completed execution5 returns without their functions and15 RIGHT suboperation bundles, including scalar52–54 input/returns and failed input55. Preserve1311-file tree and all1441 execution6 files in place. Only RIGHT source/current/face/selected16/17 work is new. Exact epsilon-only Add/Mul/nonnegative-integer-power coefficient convolution replaces new Poly domain inference; epsilon-independent expressions stay opaque, eta/sigma live, full reconstruction/current guards unchanged. Both old and fresh JSON summaries and differing paths persist before any remaining mismatch stops. Verify source/symbol/argument joins; do not interpret an integrity exception as physics.\n\nSame omega1,tangents(1/5,1/10),LAB_HELD/RHO4_CONSTANT.16GiB cgroup AND native address-space cap,zero swap,one CPU,32tasks,one thread,4GiB host reserve,desktop-managed priority,no scheduler changes,global lock/overlap refusal and existing normalization supervisor. Shared guard unchanged. No total deadline;3600s without distinct completed saved result stops,outer60s grace. Verify actual systemd RuntimeMaxUSec=infinity/Restart=no. Input receipts or CPU activity are not progress. No unguarded fallback,restart or retry.\n\nReport actual right-end result,costs,coverage and STOP. If RIGHT does not complete, retain selected REFERENCE/LEFT support and RIGHT unresolved; no more tooling cycle. Both selected RIGHT seeds16/17 required for selected support, c2 remains limited velocity/index evidence. No all-grade matched lossless ends,total loss,omega3,forced field,power,calibrated analog-light,Green/FORM/A11/A12 or radiating-witness claim. Prior finite omega1 current deficit unresolved at1e-6 is not a physical bound. Preserve historical failures/review debts/incident records,Lean/S11_lean/shared guard and suffixf01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2. Scientific restoration stays contained; completion inspection JSON/hash/source only. Scratch never committed. Hook targets01a0e01b-ef84-7192-817f-584cda5d339b. No model polling or recurring task.\n'


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
    repair=read(REPAIR);minor=read(MINOR)
    assert gate['status']=='READY_FOR_ONE_FINAL_RIGHT_END_CONTINUATION'
    assert gate['runtimeApproval']==gate['userScienceApprovalRecord']==route(APPROVAL)
    assert approval['status']=='USER_DIRECTED_FINAL_RIGHT_END_CONTINUATION'
    assert approval['scientificExecutionOrdinal']==7 and approval['originalExecutionsUsed']==6
    assert approval['maximumNewScientificExecutions']==1 and approval['stopAfterResults'] is True
    assert approval['finalAttempt'] is True and approval['maximumReviewRounds']==1
    assert (approval['seconds'],approval['nativeSeconds'],approval['progressStallSeconds'],approval['guardBookkeepingGraceSeconds'])==(0,0,3600,60)
    assert approval['memoryBytes']==16*1024**3 and approval['priorityPolicy']=='DESKTOP_MANAGED'
    assert approval['automaticRetry'] is False
    pins=dict(workerSha256=sha(WORKER),inputManifestSha256=sha(INPUT),
              scopeSha256=sha(SCOPE),authorizationSha256=sha(APPROVAL),guardSha256=sha(GUARD),
              sharedGuardSha256=sha(SHARED_GUARD),supervisorSha256=sha(SUPERVISOR))
    assert all(gate[k]==v for k,v in pins.items())
    assert all(review[k]==v for k,v in pins.items() if k!='workerSha256')
    assert gate['independentBuildClearance'] is False and review['independentBuildClearance'] is False
    assert gate['executionApproval']=='USER_AUTHORIZED_MINOR_SERIALIZATION_REPAIR'
    assert gate['minorRepairRecord']==route(REPAIR) and gate['minorFixAuthority']==route(MINOR)
    assert repair['status']=='MINOR_SERIALIZATION_REPAIR_VERIFIED_USER_AUTHORIZED'
    assert repair['reviewedWorkerSha256']==review['workerSha256'] and repair['workerSha256']==pins['workerSha256']
    assert repair['onlyCodecAndAuthorizationChanged'] is True and repair['independentBuildClearance'] is False
    assert repair['reviewRecord']==route(REVIEW) and repair['authority']==route(MINOR)
    assert minor['userReply']=='If it only needs minor fixes then fix it then run it'
    assert minor['reviewedPacketSha256']==review['packetSha256']
    assert minor['maximumNewScientificExecutions']==1 and minor['stopAfterResults'] is True
    for key in ('testRecord','testSource','sourceEvidence'):assert route(repair[key]['path'])==repair[key]
    assert read(Path(repair['testRecord']['path']))['workerSha256']==pins['workerSha256']
    for item in read(Path(repair['sourceEvidence']['path']))['sourceExcerpts']:
        assert route(item['source']['path'])==item['source']
    assert gate['scienceJobOrdinal']==gate['scientificExecutionOrdinal']==7 and gate['finalAttempt'] is True
    assert gate['maximumNewScientificExecutions']==1 and gate['reviewRecord']==route(REVIEW)
    assert len(review['reviews'])==2 and {leg['engine'] for leg in review['reviews']}=={'claude','grok'}
    for leg in review['reviews']:
        assert leg['verdict']=={'claude':'NEEDS REVISION','grok':'CLEAR'}[leg['engine']]
        assert route(leg['output']['path'])==leg['output']
    assert gate['launcherSha256']==sha(Path(__file__))
    assert gate['completionMessageSha256']==hashlib.sha256(MESSAGE.encode()).hexdigest()
    assert gate['completionHookSha256']==sha(HOOK)
    assert (gate['seconds'],gate['nativeSeconds'],gate['progressStallSeconds'],gate['guardBookkeepingGraceSeconds'],gate['maximumEndBindings'])==(0,0,3600,60,1)
    assert (gate['memoryGiB'],gate['memoryBytes'],gate['tasksMax'],gate['swapMax'],gate['cpuCount'],gate['nativeThreads'])==(16,16*1024**3,32,0,1,1)
    assert gate['priorityPolicy']=='DESKTOP_MANAGED'
    assert gate['automaticRetry'] is False and gate['progressDependentRuntimeExplicitlyApproved'] is True
    for item in read(PREPARATION)['sourceRecords']:
        path=ROOT/item['path']
        if path==WORKER:assert sha(path)==repair['workerSha256'] and item['sourceSha256']==repair['reviewedWorkerSha256']
        else:assert sha(path)==item['sourceSha256'],item['path']
    assert spec['scientificExecutionOrdinal']==spec['scientificExecutionCeiling']==7
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
            '--stage', 'transverse_face_right_final', '--', '/usr/bin/python3', '-u', str(WORKER),
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
                  MINOR, REPAIR, M/(PREFIX+'_minor_repair_checks.json'),
                  M/(PREFIX+'_minor_repair_test.py'), M/(PREFIX+'_minor_repair_evidence.json'),
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
