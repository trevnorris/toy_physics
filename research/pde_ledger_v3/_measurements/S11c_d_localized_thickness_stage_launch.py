"""Launch one separately approved corrected thickness stage with its hook armed first."""
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
BASE = ROOT/'_scratch/s11c/s11c-d-localized-thickness-20260928'
RUN = BASE/'production'
M = ROOT/'research/pde_ledger_v3/_measurements'
GATE = M/'S11c_d_localized_thickness_gate.json'
INPUT = M/'S11c_d_localized_thickness_corrected_inputs.json'
WORKER = M/'S11c_d_localized_thickness_response.py'
THREAD = '01a0e01b-ef84-7192-817f-584cda5d339b'
PROPOSAL = M/'S11c_d_localized_thickness_proposal.json'
APPROVAL = M/'S11c_d_localized_thickness_authorization.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    with path.open('x') as f:
        json.dump(value, f, indent=2)
        f.write('\n')
        f.flush()
        os.fsync(f.fileno())


def preflight():
    gate = json.loads(GATE.read_text())
    assert gate['status'] == 'READY_FOR_ONE_GUARDED_LOCALIZED_THICKNESS_RESPONSE'
    assert gate['workerSha256'] == sha(WORKER)
    assert gate['inputManifestSha256'] == sha(INPUT)
    assert gate['adjudicationSha256'] == sha(Path(gate['adjudicationPath']))
    adjudication=json.loads(Path(gate['adjudicationPath']).read_text())
    assert adjudication['status']=='BOUNDED_FINDINGS_LOCALLY_DISPOSED_PENDING_STAGE_APPROVAL'
    assert adjudication['substantiveFindingsClosed'] is True
    assert adjudication['correctedWorkerSha256']==sha(WORKER)
    assert adjudication['correctedInputManifestSha256']==sha(INPUT)
    assert gate['scopeExplicitlyApproved'] is True and gate['substantiveReviewFindingsClosed'] is True
    assert gate['scopeApprovalSha256'] == sha(APPROVAL)
    assert gate['pendingGateSha256'] == sha(PROPOSAL)
    assert gate['launcherSha256'] == sha(Path(__file__))
    assert gate['completionMessageSha256'] == sha(M/'S11c_d_localized_thickness_completion_message.txt')
    assert gate['sharedGuardSha256'] == sha(ROOT/'scripts/s11c_guarded_run.py')
    assert gate['supervisorSha256'] == sha(M/'S11c_d_end_normalization_run.py')
    assert gate['completionHookSha256'] == sha(ROOT/'scripts/codex_job_watch.py')
    assert gate['reviewDispositionSha256'] == sha(Path(gate['reviewDispositionPath']))
    assert gate['nativeSeconds'] == 840 and gate['automaticRetry'] is False
    approval = json.loads(APPROVAL.read_text())
    assert approval['status'] == 'AUTHORIZED_ONE_GUARDED_LOCALIZED_THICKNESS_RESPONSE'
    assert approval['pendingGateSha256'] == sha(PROPOSAL)
    assert approval['workerSha256'] == sha(WORKER) and approval['inputManifestSha256'] == sha(INPUT)
    assert (gate['seconds'], gate['memoryGiB'], gate['tasksMax']) == (900, 2, 32)
    assert (gate['swapMax'], gate['cpuCount'], gate['nice'], gate['nativeThreads']) == (0, 1, 15, 1)
    for record in json.loads(INPUT.read_text())['inputs'].values():
        p = Path(record['path'])
        assert p.stat().st_size == record['bytes'] and str(p.resolve()) == record['canonicalPath']
        assert sha(p) == record['sha256'], str(p)
    return {'workerSha256': sha(WORKER), 'inputManifestSha256': sha(INPUT),
            'gateReceiptSha256': sha(GATE), 'launcherSha256': sha(Path(__file__))}


def command():
    return ['/usr/bin/python3', str(ROOT/'scripts/s11c_guarded_run.py'),
        '--log-directory', str(RUN/'resource-guard'), '--seconds', '900',
        '--memory-gib', '2', '--tasks-max', '32', '--',
        '/usr/bin/python3', str(M/'S11c_d_end_normalization_run.py'),
        '--run-root', str(RUN), '--stage', 'localized_thickness_response', '--',
        '/usr/bin/python3', '-u', str(WORKER), '--input-manifest', str(INPUT),
        '--gate-receipt', str(GATE), '--run-directory', str(RUN/'complete')]


def coordinate(descriptor):
    try:
        ready, _, _ = select.select([descriptor], [], [], 30)
        if not ready or os.read(descriptor, 1) != b'1':
            raise RuntimeError('Completion hook did not arm; no workload launched')
        os.close(descriptor)
        pins = preflight()
        assert pins == json.loads((RUN/'launch.json').read_text())['pins']
        started = time.monotonic()
        save(RUN/'stage-start.json', {'startedUtc': datetime.now(timezone.utc).isoformat(),
            'command': command(), 'pins': pins, 'automaticRetry': False})
        # Local blocking process wait; no model polling or retry.
        with (RUN/'guard-launch.stdout').open('x') as out, (RUN/'guard-launch.stderr').open('x') as err:
            result = subprocess.run(command(), cwd=ROOT, stdin=subprocess.DEVNULL,
                                    stdout=out, stderr=err)
        save(RUN/'coordinator-outcome.json', {'exitCode': result.returncode,
            'wallSeconds': time.monotonic()-started, 'scientificAcceptance': 'PENDING',
            'finishedUtc': datetime.now(timezone.utc).isoformat()})
        if result.returncode:
            (RUN/'issues.log').write_text('Guarded stage failed; inspect preserved records. No retry.\n')
        return result.returncode
    except BaseException:
        failure=traceback.format_exc()
        if not (RUN/'coordinator-outcome.json').exists():
            save(RUN/'coordinator-outcome.json',{'status':'COORDINATOR_ERROR_NO_RETRY','traceback':failure,
                'scientificAcceptance':'UNACCEPTED','finishedUtc':datetime.now(timezone.utc).isoformat()})
        (RUN/'issues.log').write_text(failure)
        raise


def launch():
    pins = preflight()
    assert shutil.which('codex'), 'Completion queue unavailable'
    RUN.mkdir(exist_ok=False)
    snapshot = RUN/'source';snapshot.mkdir()
    for path in (WORKER,INPUT,GATE,PROPOSAL,APPROVAL,Path(__file__),
                 M/'S11c_d_localized_thickness_completion_message.txt',
                 ROOT/'scripts/s11c_guarded_run.py',M/'S11c_d_end_normalization_run.py',
                 ROOT/'scripts/codex_job_watch.py',
                 ROOT/'research/pde_ledger_v3/directives/S11c_d_localized_thickness_implementation.md',
                 M/'S11c_d_localized_thickness_review_disposition.md',
                 M/'S11c_d_localized_thickness_review_adjudication.json',
                 M/'S11c_d_localized_thickness_phase_source_evidence.json',
                 M/'S11c_d_localized_thickness_claude_review.md',
                 M/'S11c_d_localized_thickness_grok_review.md'):
        target=snapshot/path.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(path,target)
        assert sha(target)==sha(path)
    read_fd, write_fd = os.pipe()
    with (RUN/'coordinator.stdout').open('x') as out, (RUN/'coordinator.stderr').open('x') as err:
        job = subprocess.Popen(['/usr/bin/python3', str(Path(__file__)), '--coordinator', str(read_fd)],
            cwd=ROOT, stdin=subprocess.DEVNULL, stdout=out, stderr=err,
            pass_fds=(read_fd,), start_new_session=True)
    os.close(read_fd)
    try:
        with (RUN/'watcher.stdout').open('x') as out, (RUN/'watcher.stderr').open('x') as err:
            watcher = subprocess.Popen(['/usr/bin/python3', str(ROOT/'scripts/codex_job_watch.py'),
                '--pid', str(job.pid), '--thread', THREAD,
                '--directory', str(RUN/'completion-watcher'),
                '--message-file', str(M/'S11c_d_localized_thickness_completion_message.txt'),
                '--error-log', str(RUN/'issues.log')], cwd=ROOT,
                stdin=subprocess.DEVNULL, stdout=out, stderr=err, start_new_session=True)
        record = {'jobPid': job.pid, 'watcherPid': watcher.pid, 'pins': pins,
            'command': command(), 'thread': THREAD, 'maximumScientificWorkers': 1,
            'completionMessageSha256': sha(M/'S11c_d_localized_thickness_completion_message.txt')}
        save(RUN/'launch.json', record)
        for _ in range(50):
            path = RUN/'completion-watcher/state.json'
            if path.exists() and json.loads(path.read_text())['status'] == 'waiting':
                os.write(write_fd, b'1')
                print(json.dumps({'jobPid': job.pid, 'watcherPid': watcher.pid,
                    'hookStatus': 'waiting', 'runDirectory': str(RUN)}))
                return
            if watcher.poll() is not None or job.poll() is not None:
                raise RuntimeError('Coordinator or hook exited before arming')
            time.sleep(.1)
        raise RuntimeError('Completion hook did not arm')
    finally:
        # EOF releases a waiting coordinator without launching if arming failed.
        os.close(write_fd)



def activate_gate():
    # This file is written only after a new explicit science-stage approval.
    # Missing approval stops before creating the gate, run directory or child.
    approval=json.loads(APPROVAL.read_text())
    proposal=json.loads(PROPOSAL.read_text())
    assert not GATE.exists() and not RUN.exists(), 'no retry or overwrite'
    assert approval['status']=='AUTHORIZED_ONE_GUARDED_LOCALIZED_THICKNESS_RESPONSE'
    assert approval['pendingGateSha256']==sha(PROPOSAL)
    assert approval['workerSha256']==proposal['workerSha256']==sha(WORKER)
    assert approval['inputManifestSha256']==proposal['inputManifestSha256']==sha(INPUT)
    assert proposal['launcherSha256']==sha(Path(__file__))
    assert proposal['status']=='PENDING_EXPLICIT_APPROVAL_FOR_ONE_GUARDED_LOCALIZED_THICKNESS_RESPONSE'
    assert proposal['scopeExplicitlyApproved'] is False and proposal['substantiveReviewFindingsClosed'] is True
    gate={**proposal,'status':'READY_FOR_ONE_GUARDED_LOCALIZED_THICKNESS_RESPONSE',
          'scopeExplicitlyApproved':True,'scopeApprovalSha256':sha(APPROVAL),'scopeApprovalPath':str(APPROVAL),
          'pendingGateSha256':sha(PROPOSAL),'activatedUtc':datetime.now(timezone.utc).isoformat()}
    save(GATE,gate)
    preflight()


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--coordinator':
        sys.exit(coordinate(int(sys.argv[2])))
    if len(sys.argv) != 1:
        raise ValueError('Unexpected launch arguments')
    activate_gate()
    launch()
