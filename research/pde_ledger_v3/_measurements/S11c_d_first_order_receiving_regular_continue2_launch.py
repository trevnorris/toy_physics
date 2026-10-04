#!/usr/bin/env python3
"""One standing-authorized source-projection and real-axis receiving prerequisite; hook first, no retry."""
from datetime import datetime, timezone
import hashlib,json,os,select,shutil,subprocess,sys,time,traceback,runpy
from pathlib import Path
ROOT=Path('/var/projects/toy_physics')
M=ROOT/'research/pde_ledger_v3/_measurements'
PREFIX='S11c_d_first_order_receiving_regular_continue2'
GATE=M/(PREFIX+'_gate.json')
MANIFEST=M/(PREFIX+'_inputs.json')
RUN=Path('/var/projects/toy_physics/_scratch/s11c/s11c-defect-packet-action-20261003/first-order-receiving-regular-continuation-02')
HOOK=ROOT/'scripts/codex_job_watch.py'
THREAD='01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE='The user-authorized saved-prefix REAL-AXIS RECEIVING polynomial-representation continuation has completed or reported an error.\nRoot: /var/projects/toy_physics/_scratch/s11c/s11c-defect-packet-action-20261003/first-order-receiving-regular-continuation-02.\nInspect actual coordinator/guard/supervisor commands, no-deadline resources, strict stderr/stdout-checks, all copied/prior/source posthashes, complete original and restored arguments and every new polynomial/GCD/Bezout/Sturm/domain/determinant/adjugate/growth/control return. JSON/source/hash only outside containment; no scientific restoration, validator or automatic retry. Exit is not acceptance.\nOriginal failure211ec378 and first continuation failure678dd6f0 stay immutable. Latest prefix has99artifacts/12newliteralBezoutzeros and12constant-ray exclusions for entries0/1/2, no determinant;12.918s,peak211738624B,zeroevents/swap,11203posthashes intact. Original639zeros/171inheritedzeros/2controls/source/movingbasis/nativecrosscurrents remain saved. All4761latest-treefiles copied. No completed scientific prefix replay.\nThe task-local correction saves exact polynomial original/native coefficients/reconstruction/expanded operands then checks an exact zero residual instead of structurally comparing grouped complex coefficients with expanded terms. No tolerance or parameter/equation/ray/Sturm/root decision change. Actual unpublished failed reconstruction is not claimed as restored or passed. Source data determine the result. Restore published Cq/40forcing factors and all12completedray/proof operands without calculation; new work starts at entry3base0 representation check. Complete remaining scientific tail unchanged. independentBuildClearance=false records tested local-tooling authority; historical ClaudeCLEAR/partial coverage is not fresh clearance.\nUnknown/nonzero/pole stops remain. A pass gives scoped C3 real-axis regularity/conditional weightedL2 and L1decay, not a full five-field inverse, far-bulk flux or leakage. Reflection is survival. Nonuniform physical work including exterior acoustic flux, memory, mass-rate/chemical and LAB_HELD work remains required next. Mixed recovery PARKED; no general verification expansion.\nStanding Claude-only user authority covers faithful bounded physics continuation without routine questions. Same realomega3 restbulkLAB_HELD/RHO4,savededges,cs=sqrt6/2,arbitrary doublet. No field/Fourier/power integral,current/leakage/loss/sweep/drain/calibration. Shared no-deadline guard4GiBnative/cgroup within16GiBpool,zero swap,oneCPU/thread,32tasks,4GiBreserve,desktoppriority,RuntimeMaxUSec=infinity/Restart=no. HookFIRST01a0e01b-ef84-7192-817f-584cda5d339b. No fallback/schedulerchange/automaticretry/modelpolling/recurringtask. Preserve allhistory/Lean/S11_lean/sharedguard/protectedbuilder;scratchnevercommitted.\n'

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def save(path,value):
    with path.open('x') as f:
        json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
def verify():
    gate,manifest=read(GATE),read(MANIFEST)
    assert gate['status']=='READY_FOR_ONE_FIRST_ORDER_RECEIVING_REGULAR_CONTINUATION2'
    assert gate['pooledExecution'] is True  # Exact assessment/repair is checked by verify_gate below.
    assert gate['launcher']==manifest['launcher']==str(Path(__file__).resolve())
    assert gate['sharedGuard']==str(ROOT/'scripts/s11c_guarded_run.py')
    assert gate['supervisor']==str(M/'S11c_d_end_normalization_run.py')
    assert gate['launcherSha256']==sha(__file__)
    assert gate['sourcePins']==manifest['sourcePins']
    assert gate['completionMessageSha256']==hashlib.sha256(MESSAGE.encode()).hexdigest()
    assert gate['outputDirectory']==str(RUN/'complete')
    assert manifest['resources']['memoryBytes']==4*1024**3 and manifest['resources']['durationLimits'] is None
    # Standard-library gate definitions only; no scientific import or main call.
    namespace=runpy.run_path(str(M/(PREFIX+'.py')),run_name='launch_gate_only')
    namespace['verify_gate'](GATE,MANIFEST,manifest)
    assert gate['command']==[sys.executable,str(ROOT/'scripts/s11c_guarded_run.py'),
        '--pool','s11c-near-unity','--memory-gib','4','--tasks-max','32',
        '--log-directory',str(RUN/'resource-guard'),'--',sys.executable,
        str(M/'S11c_d_end_normalization_run.py'),'--parallel-prerequisite-read',
        '--run-root',str(RUN),'--stage','first_order_receiving_regular_continue2','--',sys.executable,'-u',
        str(M/(PREFIX+'.py')),'--out',str(RUN/'complete'),'--inputs',str(MANIFEST),'--gate',str(GATE)]
    return gate,manifest

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
