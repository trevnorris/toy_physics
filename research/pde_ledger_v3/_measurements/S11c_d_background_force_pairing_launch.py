#!/usr/bin/env python3
"""One standing-authorized background force-pairing instrument; hook first, no retry."""
from datetime import datetime, timezone
import hashlib,json,os,select,shutil,subprocess,sys,time,traceback,runpy
from pathlib import Path
ROOT=Path('/var/projects/toy_physics')
M=ROOT/'research/pde_ledger_v3/_measurements'
PREFIX='S11c_d_background_force_pairing'
GATE=M/(PREFIX+'_gate.json')
MANIFEST=M/(PREFIX+'_inputs.json')
RUN=Path('/var/projects/toy_physics/_scratch/s11c/s11c-defect-packet-action-20261003/background-force-pairing-01')
HOOK=ROOT/'scripts/codex_job_watch.py'
THREAD='01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE='The user-authorized bounded BACKGROUND FORCE-PAIRING instrument has completed or reported an error.\nRoot: /var/projects/toy_physics/_scratch/s11c/s11c-defect-packet-action-20261003/background-force-pairing-01.\nInspect actual guard/supervisor/coordinator/no-deadline resources, strict stderr/stdout-checks, every source/copy/posthash, full registry/face/zero-unit/order/area/rate/pullback/control operands and returns. Exit is not acceptance; source/JSON/hash inspection only outside containment, no automatic validator/restoration/retry.\nThis first finite screen compares unbound background force representations and source-order ledgers only. b traction measure and total R0_c remain UNSUPPLIED; no support law or centre closure follows. Exact versus stored normals, saved versus exact velocity normalization, raw/retained area and sigma2 products remain separate. Legacy logical residual is literal only; residual zero face labels are not used, zero-slot dimensions are vacuous. Three controls are formal coefficient sensitivity, not nonzero physical power. No original producer/profile/grade/current/receiving/integral replay.\nFull selected incremental centre drive, physical nonuniform work, missing20/02,end-current and thresholdremainder remain open. Reflection is survival. No field/power integral,leakage/current/loss/inverse/sweep/drain/calibration. Mixed recovery PARKED; no general verifier expansion. Claude-only standing user authority covers faithful bounded physics work after inspection without routine questions.\nSame realomega3 restbulkLAB_HELD/RHO4,savededges,cs=sqrt6/2,arbitrarydoublet; this instrument leaves native coefficients unbound. Shared guard4GiBnative/cgroup within16GiBpool,zero swap,oneCPU/thread,32tasks,4GiBreserve,desktoppriority,no deadlines. HookFIRST01a0e01b-ef84-7192-817f-584cda5d339b. No scheduler/fallback/automaticretry/modelpolling/recurringtask. Preserve allhistory/Lean/S11_lean/sharedguard/protectedbuilder;scratchnevercommitted.\n'

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def save(path,value):
    with path.open('x') as f:
        json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
def verify():
    gate,manifest=read(GATE),read(MANIFEST)
    assert gate['status']=='READY_FOR_ONE_BACKGROUND_FORCE_PAIRING_INSTRUMENT'
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
        '--run-root',str(RUN),'--stage','background_force_pairing','--',sys.executable,'-u',
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
