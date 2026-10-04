#!/usr/bin/env python3
"""One standing-authorized source-projection and real-axis receiving prerequisite; hook first, no retry."""
from datetime import datetime, timezone
import hashlib,json,os,select,shutil,subprocess,sys,time,traceback,runpy
from pathlib import Path
ROOT=Path('/var/projects/toy_physics')
M=ROOT/'research/pde_ledger_v3/_measurements'
PREFIX='S11c_d_first_order_receiving_regular'
GATE=M/(PREFIX+'_gate.json')
MANIFEST=M/(PREFIX+'_inputs.json')
RUN=Path('/var/projects/toy_physics/_scratch/s11c/s11c-defect-packet-action-20261003/first-order-receiving-regular-01')
HOOK=ROOT/'scripts/codex_job_watch.py'
THREAD='01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE='The user-authorized bounded SOURCE-PROJECTION / CROSS-FLUX / REAL-AXIS RECEIVING prerequisite has completed or reported an error.\nRoot: /var/projects/toy_physics/_scratch/s11c/s11c-defect-packet-action-20261003/first-order-receiving-regular-01.\nInspect actual coordinator/guard/supervisor/no-deadline resources, strict stderr/stdout-checks identity, every copied/source posthash, inherited proof arguments and all new full-source projections, native affine pressure joins, moving-projector distributions, opposite-leg current, original denominator, exact determinant/Bezout/Sturm/adjugate/growth/control operands and returns. Source/JSON/hash inspection only outside containment; no scientific restoration, automatic validator or retry. Preserve any pole, nonzero residue/cross current or unknown decision. Exit is not acceptance.\nMethod literal Claude CLEAR is source-only and partial; concrete build and gate require separate inspection. Full five rows precede R3 indices2/3/4. Native UB full three-row source zero is not zero leakage; nonlocalized transverse step stays matched. Q multiplies the whole step before the receiving cusp. Native pressure/chemical normalization enters once and both lab normals remain. C3 root exclusion covers both outgoing rays at fixed physical inputs; whole five-field operator still has transverse poles. Conditional weightedL2 uniqueness/L1 decay includes E3 growth, not decay rate, general distributional inverse or far-bulk flux. All old source/profile/current/moment/receiving constructors/proofs remain saved without replay.\nPhysical nonuniform work identity including exterior acoustic flux, memory dissipation, mass-rate/chemical and LAB_HELD work is REQUIRED next before quadratic leakage. No field/transform/power integral, leakage/current value, sweep, drain or calibration. Mixed numerical recovery PARKED; no general verifier expansion. Claude-only standing authority permits faithful bounded next physics work after actual inspection without routine permission. Same realomega3 restbulkLAB_HELD/RHO4,savededges,cs=sqrt6/2,arbitrary incident doublet. Guard4GiBnative/cgroup16GiBpool,zero swap,oneCPU/thread,32tasks,4GiBreserve,desktoppriority,no deadlines,RuntimeMaxUSec=infinity/Restart=no;hookFIRST01a0e01b-ef84-7192-817f-584cda5d339b. No fallback/schedulerchange/automaticretry/modelpolling/recurringtask. Preserve allhistory/Lean/S11_lean/sharedguard/protectedbuilder;scratch nevercommitted.\n'

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def save(path,value):
    with path.open('x') as f:
        json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
def verify():
    gate,manifest=read(GATE),read(MANIFEST)
    assert gate['status']=='READY_FOR_ONE_FIRST_ORDER_RECEIVING_REGULAR_INSTRUMENT'
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
        '--run-root',str(RUN),'--stage','first_order_receiving_regular','--',sys.executable,'-u',
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
