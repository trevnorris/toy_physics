#!/usr/bin/env python3
"""One standing-authorized reference-response instrument; hook first, no retry."""
from datetime import datetime, timezone
import hashlib,json,os,select,shutil,subprocess,sys,time,traceback,runpy
from pathlib import Path
ROOT=Path('/var/projects/toy_physics')
M=ROOT/'research/pde_ledger_v3/_measurements'
PREFIX='S11c_d_defect_end_uniform'
GATE=M/(PREFIX+'_gate.json')
MANIFEST=M/(PREFIX+'_inputs.json')
RUN=Path('/var/projects/toy_physics/_scratch/s11c/s11c-defect-end-uniform-20261003/diagnostic-01')
HOOK=ROOT/'scripts/codex_job_watch.py'
THREAD='01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE="The standing-user-authorized bounded weak-end/selected-uniform comparison has completed or reported an error.\nRoot: /var/projects/toy_physics/_scratch/s11c/s11c-defect-end-uniform-20261003/diagnostic-01.\nInspect actual coordinator/launch/shared-guard/supervisor/resources/no-deadline enforcement/outcomes, strict stderr and stdout/checks bytes, failure/incomplete operation, every restored/new source/chart/grade/domain/limit/control argument and return, evidence chain, copies/snapshots and all posthashes. Exit code is not acceptance. Source/JSON/hash/opaque inspection only outside containment; no scientific restoration or automatic retry/validator.\nMethod pair literally CLEAR at e41557f0; inspect actual concrete build assessment and gate separately. Reuse the200weakendcells and original uniform source/lift/restriction/limits; no producers, selected-mode/current/limit/inventory replay. New source cyclic covariance joins weak profile axis1 to native uniform axis3; inspect native fields/rows/units/epsilon/Fourier/grade/physical inputs and actual source identity. Do not assume a permutation by labels. Failure is SOURCE_MAP_UNRESOLVED, not a physics or loss result.\nSelected A/R are dependent via inherited raw I_old; its zero return is on-wave only. Compare independent retained grades with raw-source/finite-origin joins and exact excluded remainder before attribution. FullDelta may differ, retained mismatch and finite-truncation difference are separate. Save actual domains, closed grazingR0 plus4inheritedpath joins, per-cell lift participation and changed-minus-baseline controls. Pressure/direct entries annihilated by the lift remain untested by a selected pass. No physical current or finite inverse follows.\nRealomega3 strictrestbulk LAB_HELD/RHO4,savededges,cs[1,2],selected equation correspondence only. No new roots, frequency derivative, plane-wave scattering, evaluated integral, finite solve, production overwrite, defect sweep or loss. Drain/calibration/oldfinite-route debts remain open. Preserve old late summary failure and all completed weak-end results without replay.\nStanding user requests continued clear bounded work without routine permission. Science only shared no-deadline pooled guard around existing supervisor:4GiBnative/cgroup in16GiBpool,zero swap,oneCPU/thread,tasks32,4GiBhostreserve,desktoppriority,RuntimeMaxUSec=infinity/Restart=no. Hook first to01a0e01b-ef84-7192-817f-584cda5d339b. No fallback/schedulerchanges/automatic science retry/model polling/recurringtask. Preserve Lean/S11_lean/sharedguard/protectedbuilder and all history. Scratch never committed.\n\nLiteral corrected build reviews remain Claude NEEDS REVISION/Grok CLEAR at3803e3fb. Claude identified omitted generated gamma units in the static schema. The tested local metadata repair reads the SAME four already-restored end-source/pairing unit registries, joins original object and frozen producer receipts and exact native symbols, and rejects missing/conflicting/nonrational units. The dimension predicate and comparison equations/assignments remain unchanged. No new dimension inference, source calculation or independent CLEAR. Gate records independentBuildClearance=false and standing local-tooling execution authority. Inspect actual gamma registry operands and complete native unit results. Grok's incidental zero-first-cell statement is refuted by the frozen JSON and preserved literally.\n"

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def save(path,value):
    with path.open('x') as f:
        json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
def verify():
    gate,manifest=read(GATE),read(MANIFEST)
    assert gate['status']=='READY_FOR_ONE_END_UNIFORM_COMPARISON'
    assert gate['independentBuildClearance'] is False and gate['localToolingExecutionAuthority'] is True and gate['pooledExecution'] is True
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
        '--run-root',str(RUN),'--stage','defect_end_uniform','--',sys.executable,'-u',
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
