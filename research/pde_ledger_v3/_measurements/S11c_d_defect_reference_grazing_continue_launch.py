#!/usr/bin/env python3
"""One standing-authorized saved-evidence reference continuation; hook first, no retry."""
from datetime import datetime, timezone
import hashlib,json,os,select,shutil,subprocess,sys,time,traceback,runpy
from pathlib import Path
ROOT=Path('/var/projects/toy_physics')
M=ROOT/'research/pde_ledger_v3/_measurements'
PREFIX='S11c_d_defect_reference_grazing_continue'
GATE=M/(PREFIX+'_gate.json')
MANIFEST=M/(PREFIX+'_inputs.json')
RUN=Path('/var/projects/toy_physics/_scratch/s11c/s11c-defect-reference-grazing-20261002/continuation-01')
HOOK=ROOT/'scripts/codex_job_watch.py'
THREAD='01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE='The standing-user-authorized saved-evidence reference-grazing continuation has completed or reported an error.\nRoot: /var/projects/toy_physics/_scratch/s11c/s11c-defect-reference-grazing-20261002/continuation-01.\nInspect actual coordinator/launch/guard/supervisor logs, no-deadline resources/outcomes, strict scientific stderr, stdout/checks bytes, failure/incomplete operation, every prior-copy/hash/restored operand-return and context join, exact constant certificates and all new bounds/limits/units/controls. Exit code is not acceptance. JSON/source/hash/opaque inspection only outside containment; no scientific restoration or automatic retry/validator.\nOriginal failed result8da2de5f is preserved:454files/1251934bytes,129literal zero returns,427artifacts,15.568s,peak78917632bytes,zero swap/events,59posthashes/38snapshots/21savedcopies intact. Failure was exact structural zero refusal of saved Add(Integer(-9),Integer(9)) at new-tail-root-gap-reconstruction. The new exact rational constructor decision accepts only Integer/Rational/Add/Mul/bounded integer Pow with no symbols/floats/unknowns; requires exact zero, persists full operands and arithmetic steps. Cause of old unevaluated Add remains unestablished. Do not label prior top-level stage complete.\nVerify all454 prior files copied byte-identically,129published zero returns and operands restored without their functions, failed actual residual/operands preserved, and new context uses actual saved coefficients, symbols and assumptions. Only small scalar/argument/lambda context is reconstructed. The70-statement original unfinished tail starts at new-tail-polynomial and is AST-identical; no native source/closure/trace/contact/PV/sheet prefix replay. Four controls must now run and respond. Bounds and analytic limits remain assessed mathematics, not machine theorems. No integral, fullsourceconsumer, matrix, production overwrite, defectsweep or loss.\nOriginal method/build pairs CLEAR remain source-only for original worker. This saved-evidence continuation is a tested local representation/checkpoint repair, not fresh independent build clearance. Gate records continuationIndependentBuildClearance=false and standing local-tooling authority. Preserve every literal report and failed source. Scientific work only under unchanged shared no-deadline pooled guard around existing supervisor:4GiB native/cgroup within16GiBpool,zero swap,oneCPU/thread,tasks32,4GiBhostreserve,desktoppriority,RuntimeMaxUSec=infinity/Restart=no. No scheduler changes/fallback. Hook first to01a0e01b-ef84-7192-817f-584cda5d339b. No model polling/recurring task.\nContinue clear authorized bounded next steps after actual inspection without routine permission. Preserve Lean/S11_lean,sharedguard,protected builder and all historical artifacts. Scratch never committed; no immediate defect sweep. Limits restbulk LAB_HELD/RHO4_CONSTANT,omega3,saved edges,cs[1,2],k/l[-3,3]; old cutoffs4/6,fulloperator,finiteinverse,drain,calibration,kappa0/beta0 remain excluded.\n'

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def save(path,value):
    with path.open('x') as f:
        json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
def verify():
    gate,manifest=read(GATE),read(MANIFEST)
    assert gate['status']=='READY_FOR_ONE_SAVED_REFERENCE_CONTINUATION'
    assert gate['continuationIndependentBuildClearance'] is False and gate['localToolingRepairAuthorized'] is True and gate['pooledExecution'] is True
    assert gate['launcher']==manifest['launcher']==str(Path(__file__).resolve())
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
        '--run-root',str(RUN),'--stage','defect_reference_grazing_continue','--',sys.executable,'-u',
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
