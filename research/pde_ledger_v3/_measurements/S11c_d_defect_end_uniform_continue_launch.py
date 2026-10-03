#!/usr/bin/env python3
"""One standing-authorized reference-response instrument; hook first, no retry."""
from datetime import datetime, timezone
import hashlib,json,os,select,shutil,subprocess,sys,time,traceback,runpy
from pathlib import Path
ROOT=Path('/var/projects/toy_physics')
M=ROOT/'research/pde_ledger_v3/_measurements'
PREFIX='S11c_d_defect_end_uniform_continue'
GATE=M/(PREFIX+'_gate.json')
MANIFEST=M/(PREFIX+'_inputs.json')
RUN=Path('/var/projects/toy_physics/_scratch/s11c/s11c-defect-end-uniform-20261003/continuation-01')
HOOK=ROOT/'scripts/codex_job_watch.py'
THREAD='01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE='The standing-user-authorized saved-evidence RIGHT end/uniform continuation has completed or reported an error.\nRoot: /var/projects/toy_physics/_scratch/s11c/s11c-defect-end-uniform-20261003/continuation-01.\nInspect actual coordinator/launch/shared-guard/supervisor commands, duration/resources/outcomes, strict stderr, stdout/checks identity, failure/incomplete operation, every prior-copy/source posthash, chain, restored operands/returns, typed-operand selectors and new RIGHT grazing/control evidence. Exit code is not scientific acceptance. Source/JSON/hash/opaque-byte inspection only outside containment; no scientific restoration, automatic retry or extra validator.\nOriginal failure4c33c471 remains:6885files/540133852bytes,6796chained artifacts,2281literal zeros,18restoredoldoperations,0newtopcompleteoperations. Both ends saved selectedA/R and all4retainedgrade agreement on declared wave/denominator domains. LEFT completed2grazing targets/4path joins/3responsive controls. RIGHT stopped BEFORE publishing native-minus-pairing join in JSON allocation MemoryError; its result remains unresolved until this actual continuation.4916.114s,peak4236333056bytes,zero swap/events,no timeout or cgroup kill. All old files are copied byte-identically; no LEFT/source/chart/unit/normalization/grade/selected comparison replay.\nThe local tooling repair streams JSON and refers large exact operands to full immutable raw blob receipts plus literal selectors. The unchanged exact_structure predicate decides equality; pointers or hashes alone do not. Full join operand receipts precede comparison and return; an identical completed join may reuse its return only for the same live operand objects. No printed summary comparison, numerical tolerance or equation change. Inspect actual pointers/joins; missing provenance or false equality is fatal. The original RIGHT minus tail starts at failed pairing join; remaining plus grazing and controls use AST-identical source. Only unsaved small matrix/argument/participation context is reconstructed from published cells/operands. Current build history remains Claude NEEDS REVISION/Grok CLEAR at3803e3fb and original gamma metadata repair; continuationIndependentBuildClearance=false, no new review or authorCLEAR.\nInspect both RIGHT closedR0 targets, all4savedradiating/evanescentpaths, changed-minus-baseline omit-cell/lift/sheet controls, original domain certificates and all200cell participation. A/R are dependent through rawI_old; selected zero-weighted pressure/direct entries stay untested. No source/limit/producer functions, roots, current, response integral, finite solve, production overwrite, defect sweep or loss.\nStanding user directs continued clear bounded work without routine permission. Same realomega3 strictrestbulk LAB_HELD/RHO4,savededges,cs[1,2],independentrectangle,selected equation correspondence only. Preserve every historical failure/report/result,Lean/S11_lean,sharedguard,protectedbuilder. No deadline; shared guard around existing supervisor,4GiBnative/cgroup within16GiBpool,zero swap,oneCPU/thread,tasks32,4GiBhostreserve,desktop priority,RuntimeMaxUSec=infinity/Restart=no. Hook first tosession01a0e01b-ef84-7192-817f-584cda5d339b. Scratch never committed. No scheduler/fallback/modelpolling/recurringtask. Drain/calibration/oldfinite-route debts stay open.\n'

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def save(path,value):
    with path.open('x') as f:
        json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
def verify():
    gate,manifest=read(GATE),read(MANIFEST)
    assert gate['status']=='READY_FOR_ONE_SAVED_END_UNIFORM_CONTINUATION'
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
        '--run-root',str(RUN),'--stage','defect_end_uniform_continue','--',sys.executable,'-u',
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
