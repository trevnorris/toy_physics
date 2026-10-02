#!/usr/bin/env python3
"""One standing-authorized guarded raw increment; hook first, no retry."""
from datetime import datetime, timezone
import hashlib,json,os,select,shutil,subprocess,sys,time,traceback,runpy
from pathlib import Path
ROOT=Path('/var/projects/toy_physics')
M=ROOT/'research/pde_ledger_v3/_measurements'
PREFIX='S11c_d_defect_raw_increment_continue'
GATE=M/(PREFIX+'_gate.json')
MANIFEST=M/(PREFIX+'_inputs.json')
RUN=Path('/var/projects/toy_physics/_scratch/s11c/s11c-defect-raw-increment-20261001/continuation-01')
HOOK=ROOT/'scripts/codex_job_watch.py'
THREAD='01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE="The user-renewed saved-evidence raw-increment CONTROL continuation has completed or reported an error.\nRoot: /var/projects/toy_physics/_scratch/s11c/s11c-defect-raw-increment-20261001/continuation-01.\nInspect coordinator/launch/stage-start, shared-guard duration/enforced resources/outcomes/logs, supervisor defect_raw_increment_continue invocation/active, strict stderr, stdout/checks identity, checks/failure/incomplete operation, all journal input/return/control records, prior-copy index and source/copy posthashes. Process exit is not scientific acceptance. Source/JSON/hash/opaque-byte inspection only outside containment; no scientific restoration or additional validator/retry.\nUser said Continue after preserved formatter failure fcc3ee32. Preserve all306 prior files/10406596bytes and four completed nested operations. Restore lower-boundary, physical-factorization, wrong-height and omitted-slope returns without calling functions. Restore saved row/closure/source evidence and join actual argument/symbol/grade/physical-input context. No boundary, factorization, source/grade/census/row construction or completed-function replay. Only small exact callable/argument context and the original unfinished control point declarations are reconstructed.\nThe only changed original control-tail AST is depth-map evidence serialization: symbolic keys become strings with collision refusal; all15 mathematical assignment ASTs and acceptance rules are unchanged.22stdlib tests pass; no science was restored during preparation. Literal build reports remain Claude CLEAR FOR THIS BOUNDED RAW-INCREMENT BUILD/Grok NEEDS REVISION, not a fresh independent CLEAR. Original failed evidence, worker and review packets unchanged. Gate records renewed continuation and local-tooling authority honestly.\nInspect actual height/slope dispersion joins; wrong-height/slope/sheet nonzero responses; THETA/E_W actual lower-row omission/doubling; lower-jet sign and final constant-end checks. All must have actual saved evidence. Remaining nonzero/unknown or integrity failure is preserved and stops, not softened. Eight response controls are sensitivity checks, not independent correctness. Nongrazing both-face raw operator only; exact-match closed limits/local bound unresolved. No evaluated integral, finite matrix, production change, defect sweep or loss claim.\nNo deadline. Shared guard around existing supervisor;4GiB native/cgroup in16GiB pool,zero swap,one assignedCPU/thread,tasks32,4GiB host reserve,desktop priority,RuntimeMaxUSec=infinity/Restart=no. No scheduler changes, fallback, automatic retry, model polling or recurring task. Hook to01a0e01b-ef84-7192-817f-584cda5d339b. Preserve Lean/S11_lean,shared guard,protected suffixf01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2 and all result/failure/review/incident history. Scratch never committed. Stop with this continuation's actual evidence; no immediate defect sweep.\n"

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def save(path,value):
    with path.open('x') as f:
        json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
def verify():
    gate,manifest=read(GATE),read(MANIFEST)
    assert gate['status']=='READY_FOR_ONE_SAVED_RAW_CONTROL_CONTINUATION'
    assert gate['independentBuildClearance'] is False and gate['pooledExecution'] is True
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
        '--run-root',str(RUN),'--stage','defect_raw_increment_continue','--',sys.executable,'-u',
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
