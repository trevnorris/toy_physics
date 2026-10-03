#!/usr/bin/env python3
"""One standing-authorized packet adapter/tail preflight; hook first, no retry."""
from datetime import datetime, timezone
import hashlib,json,os,select,shutil,subprocess,sys,time,traceback,runpy
from pathlib import Path
ROOT=Path('/var/projects/toy_physics')
M=ROOT/'research/pde_ledger_v3/_measurements'
PREFIX='S11c_d_defect_packet_fourier'
GATE=M/(PREFIX+'_gate.json')
MANIFEST=M/(PREFIX+'_inputs.json')
RUN=Path('/var/projects/toy_physics/_scratch/s11c/s11c-defect-packet-action-20261003/fourier-01')
HOOK=ROOT/'scripts/codex_job_watch.py'
THREAD='01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE='The standing-user-authorized bounded packet FOURIER request bank has completed or reported an error.\nRoot: /var/projects/toy_physics/_scratch/s11c/s11c-defect-packet-action-20261003/fourier-01.\nInspect actual launch/coordinator/guard/supervisor commands and enforced no-deadline resources, strict stderr/stdout-checks identity, failed/incomplete request, all original saved/copy/source hashes and SQLite FULL-transaction operands/returns, every rule/node/panel/tail/radius/absolute comparison. Process exit is not scientific acceptance. JSON/source/hash/opaque-byte inspection only outside containment; no scientific restoration or automatic retry.\nThe preflight ccb35e8c is reused, not replayed: 544 selected addresses, 34 coefficient vectors, all native jets and original source proofs, K27/U75/T122. This instrument computes only190 declared transforms:13 actual source-jet families plus3 tests and3 test derivatives, two carriers, five signed physical transform arguments. Inspect exact coefficients, native derivatives before source multiplication, X/Y signs and normalization, inherited unit context and actual native arguments. Separate shifted-contour compositeGL24/48 at30digits and real-axis adaptiveG7/K15 at50digits; each request gets Gaussian-moment tail radius first. Exact Stieltjes/rule moments, all open points/weights/values and subdivisions are persisted. Constant fields have independent analytic Gaussian references. Absolute1e-12 miss stops, no higher-order retry. Adaptive error is empirical, not an analytic quadrature proof.\nNo response/H/J/D integral, complete packet action, full summand-unit proof, uniform transform-error bound, scattering/current/loss follows. Future evaluator still must compare every actual request and complete independent collision-aware action routes/local tails/physical profile and paired-PV joins. Reuse exact completed requests only with actual argument identity.\nStanding user/AGENTS authority covers continued clear bounded work and necessary established-Claude/Grok reviews without routine questions. Science only shared no-deadline pooled guard around existing supervisor,4GiBnative/cgroup in16GiBpool,zero swap,oneCPU/thread,tasks32,4GiBhostreserve,desktoppriority,RuntimeMaxUSec=infinity/Restart=no. Hook first tosession01a0e01b-ef84-7192-817f-584cda5d339b. No fallback/scheduler changes/automatic scientific retry/modelpolling/recurringtask. Preserve allhistory/Lean/S11_lean/sharedguard/protectedbuilder. Scratch nevercommitted. RestbulkLAB_HELD/RHO4,realomega3,savededges,one effectivecs=sqrt6/2. No sweep,drain/calibration or loss.\n'

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def save(path,value):
    with path.open('x') as f:
        json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
def verify():
    gate,manifest=read(GATE),read(MANIFEST)
    assert gate['status']=='READY_FOR_ONE_PACKET_FOURIER_BANK'
    assert gate['independentBuildClearance'] is True and gate['pooledExecution'] is True
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
        '--run-root',str(RUN),'--stage','defect_packet_fourier','--',sys.executable,'-u',
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
        relative = path.relative_to(ROOT) if path.is_relative_to(ROOT) else Path('external-runtime') / hashlib.sha256(str(path).encode()).hexdigest() / path.name
        target = RUN / 'source' / relative
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
