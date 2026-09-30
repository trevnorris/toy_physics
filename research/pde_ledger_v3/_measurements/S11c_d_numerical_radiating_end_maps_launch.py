#!/usr/bin/env python3
"""Single hook-first guarded launch, no scientific imports, retries or polling."""
from datetime import datetime,timezone
import hashlib,json,os,select,shutil,subprocess,sys,time,traceback
from pathlib import Path
ROOT=Path('/var/projects/toy_physics');M=ROOT/'research/pde_ledger_v3/_measurements'
PREFIX='S11c_d_numerical_radiating_end_maps'
GATE=M/(PREFIX+'_gate.json');MANIFEST=M/(PREFIX+'_inputs.json')
RUN=ROOT/'_scratch/s11c/s11c-d-numerical-radiating-20260929/end-maps-01'
HOOK=ROOT/'scripts/codex_job_watch.py';THREAD='01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE='''The user-authorized numerical radiating current/face/end-map job has completed or reported an error.
Root: /var/projects/toy_physics/_scratch/s11c/s11c-d-numerical-radiating-20260929/end-maps-01.
Inspect actual coordinator/launch/stage-start, guard duration/effective-limit/resource/outcome records and logs, supervisor numerical_radiating_end_maps invocation/active, strict scientific stderr, stdout/checks byte identity, checks/failure/failed-check evidence, all journal operands/returns, result-artifact receipts, source snapshots and posthashes. Exit status is not science acceptance. JSON/source/hash/opaque bytes are lightweight; scientific restoration stays inside containment. Preserve every completed and failed operation; no automatic scientific retry.
The preceding central boundary continuation is accepted only for candidate continuation: both ends, all18 candidates, four contrasts, 30131 complete operations, 48.8min, peak555.7MiB, all inputs intact. This worker restores its saved final states and original accepted pairing/frequency/unit/endpoint packets. It binds omega3 and materials before evaluating physical current and saved face rows; uses epsilon-only homogeneity plus rescaling rather than Poly extraction. Inspect actual pencil/physical-q joins, both harmonic face legs, transverse drives/current normalization, native velocity omission, off-diagonal/sign controls, full cross-current forms and five-outgoing/two-incoming trace rank/residuals. Identical source/state joins can reuse a completed map. There is no root/path/producer replay or field/integral solve.
Standing user policy removes all wall/native/CPU-time/inactivity limits. Canonical shared scripts/s11c_guarded_run.py now enforces RuntimeMaxUSec=infinity/Restart=no with2GiB native/cgroup,zero swap,one CPU,32tasks,one thread,4GiB host reserve,locking and overlap refusal around the existing supervisor. Prior pinned guard snapshots are preserved. No scheduler changes,fallback or retry.
User explicitly directed no further Grok pass, minor corrections and run. Literal history remains round2 Claude CLEAR with five minor corrections; Grok no final report. No paired clearance or new review is claimed. If actual central end-map/current/face checks pass, continue the already approved numerical pilot through endpoint integration, finite assembly and uniform/scaling/refinement/domain/regulator/sign controls without repetitive stage permission. A substantive central method failure stops neighbors; do not silently change physics. Scope remainsLAB_HELD/RHO4_CONSTANT,saved tangents,omega3 then2.3/4,max30 finite solves. No physical leakage,loss magnitude,calibration,Green/FORM/A11/A12 acceptance follows from an end map. Preserve all historical results/failures/review debt/incident records,Lean/S11_lean and protected builder suffixf01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2. Scratch never committed. No model polling or recurring task. Hook targets01a0e01b-ef84-7192-817f-584cda5d339b.
'''
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def save(p,d):
 with p.open('x') as f:json.dump(d,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
def verify():
 g=read(GATE);m=read(MANIFEST)
 assert g['status']=='READY_FOR_AUTHORIZED_NUMERICAL_END_MAPS'
 assert g['independentMethodClearance'] is False
 assert g['proceedAuthority']=='EXPLICIT_USER_NO_FURTHER_GROK_MINOR_CORRECTIONS_AND_RUN'
 assert g['launcherSha256']==sha(__file__) and g['manifestSha256']==sha(MANIFEST)
 assert g['completionMessageSha256']==hashlib.sha256(MESSAGE.encode()).hexdigest()
 for path,h in g['sourcePins'].items():assert sha(path)==h,path
 for item in m['packets'].values():assert sha(item['path'])==item['sha256'],item['path']
 assert str(RUN)==m['runRoot'] and str(RUN/'complete')==m['resultDirectory']
 assert g['command'][1]==str(ROOT/'scripts/s11c_guarded_run.py')
 assert sha(m['boundary']['database'])==m['boundary']['sha256']
 return g,m

def coordinate(fd):
 try:
  ready,_,_=select.select([fd],[],[],30)
  if not ready or os.read(fd,1)!=b'1':raise RuntimeError('Hook not armed; no scientific launch')
  os.close(fd);g,m=verify();start=time.monotonic()
  save(RUN/'stage-start.json',{'utc':datetime.now(timezone.utc).isoformat(),'gateSha256':sha(GATE),'command':g['command']})
  with (RUN/'guard-launch.stdout').open('x') as out,(RUN/'guard-launch.stderr').open('x') as err:
   result=subprocess.run(g['command'],cwd=ROOT,stdin=subprocess.DEVNULL,stdout=out,stderr=err)
  save(RUN/'coordinator-outcome.json',{'exitCode':result.returncode,'wallSeconds':time.monotonic()-start,'scientificAcceptance':'PENDING_ACTUAL_INSPECTION'})
  return result.returncode
 except BaseException:
  save(RUN/'coordinator-failure.json',{'traceback':traceback.format_exc(),'automaticRetry':False});return 1

def launch():
 g,m=verify();assert shutil.which('codex'),'completion queue unavailable'
 RUN.mkdir(exist_ok=False)
 for source,h in {**g['sourcePins'],str(GATE):sha(GATE),str(MANIFEST):sha(MANIFEST)}.items():
  path=Path(source);target=RUN/'source'/path.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True)
  shutil.copyfile(path,target);assert sha(target)==h
 save(RUN/'source-snapshots.json',{p:sha(RUN/'source'/Path(p).relative_to(ROOT)) for p in g['sourcePins']})
 (RUN/'completion-message.txt').write_text(MESSAGE)
 r,w=os.pipe()
 with (RUN/'coordinator.stdout').open('x') as out,(RUN/'coordinator.stderr').open('x') as err:
  job=subprocess.Popen(['/usr/bin/python3',str(Path(__file__).resolve()),'--coordinator',str(r)],cwd=ROOT,stdin=subprocess.DEVNULL,stdout=out,stderr=err,pass_fds=(r,),start_new_session=True)
 os.close(r)
 try:
  with (RUN/'watcher.stdout').open('x') as out,(RUN/'watcher.stderr').open('x') as err:
   watcher=subprocess.Popen(['/usr/bin/python3',str(HOOK),'--pid',str(job.pid),'--thread',THREAD,'--directory',str(RUN/'completion-watcher'),'--message-file',str(RUN/'completion-message.txt')],cwd=ROOT,stdin=subprocess.DEVNULL,stdout=out,stderr=err,start_new_session=True)
  save(RUN/'launch.json',{'jobPid':job.pid,'watcherPid':watcher.pid,'gateSha256':sha(GATE),'command':g['command'],'thread':THREAD,'automaticRetry':False})
  for _ in range(50):
   path=RUN/'completion-watcher/state.json'
   if path.exists() and read(path)['status']=='waiting':
    os.write(w,b'1');print(json.dumps({'runDirectory':str(RUN),'jobPid':job.pid,'watcherPid':watcher.pid,'hookStatus':'waiting'}));return
   if job.poll() is not None or watcher.poll() is not None:raise RuntimeError('Coordinator/watcher exited before arming')
   time.sleep(.1)
  raise RuntimeError('Hook did not arm')
 finally:os.close(w)
if __name__=='__main__':
 if len(sys.argv)==3 and sys.argv[1]=='--coordinator':sys.exit(coordinate(int(sys.argv[2])))
 if len(sys.argv)!=1:raise ValueError('Unexpected launch arguments')
 launch()
