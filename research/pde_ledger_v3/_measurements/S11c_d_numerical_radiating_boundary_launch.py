#!/usr/bin/env python3
"""Single hook-first guarded launch, no scientific imports, retries or polling."""
from datetime import datetime,timezone
import hashlib,json,os,select,shutil,subprocess,sys,time,traceback
from pathlib import Path
ROOT=Path('/var/projects/toy_physics');M=ROOT/'research/pde_ledger_v3/_measurements'
PREFIX='S11c_d_numerical_radiating_boundary'
GATE=M/(PREFIX+'_gate.json');MANIFEST=M/(PREFIX+'_inputs.json')
RUN=ROOT/'_scratch/s11c/s11c-d-numerical-radiating-20260929/boundary-01'
HOOK=ROOT/'scripts/codex_job_watch.py';THREAD='01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE='''The user-authorized numerical radiating pilot central candidate-continuation gate has completed or reported an error.
Root: /var/projects/toy_physics/_scratch/s11c/s11c-d-numerical-radiating-20260929/boundary-01.
Inspect actual coordinator/launch, guard invocation/effective limits/limit validation/resource samples/child outcome/outcome/logs, supervisor numerical_radiating_boundary invocation/active, strict worker stderr, stdout/checks identity, checks/failure and incomplete operation. Inspect complete journal input/return receipts, actual point/candidate/refinement/reference/sheet JSON, source snapshots and posthashes. Exit code is not scientific acceptance. Scientific payload restoration only under containment; lightweight JSON/hash/source inspection is allowed. Preserve failed and completed work; no scientific retry.
The user removed the automatic two-day authoring cap and will decide the time budget. They expressly directed no further Grok pass, minor corrections and execution. Round2 Claude literally CLEAR with five minor pre-implementation corrections; Grok delivered no report. No paired independent clearance is claimed. Reviewed draft and raw/log-recovered reports remain at ba15c87b. No external review or exact omega1 restart.
This worker is the pilot's first required boundary gate: original full end pencils and all18 accepted candidates, contrast continuation then upper-half frequency paths to omega3, original-pencil/nullity/denominator/projector/refinement/independent joint sheet checks. Restores four pinned accepted packets; no producer or accepted-root replay. Literal contrast-independent sources may reuse completed modes/paths with their actual source joins. This is not a field solve, current/face premise, final modal map or deficit. All18 failures remain unresolved, never absent. Central method/control failure stops neighbors.
If the gate passes, the standing user goal authorizes continued implementation of the remaining numerical current/face checks, endpoint integration, finite assembly and controls, reusing completed returns. Do not ask repetitive science-stage approval. Actual readiness and ordinary shared guard/supervisor/native limits still apply. If a substantive boundary gate fails, report the exact evidence and stop; do not silently alter physics or expand the boundary method.
One900s outer/840s native,2GiB/zero swap/one CPU/nice15/32tasks/one thread via unchanged shared guard around existing supervisor. No inherited exact-track exceptions, overlap, fallback, scheduler changes or automatic retry. User goal remains LAB_HELD/RHO4_CONSTANT, saved tangents, frequencies3 then2.3/4, at most30 finite solves with uniform/scaling/resolution/domain/regulator/negative-deficit controls. Finite-model deficit is not automatically physical leakage. Total loss,analog-light calibration,Green/FORM/A11/A12 remain open. Preserve all history, Lean/S11_lean, shared guard and protected suffix f01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2. Canonical progress outside scratch; scratch never committed. No model polling or recurring task.
'''
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def save(p,d):
 with p.open('x') as f:json.dump(d,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
def verify():
 g=read(GATE);m=read(MANIFEST)
 assert g['status']=='READY_FOR_USER_AUTHORIZED_NUMERICAL_BOUNDARY_GATE'
 assert g['independentMethodClearance'] is False
 assert g['proceedAuthority']=='EXPLICIT_USER_NO_FURTHER_GROK_MINOR_CORRECTIONS_AND_RUN'
 assert g['launcherSha256']==sha(__file__) and g['manifestSha256']==sha(MANIFEST)
 assert g['completionMessageSha256']==hashlib.sha256(MESSAGE.encode()).hexdigest()
 for path,h in g['sourcePins'].items():assert sha(path)==h,path
 for item in m['packets'].values():assert sha(item['path'])==item['sha256'],item['path']
 assert str(RUN)==m['runRoot'] and str(RUN/'complete')==m['resultDirectory']
 assert g['command'][1]==str(ROOT/'scripts/s11c_guarded_run.py')
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
