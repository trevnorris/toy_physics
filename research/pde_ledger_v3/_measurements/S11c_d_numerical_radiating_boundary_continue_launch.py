#!/usr/bin/env python3
"""Single hook-first guarded launch, no scientific imports, retries or polling."""
from datetime import datetime,timezone
import hashlib,json,os,select,shutil,subprocess,sys,time,traceback
from pathlib import Path
ROOT=Path('/var/projects/toy_physics');M=ROOT/'research/pde_ledger_v3/_measurements'
PREFIX='S11c_d_numerical_radiating_boundary_continue'
GATE=M/(PREFIX+'_gate.json');MANIFEST=M/(PREFIX+'_inputs.json')
RUN=ROOT/'_scratch/s11c/s11c-d-numerical-radiating-20260929/boundary-continuation'
HOOK=ROOT/'scripts/codex_job_watch.py';THREAD='01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE='''The user-authorized no-time-limit saved-return numerical radiating boundary continuation has completed or reported an error.
Root: /var/projects/toy_physics/_scratch/s11c/s11c-d-numerical-radiating-20260929/boundary-continuation.
Inspect actual coordinator/launch, task-local guard timeAuthorization/duration-validation/effective limits/resource samples/outcomes/logs, existing supervisor numerical_radiating_boundary_continue invocation/active, strict scientific stderr, stdout/checks identity, checks/failure/incomplete operation, all source and prior posthashes, source snapshots/copy index, journal records and actual candidate/path/reference/sheet checks. Exit code is not scientific acceptance. Only JSON/source/hash/opaque-byte inspection outside containment; scientific restoration requires guard.
User said: Stop it with the time limits. We're losing work because of that crap. Seriously. Yes run a continuation. No wall, native or inactivity deadline is imposed. Verify actual RuntimeMaxUSec=infinity and Restart=no. Task-local copy preserves2GiB native/cgroup,zero swap,one CPU,32tasks,one thread,host4GiB reserve,locking and overlap refusal around the existing supervisor; shared guard unchanged. No scheduler changes, unguarded fallback or automatic scientific retry. The former900/840s limits and progress-stall timers do not apply to this continuation.
Restore all5001 complete old operations without calling their functions, with exact argument joins. Copy all36 old complete files byte-identically; open their ZIP once. Reconstitute only ephemeral Table callables from saved coefficients and the explicitly recorded unsaved runtime setup, not producer/root/Newton/derivative replay. First new operation is LEFT/omega-3/a-1/candidate-17/height-0.025/across/point-0035-joint-sheet, joined to the saved pending input. New immutable raw blobs use SQLite durable transactions; index distinguishes RESTORED_PRIOR_COMPLETE_RETURN and COMPLETE. No printed symbolic summary comparison. An argument/integrity mismatch is not physics; save evidence and stop. Both storage and restoration synthetic tests were local pure-tooling tests, not scientific acceptance.
Original boundary-01 remains preserved:840.270613worker seconds,251.9MiB peak,5001 complete operations,17/18 LEFT target candidates;RIGHT unvisited,zero finite solves. The prior timeout occurred in repeated ZIP-index loading; whole-runtime cost was not profiled. All71 old input hashes intact. No paired method clearance: round2 Claude CLEAR with five minor corrections, Grok no final report; user's explicit no-further-Grok execution authority remains. No new external review or exactomega1 restart.
If the actual boundary gate passes, continue the already-approved numerical pilot through current/face checks, endpoint integration, finite assembly and required controls, without repetitive stage permission. Preserve completed returns; do not relaunch completed calculations. A substantive central boundary/method failure stops neighbors and is reported without silently changing physics. Scope remainsLAB_HELD/RHO4_CONSTANT,saved tangents,omega3 then2.3/4,max30 finite solves. This boundary continuation alone supplies no face/current premise,finite deficit or physical loss. Total loss,calibration,Green/FORM/A11/A12 stay open.
Keep Lean/S11_lean,shared guard,protected builder suffixf01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2 and all accepted/failure/incident/review history. Canonical records outside scratch; scratch never committed. No model polling or recurring task. Hook targets01a0e01b-ef84-7192-817f-584cda5d339b.
'''
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def save(p,d):
 with p.open('x') as f:json.dump(d,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
def verify():
 g=read(GATE);m=read(MANIFEST)
 assert g['status']=='READY_FOR_USER_AUTHORIZED_NO_DEADLINE_BOUNDARY_CONTINUATION'
 assert g['independentMethodClearance'] is False
 assert g['proceedAuthority']=='EXPLICIT_USER_NO_TIME_LIMIT_SAVED_RETURN_CONTINUATION'
 assert g['launcherSha256']==sha(__file__) and g['manifestSha256']==sha(MANIFEST)
 assert g['completionMessageSha256']==hashlib.sha256(MESSAGE.encode()).hexdigest()
 for path,h in g['sourcePins'].items():assert sha(path)==h,path
 for item in m['packets'].values():assert sha(item['path'])==item['sha256'],item['path']
 assert str(RUN)==m['runRoot'] and str(RUN/'complete')==m['resultDirectory']
 assert g['command'][1]==str(M/(PREFIX+'_guard.py'))
 for rel,item in m['priorComplete']['files'].items():assert sha(Path(m['priorComplete']['path'])/rel)==item['sha256'],rel
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
