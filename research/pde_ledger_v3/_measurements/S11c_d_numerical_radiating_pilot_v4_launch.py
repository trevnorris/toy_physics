#!/usr/bin/env python3
"""Single hook-first guarded launch, no scientific imports, retries or polling."""
from datetime import datetime,timezone
import hashlib,json,os,select,shutil,subprocess,sys,time,traceback
from pathlib import Path
ROOT=Path('/var/projects/toy_physics');M=ROOT/'research/pde_ledger_v3/_measurements'
PREFIX='S11c_d_numerical_radiating_pilot_v4'
GATE=M/(PREFIX+'_gate.json');MANIFEST=M/(PREFIX+'_inputs.json')
RUN=ROOT/'_scratch/s11c/s11c-d-numerical-radiating-20260929/central-pilot-v4'
HOOK=ROOT/'scripts/codex_job_watch.py';THREAD='01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE='The user-authorized central omega3 numerical radiating pilot v4 has completed or reported an error.\nRoot: /var/projects/toy_physics/_scratch/s11c/s11c-d-numerical-radiating-20260929/central-pilot-v4.\nInspect actual coordinator/launch/stage-start; shared guard duration/effective limits/resource samples/outcomes/logs; supervisor numerical_radiating_pilot_v4 invocation/active; strict scientific stderr; stdout/checks identity; checks/failure/incomplete operation; all journal operands/returns and result/partial-array receipts; snapshots and source posthashes. Exit code is not scientific acceptance. Outside containment use only JSON/source/hash/opaque bytes; no scientific restoration.\nThis is the same approved LAB_HELD/RHO4_CONSTANT central omega3 pilot at saved tangents. Accepted end maps and candidate paths are reused; no producer/root/path/current/face replay. Restore all twelve completed v2 returns without calling their functions, copy its twenty-five files exactly, join actual inputs and resume unfinished integration/transformed-measure. V3 completed only eight restorations, then plain nested-dictionary equality attempted a Boolean conversion of a NumPy array. No new numerical work or finite solve. All 39 files/53,363,533 bytes, 69 input pins, 40 snapshots and 25 copied files were checked and remain in place; see pilot_array_disposition.json. Verify their hashes during completion inspection. V4 only imports and uses existing exact_structure for the rows/jets/source reuse comparison: container types/keys/lengths, array dtype/shape/exact bytes and scalar equality. No tolerance or scientific equation change. Six synthetic stdlib regression/reuse tests pass; actual runtime identity remains required. Both numerical libraries are unchanged. Earlier initialization and floating-panel failures stay preserved.\nInspect actual 375 source joins, 160 native occurrences, 80 endpoint bounds, physical-root joins, weighted endpoint values, independent middle-leg quadrature, responsive wrong-sheet/Jacobian controls and uniform Gaussian operator versus original REFERENCE pencil before accepting finite assembly. At most ten central solves retain full current cross terms and uniform/contrast-scaling/refinement/domain/regulator/negative-deficit controls. A substantive scientific central check failure stops neighbors; preserve evidence, no automatic scientific retry or silent method change. If central controls pass, the standing user goal covers omega2.3/4 after actual readiness, max30 finite solves; no repetitive stage permission.\nNo wall/native/CPU-time/inactivity deadline: shared scripts/s11c_guarded_run.py around existing normalization supervisor, RuntimeMaxUSec=infinity/Restart=no,2GiB native/cgroup,zero swap,one CPU,32tasks,one thread,4GiB host reserve,locking/overlap refusal. No scheduler change or fallback. User expressly authorized minor corrections and execution with no further Grok pass. Literal method history remains round2 Claude CLEAR with5minor corrections, Grok no final report; no paired clearance. Finite-model deficit is not automatically physical leakage. Retained-order loss interpretation, analog-light calibration, Green/FORM/A11/A12 remain open. Exact omega1 stays parked.\nPreserve every accepted/failed result, incident/review history, Lean/S11_lean and protected builder suffix f01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2. Canonical records outside scratch; scratch never committed. Hook targets session01a0e01b-ef84-7192-817f-584cda5d339b. No model polling or recurring task.\n'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def save(p,d):
 with p.open('x') as f:json.dump(d,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
def verify():
 g=read(GATE);m=read(MANIFEST)
 assert g['status']=='READY_FOR_AUTHORIZED_CENTRAL_NUMERICAL_PILOT'
 assert g['independentMethodClearance'] is False
 assert g['proceedAuthority']=='EXPLICIT_USER_NO_FURTHER_GROK_MINOR_CORRECTIONS_AND_RUN'
 assert g['launcherSha256']==sha(__file__) and g['manifestSha256']==sha(MANIFEST)
 assert g['completionMessageSha256']==hashlib.sha256(MESSAGE.encode()).hexdigest()
 for path,h in g['sourcePins'].items():assert sha(path)==h,path
 for item in m['packets'].values():assert sha(item['path'])==item['sha256'],item['path']
 assert str(RUN)==m['runRoot'] and str(RUN/'complete')==m['resultDirectory']
 assert g['command'][1]==str(ROOT/'scripts/s11c_guarded_run.py')
 assert sha(m['endMaps']['database'])==m['endMaps']['sha256']
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
