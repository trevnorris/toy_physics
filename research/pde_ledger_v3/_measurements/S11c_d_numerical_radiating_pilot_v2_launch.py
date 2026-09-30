#!/usr/bin/env python3
"""Single hook-first guarded launch, no scientific imports, retries or polling."""
from datetime import datetime,timezone
import hashlib,json,os,select,shutil,subprocess,sys,time,traceback
from pathlib import Path
ROOT=Path('/var/projects/toy_physics');M=ROOT/'research/pde_ledger_v3/_measurements'
PREFIX='S11c_d_numerical_radiating_pilot_v2'
GATE=M/(PREFIX+'_gate.json');MANIFEST=M/(PREFIX+'_inputs.json')
RUN=ROOT/'_scratch/s11c/s11c-d-numerical-radiating-20260929/central-pilot-v2'
HOOK=ROOT/'scripts/codex_job_watch.py';THREAD='01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE='The user-authorized central omega3 numerical radiating integration/finite-pilot job has completed or reported an error.\nRoot: /var/projects/toy_physics/_scratch/s11c/s11c-d-numerical-radiating-20260929/central-pilot-v2.\nInspect actual coordinator/launch/stage-start, shared guard duration/effective limits/resource samples/outcomes/logs, supervisor numerical_radiating_pilot_v2 invocation/active, strict scientific stderr, stdout/checks byte identity, checks/failure/failed-check evidence and incomplete operation, every complete journal input/return and result/partial-array receipt, source snapshots and all posthashes. Exit code is not science acceptance. Scientific restoration only inside containment; JSON/hash/source/opaque SQLite bytes may be inspected lightly. Preserve complete/failed work, no automatic scientific retry.\nThe preceding omega3 end maps passed for both ends/all four contrasts:101.309s,110.375MiB,79operations,40intact inputs. This worker restores those results, frequency-live original source records, accepted native domain binding and continuum grade/occurrence records; no producer/root/path/current/face replay. It binds omega3 and contrasts, joins all375source addresses/units/160native occurrences, checks actual source-root lifts and endpoint denominator orders, two-sided weighted endpoint values, independent full middle-leg physical-k quadrature, responsive wrong-sheet/Jacobian controls and uniform eW/theta Gaussian operator versus original REFERENCE pencil. Only after these pass does it assemble the declared finite Chebyshev systems and solve at most10central cases, retaining exact current cross terms. The four numerical settings cover contrast scaling, uniform background, refinement, larger domain, regulator and negative-deficit controls. Inspect actual operands/values; source hashes or completed status alone are insufficient. All previous accepted end states remain unchanged.\nNo wall/native/CPU-time/inactivity deadline. Shared scripts/s11c_guarded_run.py around existing normalization supervisor, RuntimeMaxUSec=infinity/Restart=no,2GiB native/cgroup,zero swap,one CPU,32tasks,one thread,4GiBhost reserve,locking/overlap refusal. No scheduler changes, fallback or retry. Immutable partial matrix results are stored after4096batches; CPU work is not acceptance. Matrix/source array sharing requires actual identical operands and bounds peak memory.\nUser authorized the numerical goal, removed time caps and explicitly directed no further Grok pass, minor corrections and run. Literal round2method history remains Claude CLEAR with5minor corrections; Grok no final report. No paired method/build/result clearance or new review is claimed. If substantive central checks fail, preserve the exact evidence and stop before neighbors; do not silently change physics. If the central controls support proceeding, the standing goal covers omega2.3/4 after actual readiness, reusing saved results, with max30finite solves total. Do not ask repetitive stage permission. Finite-model deficit is not automatically physical leakage; retained-order loss interpretation, calibration, Green/FORM/A11/A12 remain open. Exact omega1 stays parked.\nPreserve historical results/failures/incidents/review debt, Lean/S11_lean and protected builder suffixf01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2. Canonical progress outside scratch; scratch never committed. No model polling or recurring tasks. Hook targets01a0e01b-ef84-7192-817f-584cda5d339b.\n\nStartup v1 is preserved: four complete input restorations, zero numerical operations/finite solves; failed on initializing the integration helper namespace after context inspection. New v2 initializes before first use and restores all four exact saved returns without their functions, joins actual input bytes and copies all11prior files/22730208bytes unchanged. The consumed non-journal end-map evidence is reused byte-identically. Both numerical libraries and the entire run body after initialization are unchanged. Three stdlib tests reproduce the old failure, pass the corrected order, verify no replay and reject an input mismatch. This is a local tooling correction, not a scientific retry or fresh method clearance. Inspect these runtime joins and all prior hashes.\n'
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
