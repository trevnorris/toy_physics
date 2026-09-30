#!/usr/bin/env python3
"""Single hook-first guarded launch, no scientific imports, retries or polling."""
from datetime import datetime,timezone
import hashlib,json,os,select,shutil,subprocess,sys,time,traceback
from pathlib import Path
ROOT=Path('/var/projects/toy_physics');M=ROOT/'research/pde_ledger_v3/_measurements'
PREFIX='S11c_d_numerical_radiating_uniform_refinement'
GATE=M/(PREFIX+'_gate.json');MANIFEST=M/(PREFIX+'_inputs.json')
RUN=ROOT/'_scratch/s11c/s11c-d-numerical-radiating-20260929/uniform-refinement-01'
HOOK=ROOT/'scripts/codex_job_watch.py';THREAD='01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE='The user-directed focused omega3 uniform Gaussian refinement has completed or reported an error.\nRoot: /var/projects/toy_physics/_scratch/s11c/s11c-d-numerical-radiating-20260929/uniform-refinement-01.\nInspect actual coordinator/launch, guard enforced duration/resources/outcomes/logs, supervisor numerical_radiating_uniform_refinement invocation/active, strict stderr, stdout/checks identity, failure/comparisons, complete journal input/return receipts, preserved source posthashes and snapshots. Metadata/hash/opaque-byte reads only outside containment; no scientific restoration. Exit status alone is not acceptance.\nUser directed focused refinement at Gaussian momentum0.6 against saved reference, unchanged equations/tolerance, reusing all completed nested integrals. This worker restores selected actual v4 returns/operands, keeps all47prior files in place pinned, and computes only five uniform row actions and two source-column assemblies at outer orders32/48/64 with512 source nodes and64 with768 source nodes. No producer/root/path/current/face/middle-integral/independent-reference replay or finite solve. Inspect actual reference/source/matrix/units/positions/settings joins, complete arrays and saved error/tolerance ratios. Unsaved ephemeral source quadrature arrays are explicitly recreated and saved. All earlier failures and36complete v4 operations remain intact.\nIf original tolerance is supported, standing user authorization covers proceeding toward omega3 finite energy balance and existing uniform/scaling/resolution/domain/regulator/sign controls, reusing completed integration work. If not, user requested a predeclared coarse-resolution balance instead of further repair cycles; do not equate operator residual to flux uncertainty. A nominal1e-5 reporting floor still requires actual finite-current controls/envelope and assessment by the requested missing independent method leg before interpreting any deficit. No physical leakage bound follows. Check replacement-method-review status/consent before export; no peer sharing or automatic new packet. A real source/integrity/method failure is preserved and reported.\nOne scientific worker under unchanged shared guard/supervisor, no wall/native/CPU/inactivity deadline,2GiB native/cgroup,zero swap,one CPU,32tasks,one thread,4GiBhost reserve,global lock/overlap refusal,RuntimeMaxUSec=infinity/Restart=no. Hook to01a0e01b-ef84-7192-817f-584cda5d339b. No polling, fallback, scientific retry or scheduler change. All prior results/incident/review history,Lean/S11_lean and protected suffixf01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2 remain intact. Scratch never committed. Exactomega1 stays parked;physical loss/calibration/Green/FORM/A11/A12 stay open.\n'
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
