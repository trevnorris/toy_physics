#!/usr/bin/env python3
"""Fixed saved-output sequence; stdlib coordinator, one guarded worker at a time."""
import argparse,fcntl,hashlib,json,os,resource,shutil,signal,subprocess,sys,time
from datetime import datetime,timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[3]
LEDGER=ROOT/'research/pde_ledger_v3'
M=LEDGER/'_measurements';STORE=ROOT/'_scratch/s11c'
FCP=M/'S11c_d_remaining_case_coordinate_output_focused.json'
WORKER=M/'S11c_d_remaining_case_coordinate_output_recover.py'
PLAN=M/'S11c_d_remaining_case_coordinate_output_staged_plan.md'
PROOF=M/'S11c_d_remaining_case_coordinate_output_staged_wiring.json'
GUARD=ROOT/'scripts/s11c_guarded_run.py'
SUPERVISOR=M/'S11c_d_end_normalization_run.py'
BASELINE='LAB_HELD__RHO4_CONSTANT'
THREADS={n:'1' for n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS')}


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024**2),b''):digest.update(block)
    return digest.hexdigest()


def save(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_name(path.name+'.new');temporary.write_text(json.dumps(value,indent=2)+'\n');temporary.replace(path)


def identical(a,b):
    if a.stat().st_size!=b.stat().st_size:return False
    with a.open('rb') as x,b.open('rb') as y:
        while True:
            left=x.read(1024**2);right=y.read(1024**2)
            if left!=right:return False
            if not left:return True


def schedule(labels):
    assert set(labels)=={a+'__'+b for a in ('LAB_HELD','MATERIAL_ADVECTED') for b in ('RHO4_CONSTANT','RHOBR_CONSTANT')}
    parts=[{'case':case,'kind':kind,'part':kind+'_'+case} for case in labels for kind in ('continuum','current','coordinate') if not(kind=='continuum' and case==BASELINE)]
    return [{'stage':'prepare'}]+[dict(v,stage='emit') for v in parts]+[dict(v,stage='replay') for v in parts]+[{'stage':'aggregate'}]


def check_schedule(entries,labels):
    assert entries==schedule(labels) and len(entries)==24,'exact ordered output sequence'
    for stage in ('emit','replay'):
        parts=[n['part'] for n in entries if n['stage']==stage]
        assert len(parts)==len(set(parts))==11
    assert {n['part'] for n in entries if n['stage']=='emit'}=={n['part'] for n in entries if n['stage']=='replay'}


def worker_command(base,focus,directory,entry):
    command=[sys.executable,'-u',str(WORKER),'--run-directory',str(base),'--stage',entry['stage'],'--phase-directory',str(directory)]
    if entry['stage']=='prepare':command+=['--resume-from',str(focus)]
    elif entry['stage'] in ('emit','replay'):command+=['--case',entry['case'],'--kind',entry['kind']]
    else:assert entry=={'stage':'aggregate'}
    return command


def inspect_guard(directory,stage):
    inv=json.loads((directory/(stage+'.invocation.json')).read_text());guard=json.loads((directory/'resource-guard/outcome.json').read_text())
    assert inv==json.loads((directory/'active.json').read_text())==json.loads((directory/'resource-guard/stdout').read_text())
    assert guard==json.loads((directory/'guard.stdout').read_text())
    assert guard['childOutcome']==json.loads((directory/'resource-guard/child-outcome.json').read_text())
    assert inv['exitCode']==inv['stderrBytes']==guard['exitCode']==guard['stderrBytes']==guard['childOutcome']['exitCode']==0
    assert guard['limitsVerified'] and guard['childOutcome']['guardReason'] is None
    assert all((directory/n).stat().st_size==0 for n in (stage+'.stderr','guard.stderr','resource-guard/stderr'))
    limits=json.loads((directory/'resource-guard/effective-limits.json').read_text())
    assert json.loads((directory/'resource-guard/limit-validation.json').read_text())=={'verified':True,'actual':limits}
    assert limits['memory.max']=='2147483648' and limits['memory.swap.max']=='0' and limits['pids.max']=='32'
    assert limits['nice']>=15 and len(limits['affinity'])==1 and limits['threads']==THREADS
    peak=0;swap=0;host=None;cap=0
    for line in (directory/'resource-guard/resource-samples.jsonl').read_text().splitlines():
        v=json.loads(line);e=dict(n.split() for n in v['memory.events'].splitlines());assert e['oom']==e['oom_kill']=='0' and e.get('oom_group_kill','0')=='0' and v['memory.swap.current']=='0'
        peak=max(peak,int(v['memory.peak']));swap=max(swap,int(v['memory.swap.current']));cap=max(cap,int(e['max']))
        host=v['hostAvailableBytes'] if host is None else min(host,v['hostAvailableBytes'])
    assert host is not None
    return {'supervisor':inv,'guard':guard,'limits':limits,'memoryPeakBytes':peak,'swapBytes':swap,'capReclaims':cap,'minimumHostAvailableBytes':host,'oomEvents':0}


def coordinate(args):
    base=args.run_directory.resolve();base.relative_to(STORE);outer=base.parent
    assert not base.exists() and not (STORE/'PAUSED_HOST_FREEZE.json').exists()
    resource.setrlimit(resource.RLIMIT_AS,(512*1024**2,512*1024**2));os.nice(max(0,15-os.getpriority(os.PRIO_PROCESS,0)))
    cp=json.loads(FCP.read_text());focus=Path(cp['runDirectory']);assert focus==args.focus.resolve()
    assert cp['status']=='ACCEPTED_CASE_MATERIAL_OUTPUT_INPUTS' and sha(focus/'checks.json')==cp['checksSha256']
    validation=cp['validation'];vr=Path(validation['runDirectory']);assert sha(vr/'checks.json')==validation['checksSha256']
    inspected=inspect_guard(vr,'validate');assert identical(vr/'checks.json',vr/'validate.stdout')
    labels=tuple(cp['cases']);entries=schedule(labels);check_schedule(entries,labels)
    proof=json.loads(PROOF.read_text());assert proof['status']=='PASSED_STATIC_GUARDED_OUTPUT_SEQUENCE'
    assert proof['coordinatorSha256']==sha(__file__) and proof['planSha256']==sha(PLAN) and proof['focusChecksSha256']==cp['checksSha256']
    pins=dict(proof['workerSourcePins']);pins.update({str(Path(__file__).resolve()):sha(__file__),str(PLAN):sha(PLAN),str(PROOF):sha(PROOF),str(FCP):sha(FCP)})
    for n,v in pins.items():assert sha(n)==v
    phases=outer/'phases';phases.mkdir(exist_ok=False);records=[];started=time.monotonic();current=[None]
    save(outer/'staged-schedule.json',{'entries':entries,'nativeWorkersAtOnce':1,'phaseSeconds':900,'automaticRetries':0,'sourcePins':pins,'focus':str(focus),'checksSha256':cp['checksSha256'],'validatorChecksSha256':validation['checksSha256']})
    def stop(signum,_frame):
        if current[0] is not None:
            current[0].send_signal(signal.SIGINT)
            current[0].wait()
        raise SystemExit(128+signum)
    signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
    for index,entry in enumerate(entries):
        for n,v in pins.items():assert sha(n)==v,('changed pinned source',n)
        assert not (STORE/'PAUSED_HOST_FREEZE.json').exists()
        name=str(index).zfill(2)+'_'+entry['stage']+('_'+entry['part'] if 'part' in entry else '')
        directory=phases/name;directory.mkdir();native=worker_command(base,focus,directory,entry)
        assert native==proof['workerCommands'][index],'exact reviewed worker call and output address'
        command=[sys.executable,str(GUARD),'--log-directory',str(directory/'resource-guard'),'--seconds','900','--',
            sys.executable,str(SUPERVISOR),'--run-root',str(directory),'--stage','output_phase','--',*native]
        save(directory/'command.json',{'entry':entry,'command':command,'workerCommand':native,'sourcePins':pins})
        if (base/'inputs.json').exists():shutil.copyfile(base/'inputs.json',directory/'input-manifest.json')
        tick=time.monotonic()
        with (directory/'guard.stdout').open('xb') as out,(directory/'guard.stderr').open('xb') as err:
            process=subprocess.Popen(command,cwd=ROOT,env=dict(os.environ,**THREADS),stdin=subprocess.DEVNULL,stdout=out,stderr=err,start_new_session=True,close_fds=True);current[0]=process
            code=process.wait();current[0]=None
        record={'index':index,'entry':entry,'directory':str(directory),'command':command,'actualGuardProcessExitCode':code,'wallSeconds':time.monotonic()-tick,'accepted':False}
        records.append(record);save(directory/'outcome.json',record);save(outer/'stage-outcomes.json',records)
        assert code==0,('phase failed; retain saved output and stop without retry',record)
        record['resourceGuard']=inspect_guard(directory,'output_phase')
        check=base/'checks.json' if entry['stage']=='aggregate' else directory/'checks.json'
        assert identical(check,directory/'output_phase.stdout'),'actual phase checks/stdout identity'
        checks=json.loads(check.read_text());assert checks['stage']==entry['stage']
        assert checks['status']==('COMPLETED_FOUR_CASE_MATERIAL_OUTPUT' if entry['stage']=='aggregate' else 'COMPLETED_MATERIAL_OUTPUT_PHASE')
        if entry['stage']=='prepare':
            assert len(checks['artifacts'])>=938 and len(checks['parts'])==1
            assert checks['completedOutputInputReuse']['checksSha256']==cp['checksSha256'] and checks['completedOutputInputReuse']['artifacts']==938
            for n,v in cp['artifacts'].items():
                target='accepted-focus-part-inventory.json' if n=='part-inventory.json' else n
                assert sha(base/target)==v['sha256']
        elif entry['stage']=='emit':
            folder=base/'parts'/entry['part'];done=json.loads((folder/'emission-complete.json').read_text())
            assert done['status']=='COMPLETED_EMISSION_AWAITING_REPLAY' and sha(folder/'full.out')==done['transcriptSha256']
            assert entry['part'] not in checks['parts']
        elif entry['stage']=='replay':
            part=checks['parts'][entry['part']];folder=base/part['directory']
            assert part['kind']==entry['kind'] and part['case']==entry['case'] and sha(folder/'full.out')==part['sha256']
            done=json.loads((folder/'emission-complete.json').read_text());assert done['transcriptSha256']==part['sha256']
            assert (folder/'emission-checks.json').is_file() and (folder/'replay-controls.json').is_file()
        else:assert len(checks['parts'])==12 and checks['aggregate'] is not None
        record.update(accepted=True,checksSha256=sha(check),stdoutSha256=sha(directory/'output_phase.stdout'))
        save(directory/'outcome.json',record);save(outer/'stage-outcomes.json',records)
        if entry['stage']!='aggregate':
            receipt=base/'phase-receipts';receipt.mkdir(exist_ok=True);save(receipt/(name+'.json'),record)
            manifest=json.loads((base/'inputs.json').read_text())
            if entry['stage']=='prepare':
                additions={}
                for p in (Path(__file__).resolve(),PLAN,PROOF):
                    n=str(p.relative_to(LEDGER));assert n not in manifest['sourceFiles']
                    path=base/'source'/n;path.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,path)
                    assert sha(path)==pins[str(p)];manifest['sourceFiles'][n]=pins[str(p)];additions[n]=pins[str(p)]
                save(base/'staged-source-joins.json',{'acceptedFocusedSources':cp['sourceFiles'],'newCoordinatorSources':additions,'sourceFiles':manifest['sourceFiles'],'noAcceptedSourceChanged':True})
                save(base/'staged-plan.json',json.loads((outer/'staged-schedule.json').read_text()))
                manifest['stagedCoordinator']={'sourceFiles':additions,'fixedPhases':24,'guardedWorkerBudgetSeconds':900,'maximumConcurrentWorkers':1,'automaticRetries':0,'noScientificCoordinatorImports':True}
                manifest['acceptedOutputCases']=cp['cases']
            for p in directory.rglob('*'):
                if p.is_file():manifest['inputPackets'][str(p)]=sha(p)
            manifest['inputPackets'][str(outer/'staged-schedule.json')]=sha(outer/'staged-schedule.json')
            save(base/'inputs.json',manifest)
        del checks
    assert len(records)==24 and all(v['accepted'] for v in records)
    for n,v in pins.items():assert sha(n)==v
    save(outer/'staged-completion.json',{'status':'COMPLETED_ALL_GUARDED_MATERIAL_OUTPUT_PHASES','phases':24,'maximumConcurrentWorkers':1,'automaticRetries':0,'wallSeconds':time.monotonic()-started,'checksSha256':sha(base/'checks.json'),'sourcePins':pins})
    with (base/'checks.json').open('rb') as source:shutil.copyfileobj(source,sys.stdout.buffer)


def supervise(args):
    outer=args.run_directory.resolve().parent;outer.relative_to(STORE)
    lock=(STORE/'material-output-coordinator.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    assert not (outer/'active.json').exists() and not (outer/'staged-completion.json').exists()
    stage='coordinate_output_construct';stdout=outer/(stage+'.stdout');stderr=outer/(stage+'.stderr')
    command=[sys.executable,'-u',str(Path(__file__).resolve()),'--run-directory',str(args.run_directory.resolve()),'--focus',str(args.focus.resolve())]
    record={'stage':stage,'command':command,'startedUtc':datetime.now(timezone.utc).isoformat(),'stdout':str(stdout),'stderr':str(stderr),'status':'running','scientificWorkersGuardedIndividually':True};started=time.monotonic()
    with stdout.open('xb') as out,stderr.open('xb') as err:
        process=subprocess.Popen(command,cwd=ROOT,stdin=subprocess.DEVNULL,stdout=out,stderr=err,start_new_session=True,close_fds=True)
        record['childPid']=process.pid;save(outer/'active.json',record)
        def stop(signum,_frame):process.send_signal(signum);process.wait();raise SystemExit(128+signum)
        signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
        code=process.wait()
    record.update(exitCode=code,status='completed' if code==0 else 'failed',stderrBytes=stderr.stat().st_size,wallSeconds=time.monotonic()-started,finishedUtc=datetime.now(timezone.utc).isoformat())
    save(outer/'active.json',record);save(outer/(stage+'.invocation.json'),record);print(json.dumps(record,indent=2))
    return code


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True);parser.add_argument('--focus',type=Path,required=True);parser.add_argument('--supervise',action='store_true');args=parser.parse_args()
    if args.supervise:sys.exit(supervise(args))
    coordinate(args)
