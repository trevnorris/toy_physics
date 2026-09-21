#!/usr/bin/env python3
"""Sequential, bounded saved-output phases with one native worker at a time."""
import argparse,ast,copy,hashlib,json,os,resource,shutil,signal,subprocess,sys,time,types
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PLAN=ROOT/'_measurements/S11c_d_remaining_case_first_jet_output_staged_plan.md'
PROOF=ROOT/'_measurements/S11c_d_remaining_case_first_jet_output_staged_focused.json'
FCP=ROOT/'_measurements/S11c_d_remaining_case_first_jet_output_focused.json'
NCP=ROOT/'_measurements/S11c_d_remaining_case_first_jet_response_checkpoint.json'
BASELINE='LAB_HELD__RHO4_CONSTANT'


def sha(path):
    result=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):result.update(block)
    return result.hexdigest()


def save(path,value):
    path.write_text(json.dumps(value,indent=2)+'\n')


def schedule(labels,parts):
    assert len(labels)==len(set(labels))==4 and BASELINE in labels
    wanted={'historical_'+BASELINE}
    for label in labels:
        assert label in ('LAB_HELD__RHO4_CONSTANT','LAB_HELD__RHOBR_CONSTANT','MATERIAL_ADVECTED__RHO4_CONSTANT','MATERIAL_ADVECTED__RHOBR_CONSTANT')
        for kind in ('continuum','current','sensitivity'):
            if kind=='continuum' and label==BASELINE:continue
            wanted.add(kind+'_'+label)
    assert len(wanted)==12 and set(parts)<=wanted
    for name,record in parts.items():
        assert record['directory']=='parts/'+name and name==record['kind']+'_'+record['case']
        assert record['kind'] in ('historical','continuum','current','sensitivity')
        if record['kind']=='historical':assert record['case']==BASELINE
    result=[{'phase':'prepare'}]
    for label in labels:
        for kind in ('continuum','current','sensitivity'):
            name=kind+'_'+label
            if kind=='continuum' and label==BASELINE or name in parts:continue
            for phase in (('emit','validate') if kind=='sensitivity' else ('native',)):
                result.append({'phase':phase,'case':label,'kind':kind,'part':name})
    result.append({'phase':'aggregate'})
    return result,wanted


def native_functions(h):
    """Extract only the existing numerical prohibition and final audit bodies."""
    main=h.h.native.function(ast.parse(Path(h.__file__).read_text()),'main')
    first=next(i for i,n in enumerate(main.body) if isinstance(n,ast.FunctionDef) and n.name=='forbidden')
    body=copy.deepcopy(main.body[first:first+2])
    node=ast.FunctionDef(name='disable_numerical',args=ast.arguments(posonlyargs=[],args=[],vararg=None,kwonlyargs=[],kw_defaults=[],kwarg=None,defaults=[]),body=body,decorator_list=[])
    assert ast.dump(ast.Module(body=body,type_ignores=[]))==ast.dump(ast.Module(body=main.body[first:first+2],type_ignores=[]))
    disable=h.h.native.compile_function(node,vars(h))
    first=next(i for i,n in enumerate(main.body) if isinstance(n,ast.For) and ast.unparse(n.iter)=="manifest['sourceFiles'].items()")
    body=copy.deepcopy(main.body[first:])
    assert ast.dump(ast.Module(body=body,type_ignores=[]))==ast.dump(ast.Module(body=main.body[first:],type_ignores=[]))
    args=('base','manifest','parts','combined','args','start')
    node=ast.FunctionDef(name='finish_all',args=ast.arguments(posonlyargs=[],args=[ast.arg(n) for n in args],vararg=None,kwonlyargs=[],kw_defaults=[],kwarg=None,defaults=[]),body=body,decorator_list=[])
    finish=h.h.native.compile_function(node,vars(h))
    return disable,finish,{'wholeNumericalProhibitionBody':True,'wholeFinalHashArtifactCheckStdoutTail':True,'sourceSha256':sha(h.__file__)}


def emission_prefix(h,label):
    """Split the original driver at entries={}, without editing the emitter."""
    main=h.h.native.function(ast.parse(Path(h.rnative.__file__).read_text()),'main')
    first=next(i for i,n in enumerate(main.body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and ast.unparse(n.value.func)=='engine.EMISSION_LINES.clear')
    stop=next(i for i,n in enumerate(main.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='entries' for t in n.targets))
    original=copy.deepcopy(main.body[first:stop]);body=copy.deepcopy(original)
    prefix='FIRST_JET_'+label.replace('__','_')+'_SENSITIVITY';key='s11cd'+prefix;changes=[]
    class Rename(ast.NodeTransformer):
        def visit_Constant(self,n):
            if n.value=='s11cdContinuumResponse':changes.append(n.value);return ast.Constant(key)
            return n
    body=[Rename().visit(n) for n in body]
    class Undo(ast.NodeTransformer):
        def visit_Constant(self,n):return ast.Constant('s11cdContinuumResponse') if n.value==key else n
    reverse=[Undo().visit(copy.deepcopy(n)) for n in body]
    assert len(changes)==1 and ast.dump(ast.Module(body=reverse,type_ignores=[]))==ast.dump(ast.Module(body=original,type_ignores=[]))
    body.extend(ast.parse('return dict(keys=keys,index=index)').body)
    node=ast.FunctionDef(name='emit_only',args=ast.arguments(posonlyargs=[],args=[ast.arg(n) for n in ('base','result','r')],vararg=None,kwonlyargs=[],kw_defaults=[],kwarg=None,defaults=[]),body=body,decorator_list=[])
    fn=h.h.native.compile_function(node,dict(vars(h.rnative),emit_result=h.emit_sensitivity,PREFIX=prefix))
    return fn,{'wholeNativeEmissionPrefix':True,'namespaceEdits':1,'emitterUnchanged':True,'emitterSha256':h.h.inputs.matrices.binding.native_body(h.emit_sensitivity),'nativeSourceSha256':sha(h.rnative.__file__)}


def worker(args):
    # Scientific modules are imported only inside the sole active bounded worker.
    import S11c_d_remaining_case_first_jet_output_finish as x
    h,f=x.h,x.f;base=args.run_directory.resolve();phase_dir=args.phase_directory.resolve()
    base.relative_to(f.STORE);assert phase_dir.parent==base.parent/'phases'
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);start=time.monotonic()
    disable,finish,native_join=native_functions(h);disable()
    proof=json.loads(PROOF.read_text());assert proof['status']=='PASSED_STAGED_SAVED_OUTPUT_WIRING' and proof['helperSha256']==sha(__file__)
    if args.worker_phase=='prepare':
        base.mkdir(exist_ok=False);manifest,labels,parts=x.load(base,args.focus.resolve())
        # Preserve the copied four-part inventory separately from the new growing one.
        copied_inventory=manifest['copiedInputs'].pop('part-inventory.json')
        assert sha(base/'part-inventory.json')==copied_inventory and not (base/'accepted-focus-part-inventory.json').exists()
        (base/'part-inventory.json').rename(base/'accepted-focus-part-inventory.json')
        manifest['copiedInputs']['accepted-focus-part-inventory.json']=copied_inventory
        save(base/'part-inventory.json',parts)
        original_pins=dict(manifest['sourceFiles'])
        for path in (Path(__file__),PLAN,PROOF):
            n=str(path.resolve().relative_to(f.ROOT));v=sha(path);manifest['sourceFiles'][n]=v;dest=base/'source'/n;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,dest);assert sha(dest)==v
        manifest['stagedOutput']={'focusDirectory':str(args.focus.resolve()),'focusChecksSha256':sha(args.focus/'checks.json'),'acceptedSourceFiles':original_pins,'focusedArtifactRoutes':{'part-inventory.json':'accepted-focus-part-inventory.json'},'workerWallBudgetSeconds':900,'maximumActiveNativeWorkers':1,'reason':'Measured saved sensitivity replay alone took 511 seconds; split new emission and replay, retain each completed part.'}
        for name in ('checks.json','check.py','validate.stdout','validate.stderr'):
            path=Path(proof['runDirectory'])/name;manifest['inputPackets'][str(path)]=sha(path)
        plan,wanted=schedule(labels,parts);save(base/'staged-plan.json',{'phases':plan,'wantedParts':sorted(wanted),'completedFocusedParts':list(parts),'nativeJoins':native_join})
        save(base/'inputs.json',manifest)
        summary={'status':'PREPARED_ACCEPTED_OUTPUT_REUSE','focusedArtifacts':len(manifest['copiedInputs']),'completedParts':list(parts),'plannedPhases':len(plan)}
    else:
        manifest=json.loads((base/'inputs.json').read_text());parts=json.loads((base/'part-inventory.json').read_text())
        for n,v in manifest['sourceFiles'].items():assert sha(f.ROOT/n)==sha(base/'source'/n)==v
        if args.worker_phase=='aggregate':
            plan=json.loads((base/'staged-plan.json').read_text());assert set(parts)==set(plan['wantedParts']) and len(parts)==12
            for name,record in parts.items():assert sha(base/record['directory']/'full.out')==record['sha256']
            combined=h.aggregate(base,parts);save(base/'aggregation-checks.json',combined)
            finish(base,manifest,parts,combined,types.SimpleNamespace(mode='construct'),start)
            return
        assert args.case in tuple(json.loads(NCP.read_text())['checks']['cases'])
        name=args.kind+'_'+args.case;assert name not in parts
        folder=base/'parts'/name
        if args.worker_phase=='native':
            assert args.kind in ('continuum','current') and not (args.case==BASELINE and args.kind=='continuum')
            record=h.emit_part(base,args.case,args.kind,manifest);parts[name]=record
            summary={'status':'COMPLETED_NATIVE_PART','part':name,'record':record}
        elif args.worker_phase=='emit':
            assert args.kind=='sensitivity';folder.mkdir(exist_ok=False)
            result=h.bundle(base,args.case,manifest);packet=folder/'sensitivity-output.pickle';f.atomic_pickle(packet,result);before=sha(packet)
            r=h.context(base,args.case,result);fn,join=emission_prefix(h,args.case);save(folder/'emitter-join.json',join)
            emitted=fn(folder,result,r);f.atomic_pickle(folder/'emission-state.pickle',emitted)
            assert sha(packet)==before and not h.engine.PHYSICAL_METADATA.dimensions.constraints
            summary={'status':'COMPLETED_EMISSION_AWAITING_REPLAY','part':name,'packetSha256':before,'transcriptSha256':sha(folder/'full.out'),'keys':len(emitted['keys']),'tags':len(h.engine.EMISSION_LINES),'emitterJoin':join}
            save(folder/'emission-complete.json',summary)
        elif args.worker_phase=='validate':
            assert args.kind=='sensitivity';emitted=json.loads((folder/'emission-complete.json').read_text());packet='sensitivity-output.pickle'
            before=sha(folder/packet);transcript=sha(folder/'full.out');assert before==emitted['packetSha256'] and transcript==emitted['transcriptSha256']
            result=f.unpickle(folder/packet);r=h.context(base,args.case,result);fn,join=x.validation_tail(args.case,folder);save(folder/'validation-join.json',join)
            check=fn(folder,result,r,manifest['sourceFiles'],manifest['inputPackets'],before)
            assert check['tags']==emitted['tags'] and len(check['keys'])==emitted['keys'] and sha(folder/'full.out')==transcript
            record=x.prior.finish_tail()(base,args.case,args.kind,manifest,folder,packet,result,join,before,check);parts[name]=record
            summary={'status':'COMPLETED_SAVED_SENSITIVITY_REPLAY','part':name,'record':record,'repeatedEmissions':0}
        else:raise ValueError(args.worker_phase)
        save(base/'inputs.json',manifest);save(base/'part-inventory.json',parts)
    summary.update(nativeJoins=native_join,sourceSha256=sha(__file__),newSolves=0,newQuadrature=0,newCurrentContractions=0,wallSeconds=time.monotonic()-start)
    save(phase_dir/'checks.json',summary);signal.alarm(0);print(json.dumps(summary,indent=2))


def coordinator(args):
    base=args.run_directory.resolve();outer=base.parent;assert not base.exists()
    cp=json.loads(FCP.read_text());assert cp['status']=='ACCEPTED_FIRST_JET_OUTPUT_FOCUS' and Path(cp['runDirectory'])==args.focus.resolve() and sha(args.focus/'checks.json')==cp['checksSha256']
    focused=json.loads((args.focus/'checks.json').read_text());numerical=json.loads(NCP.read_text());assert numerical['status']=='ACCEPTED_FIRST_JET_NUMERICAL_RESPONSES'
    planned,wanted=schedule(tuple(numerical['checks']['cases']),focused['parts'])
    assert len(focused['parts'])==4 and len(wanted-set(focused['parts']))==8 and len(planned)==13
    current=[None]
    def stop(signum,_frame):
        if current[0] is not None:
            try:os.killpg(current[0].pid,signal.SIGTERM)
            except ProcessLookupError:pass
        raise SystemExit(128+signum)
    signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
    phases=outer/'phases';phases.mkdir(exist_ok=False);records=[];begin=time.monotonic();source_sha=sha(__file__)
    save(outer/'staged-schedule.json',{'phases':planned,'wantedParts':sorted(wanted),'phaseBudgetSeconds':900,'maxNativeWorkers':1,'automaticRetries':0,'maximumChildSeconds':len(planned)*900,'sourceSha256':source_sha})
    for i,entry in enumerate(planned):
        assert sha(__file__)==source_sha
        directory=phases/(str(i).zfill(2)+'_'+entry['phase']+('_'+entry['part'] if 'part' in entry else ''));directory.mkdir()
        command=[sys.executable,'-u',str(Path(__file__).resolve()),'--worker-phase',entry['phase'],'--phase-directory',str(directory),'--run-directory',str(base),'--focus',str(args.focus.resolve())]
        if 'case' in entry:command+=['--case',entry['case'],'--kind',entry['kind']]
        if (base/'inputs.json').exists():shutil.copyfile(base/'inputs.json',directory/'input-manifest.json')
        save(directory/'command.json',{'command':command,'entry':entry,'sourceSha256':source_sha})
        threads={k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')};start=time.monotonic()
        with (directory/'stdout').open('xb') as out,(directory/'stderr').open('xb') as err:
            process=subprocess.Popen(command,cwd=ROOT,env=dict(os.environ,**threads),stdin=subprocess.DEVNULL,stdout=out,stderr=err,start_new_session=True,close_fds=True);current[0]=process
            try:code=process.wait(timeout=915)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid,signal.SIGKILL);code=process.wait()
            finally:current[0]=None
        record={'phase':entry,'directory':str(directory),'command':command,'exitCode':code,'stderrBytes':(directory/'stderr').stat().st_size,'wallSeconds':time.monotonic()-start,'sourceSha256':source_sha,'stdoutSha256':sha(directory/'stdout'),'stderrSha256':sha(directory/'stderr')}
        records.append(record);save(directory/'outcome.json',record);save(outer/'stage-outcomes.json',records)
        assert code==0 and record['stderrBytes']==0,('stopped; preserve completed outputs',record)
        check=base/'checks.json' if entry['phase']=='aggregate' else directory/'checks.json'
        assert check.read_bytes()==(directory/'stdout').read_bytes(),('phase checks/stdout identity',entry)
        if entry['phase']!='aggregate':
            receipt=base/'phase-receipts';receipt.mkdir(exist_ok=True);path=receipt/(directory.name+'.json');assert not path.exists();save(path,record)
            manifest=json.loads((base/'inputs.json').read_text())
            for p in directory.iterdir():
                if p.is_file():manifest['inputPackets'][str(p)]=sha(p)
            save(base/'inputs.json',manifest)
    assert len(records)==len(planned) and all(v['exitCode']==v['stderrBytes']==0 for v in records)
    save(outer/'staged-completion.json',{'status':'COMPLETED_ALL_SAVED_OUTPUT_PHASES','phases':len(records),'maximumConcurrentWorkers':1,'automaticRetries':0,'wallSeconds':time.monotonic()-begin,'checksSha256':sha(base/'checks.json'),'sourceSha256':source_sha})
    sys.stdout.write((base/'checks.json').read_text())


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',required=True,type=Path);ap.add_argument('--focus',required=True,type=Path)
    ap.add_argument('--worker-phase',choices=('prepare','native','emit','validate','aggregate'));ap.add_argument('--phase-directory',type=Path);ap.add_argument('--case');ap.add_argument('--kind');args=ap.parse_args()
    if args.worker_phase:worker(args)
    else:coordinator(args)
