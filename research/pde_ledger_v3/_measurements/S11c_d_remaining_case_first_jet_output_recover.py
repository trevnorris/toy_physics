#!/usr/bin/env python3
"""Resume saved output after the explicit-integer unit annotation repair."""
import ast,copy,gc,io,json,shutil,sys,contextlib
from pathlib import Path
import S11c_d_remaining_case_first_jet_output as h
f,m,engine=h.f,h.m,h.engine
ORIGIN=f.STORE/'s11c-remaining-case-first-jet-20260920/response/output/focused/complete'
REPAIR=f.M/'S11c_d_remaining_case_first_jet_output_unit_repair.json'
PLAN=f.M/'S11c_d_remaining_case_first_jet_output_recovery_plan.md'
native_load=h.load;native_emit=h.emit_part

def load(base,resume):
    proof=json.loads(REPAIR.read_text())
    f.require(f.digest(Path(h.__file__))==proof['repairedSha256'] and f.digest(ORIGIN/'source/_measurements/S11c_d_remaining_case_first_jet_output.py')==proof['originalSha256'],'explicit original/current emitter source join')
    old=json.loads((ORIGIN/'inputs.json').read_text());pins=dict(old['sourceFiles'])
    for name,sha in pins.items():
        f.require(f.digest(ORIGIN/'source'/name)==sha,'immutable original emitter source')
        expected=proof['repairedSha256'] if name=='_measurements/S11c_d_remaining_case_first_jet_output.py' else sha
        f.require(f.digest(f.ROOT/name)==expected,'only recorded integer-unit repair');pins[name]=expected
    for p in (Path(__file__),PLAN,REPAIR):pins[str(p.resolve().relative_to(f.ROOT))]=f.digest(p)
    manifest={k:old[k] for k in ('input','settings','scope','numericalOrigin','numericalChecksSha256','acceptedValidation')}
    manifest.update(runDirectory=str(base),sourceFiles=pins,inputPackets=dict(old['inputPackets']),copiedInputs={},recovery={'originalDirectory':str(ORIGIN),'repairSha256':f.digest(REPAIR),'nativeOutputMainBytecodeUnchanged':True,'newNumericalWork':0})
    cp=json.loads(h.CP.read_text());f.require(cp['status']=='ACCEPTED_FIRST_JET_NUMERICAL_RESPONSES' and cp['checksSha256']==manifest['numericalChecksSha256'],'unchanged accepted numerical origin')
    labels=tuple(cp['checks']['cases'])
    if resume:
        check=json.loads((resume/'checks.json').read_text());accepted=json.loads(h.FCP.read_text())
        f.require(accepted['status']=='ACCEPTED_FIRST_JET_OUTPUT_FOCUS' and f.digest(resume/'checks.json')==accepted['checksSha256'],'accepted recovered output focus')
        f.require(check['sourceFiles']==pins and check['status']=='COMPLETED_FIRST_JET_OUTPUT_FOCUS','exact recovered focus helper joins')
        for n,v in check['inputPackets'].items():
            f.require(n not in manifest['inputPackets'] or manifest['inputPackets'][n]==v,'same consumed output input');manifest['inputPackets'][n]=v
        for n,v in check['artifacts'].items():m.retain(resume/n,base/n,manifest,v['sha256'])
        manifest['inputPackets'][str(resume/'checks.json')]=f.digest(resume/'checks.json');manifest['inputPackets'][str(h.FCP)]=f.digest(h.FCP)
        manifest['completedEmissionReuse']={'directory':str(resume),'checksSha256':f.digest(resume/'checks.json'),'parts':list(check['parts']),'artifacts':len(check['artifacts'])};parts=check['parts']
    else:
        outcome=json.loads((ORIGIN.parent/'active.json').read_text())
        f.require(outcome==json.loads((ORIGIN.parent/'first_jet_output.invocation.json').read_text())==json.loads((ORIGIN.parent/'supervisor.stdout').read_text()) and outcome['exitCode']==1 and outcome['status']=='failed' and not (ORIGIN/'checks.json').exists(),'actual failed output outcome')
        error=(ORIGIN.parent/'first_jet_output.stderr').read_text();f.require(error.endswith("TypeError: unsupported operand type(s) for -: 'Rational' and 'bool'\n"),'actual unit-offset error')
        saved={}
        for p in sorted(ORIGIN.rglob('*')):
            if not p.is_file():continue
            n=str(p.relative_to(ORIGIN));saved[n]={'sha256':f.digest(p),'bytes':p.stat().st_size}
            target=n
            if n.startswith('source/'):target='original-output-source/'+n[len('source/'):]
            elif n=='inputs.json':target='original-output-inputs.json'
            elif n=='part-inventory.json':target='original-part-inventory.json'
            elif n=='parts/sensitivity_'+h.FOCUS+'/full.out':target='parts/sensitivity_'+h.FOCUS+'/original-partial.out'
            elif n=='parts/sensitivity_'+h.FOCUS+'/emitter-join.json':target='parts/sensitivity_'+h.FOCUS+'/original-emitter-join.json'
            m.retain(p,base/target,manifest,saved[n]['sha256'])
        for n,sha in old['copiedInputs'].items():f.require(f.digest(ORIGIN/n)==sha,'original completed numerical copies')
        for n in ('active.json','first_jet_output.invocation.json','first_jet_output.stderr','first_jet_output.stdout','supervisor.stdout','supervisor.stderr'):
            p=ORIGIN.parent/n;m.retain(p,base/'original-logs'/n,manifest)
        parts=json.loads((ORIGIN/'part-inventory.json').read_text())
        for name,record in parts.items():
            folder=ORIGIN/record['directory'];f.require(f.digest(folder/'full.out')==record['sha256'],'completed original transcript hash')
            checks=json.loads((folder/'emission-checks.json').read_text())
            if record['kind']!='historical':
                f.require(checks['tags']==record['tags'] and checks['metadataPaths']==record['metadataPaths'] and len(checks['keys'])==record['keys'],'completed original replay census')
                f.require(f.digest(folder/('continuum-response.pickle' if record['kind']=='continuum' else 'continuum-currents.pickle'))==checks['packetSha256'],'completed original emitted input')
        manifest['recovery'].update(originalOutcome=outcome,originalArtifacts=saved,completedParts=list(parts))
        f.save(base/'output-reuse.json',manifest['recovery'])
    for n,sha in pins.items():
        dest=base/'source'/n;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/n,dest);f.require(f.digest(dest)==sha,'frozen recovery source')
    f.save(base/'inputs.json',manifest);f.save(base/'part-inventory.json',parts)
    return manifest,labels,parts

def finish_tail():
    tree=ast.parse(Path(h.__file__).read_text());source=h.h.native.function(tree,'emit_part')
    start=next(i for i,n in enumerate(source.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='check' for t in n.targets))+1
    body=copy.deepcopy(source.body[start:])
    f.require(ast.dump(ast.Module(body=body,type_ignores=[]))==ast.dump(ast.Module(body=source.body[start:],type_ignores=[])),'unchanged complete individual post-emission guard/control tail')
    args=('base','label','kind','manifest','folder','packet','result','join','before','check')
    node=ast.FunctionDef(name='finish_saved_part',args=ast.arguments(posonlyargs=[],args=[ast.arg(n) for n in args],vararg=None,kwonlyargs=[],kw_defaults=[],kwarg=None,defaults=[]),body=body,decorator_list=[])
    return h.h.native.compile_function(node,vars(h))

def emit_part(base,label,kind,manifest):
    folder=base/'parts'/(kind+'_'+label)
    if kind!='sensitivity' or not (folder/'original-partial.out').exists():return native_emit(base,label,kind,manifest)
    packet='sensitivity-output.pickle';result=f.unpickle(folder/packet);before=f.digest(folder/packet)
    f.require(before==f.digest(ORIGIN/'parts'/('sensitivity_'+label)/packet),'entire saved output bundle reused without reconstruction')
    fn,join=h.custom_emitter(label);join.update(savedBundleSha256=before,originalSourceSha256=json.loads(REPAIR.read_text())['originalSha256'],integerUnitRepairSha256=f.digest(REPAIR));f.save(folder/'emitter-join.json',join)
    r=h.context(base,label,result);raw=(folder/'original-partial.out').read_text().splitlines(keepends=True);decoded=list(h.decoded_lines(folder/'original-partial.out'))
    f.require(len(raw)==len(decoded)>0,'complete original prefix codec lines');original_emit=engine.emit;count=0
    def resume_emit(name,value):
        nonlocal count
        if count>=len(raw):return original_emit(name,value)
        tag,_,body=decoded[count].rstrip('\n').partition(': ')
        f.require(tag=='PY_S11CD_'+name and h.rnative.grades._restore(body)==engine.cas(value),'actual unchanged original prefix payload')
        stream=sys.stdout;buffer=io.StringIO()
        with contextlib.redirect_stdout(buffer):original_emit(name,value)
        f.require(buffer.getvalue()==raw[count],'exact original raw prefix codec state')
        stream.write(raw[count]);count+=1
        if count==len(raw):f.save(folder/'prefix-reuse.json',{'originalPartialSha256':f.digest(folder/'original-partial.out'),'identicalDecodedPayloads':count,'identicalRawPrefixLines':count,'savedBundleSha256':before,'newNumericalWork':0})
    engine.emit=resume_emit
    try:check=fn(folder,result,r,manifest['sourceFiles'],manifest['inputPackets'],before)
    finally:engine.emit=original_emit
    f.require(count==len(raw),'all original prefix payloads replayed')
    return finish_tail()(base,label,kind,manifest,folder,packet,result,join,before,check)

if __name__=='__main__':
    h.load=load;h.emit_part=emit_part;h.main()
