#!/usr/bin/env python3
"""Validate the completed saved sensitivity stream; do not emit it again."""
import ast,copy,hashlib,json,marshal,shutil,time
from pathlib import Path
import S11c_d_remaining_case_first_jet_output_recover as prior
h=prior.h;f,m,engine=h.f,h.m,h.engine
ORIGIN=f.STORE/'s11c-remaining-case-first-jet-20260920/response/output/focused-recovery-01/complete'
PLAN=f.M/'S11c_d_remaining_case_first_jet_output_finish_plan.md'
PROOF=f.M/'S11c_d_remaining_case_first_jet_output_finish_focused.json'
native_emit=h.emit_part


def load(base,resume):
    old=json.loads((ORIGIN/'inputs.json').read_text());pins=dict(old['sourceFiles'])
    for n,v in pins.items():f.require(f.digest(f.ROOT/n)==f.digest(ORIGIN/'source'/n)==v,'exact original/current/frozen output sources')
    proof=json.loads(PROOF.read_text())
    f.require(proof['status']=='PASSED_SAVED_TRANSCRIPT_VALIDATION_WIRING' and proof['helperSha256']==f.digest(Path(__file__)) and proof['actualToolOutcome']['exitCode']==proof['actualToolOutcome']['stderrBytes']==0 and proof['actualToolOutcome']['checksStdoutIdentical'],'bounded actual saved stream/wiring regression')
    for p in (Path(__file__),PLAN,PROOF):pins[str(p.resolve().relative_to(f.ROOT))]=f.digest(p)
    manifest={k:old[k] for k in ('input','settings','scope','numericalOrigin','numericalChecksSha256','acceptedValidation')}
    manifest.update(runDirectory=str(base),sourceFiles=pins,inputPackets=dict(old['inputPackets']),copiedInputs={})
    for name in ('check.py','checks.json','validate.stdout','validate.stderr'):
        path=Path(proof['runDirectory'])/name;manifest['inputPackets'][str(path)]=f.digest(path)
    cp=json.loads(h.CP.read_text());f.require(cp['status']=='ACCEPTED_FIRST_JET_NUMERICAL_RESPONSES' and cp['checksSha256']==old['numericalChecksSha256'],'unchanged accepted numerical source')
    if resume:
        check=json.loads((resume/'checks.json').read_text());accepted=json.loads(h.FCP.read_text())
        f.require(accepted['status']=='ACCEPTED_FIRST_JET_OUTPUT_FOCUS' and f.digest(resume/'checks.json')==accepted['checksSha256'] and check['status']=='COMPLETED_FIRST_JET_OUTPUT_FOCUS','accepted finished output focus')
        f.require(check['sourceFiles']==pins,'identical validation/output sources')
        for n,v in check['inputPackets'].items():
            f.require(n not in manifest['inputPackets'] or manifest['inputPackets'][n]==v,'same exact consumed input');manifest['inputPackets'][n]=v
        for n,v in check['artifacts'].items():m.retain(resume/n,base/n,manifest,v['sha256'])
        for p in (resume/'checks.json',h.FCP):manifest['inputPackets'][str(p)]=f.digest(p)
        parts=check['parts'];manifest['completedEmissionReuse']={'directory':str(resume),'checksSha256':f.digest(resume/'checks.json'),'parts':list(parts),'artifacts':len(check['artifacts'])}
    else:
        outcome=json.loads((ORIGIN.parent/'active.json').read_text())
        f.require(outcome==json.loads((ORIGIN.parent/'first_jet_output.invocation.json').read_text())==json.loads((ORIGIN.parent/'supervisor.stdout').read_text()) and outcome['exitCode']==-14 and outcome['stderrBytes']==0 and not (ORIGIN/'checks.json').exists(),'actual timeout after saved output, not acceptance')
        f.require((ORIGIN.parent/'first_jet_output.stderr').stat().st_size==0,'empty original stderr')
        inventory={}
        for p in sorted(ORIGIN.rglob('*')):
            if not p.is_file():continue
            n=str(p.relative_to(ORIGIN));inventory[n]={'sha256':f.digest(p),'bytes':p.stat().st_size};target=n
            if n.startswith('source/'):target='timeout-source/'+n[len('source/'):]
            elif n=='inputs.json':target='timeout-inputs.json'
            elif n=='part-inventory.json':target='timeout-part-inventory.json'
            m.retain(p,base/target,manifest,inventory[n]['sha256'])
        for n,v in old['copiedInputs'].items():f.require(f.digest(ORIGIN/n)==v,'all completed original copies')
        for n in ('active.json','first_jet_output.invocation.json','first_jet_output.stderr','first_jet_output.stdout','supervisor.stdout','supervisor.stderr'):
            m.retain(ORIGIN.parent/n,base/'timeout-logs'/n,manifest)
        parts=json.loads((ORIGIN/'part-inventory.json').read_text())
        for name,record in parts.items():
            folder=base/record['directory'];f.require(f.digest(folder/'full.out')==record['sha256'],'all completed parts retained')
            check=json.loads((folder/'emission-checks.json').read_text())
            if record['kind']!='historical':
                packet='continuum-response.pickle' if record['kind']=='continuum' else 'continuum-currents.pickle'
                f.require(check['tags']==record['tags'] and len(check['keys'])==record['keys'] and check['metadataPaths']==record['metadataPaths'] and f.digest(folder/packet)==check['packetSha256'],'completed native part evidence')
        manifest['completedTranscriptReuse']={'directory':str(ORIGIN),'outcome':outcome,'artifacts':inventory,'newNumericalWork':0,'newSensitivityEmission':0}
        f.save(base/'transcript-reuse.json',manifest['completedTranscriptReuse'])
    for n,v in pins.items():
        dest=base/'source'/n;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/n,dest);f.require(f.digest(dest)==v,'frozen finish sources')
    f.save(base/'inputs.json',manifest);f.save(base/'part-inventory.json',parts)
    return manifest,tuple(cp['checks']['cases']),parts


def saved_keys_index(entries,prefix):
    tags=list(entries);base='PY_S11CD_'+prefix
    f.require(tags[-4:]==[base+'_WRITE_KEYS','PY_S11CD_METADATA_'+prefix+'_WRITE_KEYS',base+'_EMISSION_LINES','PY_S11CD_METADATA_'+prefix+'_EMISSION_LINES'],'entire saved terminal stream')
    # Reconstruct the original key algorithm and source-line map, not inferred values.
    keys={tag:'s11cd'+prefix+str(i) for i,tag in enumerate(tags[:-4]) if not tag.startswith('PY_S11CD_METADATA_')}
    f.require(engine.cas(keys)==entries[base+'_WRITE_KEYS'],'actual complete saved export keys')
    restored=h.restore_emission_index({str(k):v for k,v in entries[base+'_EMISSION_LINES']},tags[:-2])
    index=engine.emission_index(restored)
    f.require(engine.cas(index)==entries[base+'_EMISSION_LINES'],'actual complete saved line index')
    return keys,index


def validation_tail(label,folder):
    """Exact native decode/replay/metadata/hash suffix plus durable observations."""
    tree=ast.parse(Path(h.rnative.__file__).read_text());main=h.h.native.function(tree,'main')
    start=next(i for i,n in enumerate(main.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='entries' for t in n.targets))
    stop=next(i for i,n in enumerate(main.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='summary' for t in n.targets))
    original=copy.deepcopy(main.body[start:stop]);body=copy.deepcopy(original)
    # The same sole packet-address change as the accepted custom output driver.
    class Packet(ast.NodeTransformer):
        def visit_Constant(self,n):return ast.Constant('sensitivity-output.pickle') if n.value=='continuum-response.pickle' else n
    body=[Packet().visit(n) for n in body]
    before=copy.deepcopy(body)
    body.insert(2,ast.parse('keys,index=saved_keys_index(entries,PREFIX)').body[0])
    # Observers neither change operands nor replace any guard. Strip them to prove it.
    class Observe(ast.NodeTransformer):
        def visit_For(self,n):
            self.generic_visit(n)
            if ast.unparse(n.iter)=="grades.decoded_lines(base / 'full.out')":n.body.extend(ast.parse("mark('decoded',len(entries))").body)
            elif ast.unparse(n.iter)=='entries.items()':n.body.insert(0,ast.parse("mark('metadata',metadata_paths)").body[0])
            return n
        def visit_FunctionDef(self,n):
            self.generic_visit(n)
            if n.name=='replay':n.body.extend(ast.parse("mark('replayed',len(seen))").body)
            return n
    body=[Observe().visit(n) for n in body]
    out=[]
    for n in body:
        out.append(n)
        if isinstance(n,ast.For) and ast.unparse(n.iter)=="grades.decoded_lines(base / 'full.out')":out.extend(ast.parse("mark('decode_complete',len(entries),True)").body)
        if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and n.value.args and isinstance(n.value.args[-1],ast.Constant) and n.value.args[-1].value=='response payload/key census':out.extend(ast.parse("save_replay(entries,seen,keys,index)").body)
    body=out
    class Undo(ast.NodeTransformer):
        def visit_Expr(self,n):
            if isinstance(n.value,ast.Call) and isinstance(n.value.func,ast.Name) and n.value.func.id in ('mark','save_replay'):return None
            return self.generic_visit(n)
        def visit_Assign(self,n):
            if isinstance(n.value,ast.Call) and isinstance(n.value.func,ast.Name) and n.value.func.id=='saved_keys_index':return None
            return self.generic_visit(n)
    stripped=[Undo().visit(copy.deepcopy(n)) for n in body];stripped=[n for n in stripped if n is not None]
    f.require(ast.dump(ast.Module(body=stripped,type_ignores=[]))==ast.dump(ast.Module(body=before,type_ignores=[])),'whole native validation suffix, only saved key/index and progress wiring')
    reverse=copy.deepcopy(before)
    class Unpacket(ast.NodeTransformer):
        def visit_Constant(self,n):return ast.Constant('continuum-response.pickle') if n.value=='sensitivity-output.pickle' else n
    reverse=[Unpacket().visit(n) for n in reverse]
    f.require(ast.dump(ast.Module(body=reverse,type_ignores=[]))==ast.dump(ast.Module(body=original,type_ignores=[])),'whole native validation suffix packet reverse join')
    body.extend(ast.parse("return dict(tags=len(entries),keys=keys,metadataPaths=metadata_paths)").body)
    node=ast.FunctionDef(name='validate_saved',args=ast.arguments(posonlyargs=[],args=[ast.arg(n) for n in ('base','result','r','pins','operands','before')],vararg=None,kwonlyargs=[],kw_defaults=[],kwarg=None,defaults=[]),body=body,decorator_list=[])
    started=time.monotonic();last={}
    def mark(stage,count,force=False):
        bucket=count//256
        if not force and last.get(stage)==bucket:return
        last[stage]=bucket
        with (folder/'validation-progress.jsonl').open('a') as stream:stream.write(json.dumps({'stage':stage,'count':count,'seconds':time.monotonic()-started})+'\n')
    def save_replay(entries,seen,keys,index):
        f.save(folder/'completed-payload-replay.json',{'tags':len(entries),'seen':len(seen),'keys':len(keys),'transcriptSha256':f.digest(folder/'full.out'),'packetSha256':f.digest(folder/'sensitivity-output.pickle'),'nativeReplayPredicateUnchanged':True,'sourceHelperSha256':f.digest(Path(h.__file__)),'seconds':time.monotonic()-started})
        mark('replay_complete',len(seen),True)
    prefix='FIRST_JET_'+label.replace('__','_')+'_SENSITIVITY'
    fn=h.h.native.compile_function(node,dict(vars(h.rnative),emit_result=h.emit_sensitivity,PREFIX=prefix,saved_keys_index=saved_keys_index,mark=mark,save_replay=save_replay))
    join={'wholeNativeValidationSuffix':True,'packetAddressEdits':1,'savedKeysAndIndexOnly':True,'progressObserversOnly':True,'sourceSha256':f.digest(Path(h.rnative.__file__)),'emitterSha256':h.h.inputs.matrices.binding.native_body(h.emit_sensitivity),'nativeOutputMainBytecodeUnchanged':True}
    return fn,join


def emit_part(base,label,kind,manifest):
    folder=base/'parts'/(kind+'_'+label)
    if kind!='sensitivity' or not (folder/'full.out').exists():return native_emit(base,label,kind,manifest)
    packet='sensitivity-output.pickle';before=f.digest(folder/packet);transcript=f.digest(folder/'full.out')
    f.require(before==f.digest(ORIGIN/'parts'/folder.name/packet) and transcript==f.digest(ORIGIN/'parts'/folder.name/'full.out'),'exact completed saved bundle/transcript')
    result=f.unpickle(folder/packet);r=h.context(base,label,result);fn,join=validation_tail(label,folder)
    f.save(folder/'validation-join.json',join)
    check=fn(folder,result,r,manifest['sourceFiles'],manifest['inputPackets'],before)
    f.require(f.digest(folder/'full.out')==transcript,'no repeated individual emission or transcript change')
    return prior.finish_tail()(base,label,kind,manifest,folder,packet,result,join,before,check)


if __name__=='__main__':
    h.load=load;h.emit_part=emit_part;h.main()
