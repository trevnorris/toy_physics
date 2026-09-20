#!/usr/bin/env python3
"""Uniform case controls from accepted full current and modal source families."""
import argparse,ast,copy,json,resource,shutil,signal,time
from pathlib import Path
import numpy as np
import sympy as sp
import S11c_d_uniform_response as native
import S11c_d_remaining_case_boundary as boundary
import S11c_d_remaining_case_response as h
import S11c_d_remaining_case_response_finish as output
f,b,m=h.f,h.b,h.m
engine=f.engine
PLAN=f.M/'S11c_d_remaining_case_uniform_plan.md'
MCP=f.M/'S11c_d_remaining_case_modes_checkpoint.json'
UCP=f.M/'S11c_d_uniform_response_checkpoint.json'
FCP=f.M/'S11c_d_remaining_case_uniform_focused.json'
BASELINE=h.BASELINE
NEW='LAB_HELD__RHOBR_CONSTANT'
ENDS=('REFERENCE','LEFT','RIGHT')
INPUT_KEYS=('binding','physical','relation','rootPacket','pairing','acoustic','fieldUnits','frequency','depth')


def emitter(label):
    tree=ast.parse(Path(native.__file__).read_text());emit=copy.deepcopy(h.function(tree,'emit_result'));main=h.function(tree,'main')
    first=next(i for i,n in enumerate(main.body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and ast.unparse(n.value.func)=='engine.EMISSION_LINES.clear')
    stop=next(i for i,n in enumerate(main.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='summary' for t in n.targets))
    original=copy.deepcopy(main.body[first:stop]);body=copy.deepcopy(original);prefix='s11cd'+label.replace('__','_')+'UniformResponse';changes=[]
    class Keys(ast.NodeTransformer):
        def visit_Constant(self,n):
            if n.value=='s11cdUniformResponse':changes.append(True);return ast.Constant(prefix)
            return n
    body=[Keys().visit(n) for n in body];f.require(len(changes)==1,'one uniform export namespace')
    class Undo(ast.NodeTransformer):
        def visit_Constant(self,n):return ast.Constant('s11cdUniformResponse') if n.value==prefix else n
    restored=[Undo().visit(copy.deepcopy(n)) for n in body]
    f.require(ast.dump(ast.Module(body=restored,type_ignores=[]))==ast.dump(ast.Module(body=original,type_ignores=[])),'whole original uniform output and replay tail')
    result=ast.Dict(keys=[ast.Constant(k) for k in ('tags','keys','metadataPaths')],values=[ast.Call(ast.Name('len',ast.Load()),[ast.Name('entries',ast.Load())],[]),ast.Name('keys',ast.Load()),ast.Name('paths',ast.Load())])
    node=ast.FunctionDef(name='emit_case',args=ast.arguments(posonlyargs=[],args=[ast.arg(n) for n in ('base','result','r','pins','operands','before')],vararg=None,kwonlyargs=[],kw_defaults=[],kwarg=None,defaults=[]),body=body+[ast.Return(result)],decorator_list=[])
    namespace=dict(vars(native),PREFIX='UNIFORM_RESPONSE_'+label.replace('__','_'));namespace['emit_result']=h.compile_function(emit,namespace)
    return h.compile_function(node,namespace)


def load(base,resume=None):
    cp,origin=h.checked_checkpoint(MCP,'ACCEPTED_FOUR_CASE_MODE_SUBSPACES');uc,ur=h.checked_checkpoint(UCP,'PUBLISHED_ANNEX_VERIFIED')
    for accepted in (cp,uc):
        for n,v in accepted['inputPackets'].items():f.require(f.digest(Path(n))==v,'accepted uniform original input')
    pins=dict(cp['sourceFiles'])
    for n,v in uc['sourceFiles'].items():f.require(n not in pins or pins[n]==v,'same consumed uniform source');pins[n]=v
    for p in (Path(__file__).resolve(),PLAN,MCP,UCP,Path(native.__file__),Path(boundary.__file__),Path(output.__file__)):
        pins[str(p.resolve().relative_to(f.ROOT))]=f.digest(p)
    manifest={'runDirectory':str(base),'sourceFiles':pins,'inputPackets':{str(origin/'checks.json'):cp['checksSha256'],str(ur/'checks.json'):uc['checksSha256']},'copiedInputs':{},'input':cp['input'],
      'scope':'Four case-specific three-background homogeneous controls; ten original response reuses and one new RHOBR-right matching shared after full input joins. No profile quadrature, mode/current closure or continuum expansion.'}
    for accepted in (cp,uc):
        for n,v in accepted['inputPackets'].items():
            f.require(n not in manifest['inputPackets'] or manifest['inputPackets'][n]==v,'shared original uniform operand');manifest['inputPackets'][n]=v
    if resume:
        accepted=json.loads(FCP.read_text());focused=json.loads((resume/'checks.json').read_text())
        f.require(accepted['status']=='ACCEPTED_FOUR_CASE_UNIFORM_INPUTS' and accepted['checksSha256']==f.digest(resume/'checks.json'),'accepted focused uniform operands')
        f.require(focused['mode']=='focused' and focused['sourceFiles']==pins,'exact focused helper and sources')
        for n,v in focused['artifacts'].items():m.retain(resume/n,base/n,manifest,v['sha256'])
        manifest['completedFocusedReuse']={'directory':str(resume),'checksSha256':f.digest(resume/'checks.json'),'artifacts':len(focused['artifacts'])}
    else:
        for n,v in cp['artifacts'].items():m.retain(origin/n,base/n,manifest,v['sha256'])
        for n,v in uc['artifacts'].items():m.retain(ur/n,base/'accepted-uniform'/n,manifest,v['sha256'])
        m.retain(ur/'checks.json',base/'accepted-uniform/checks.json',manifest,uc['checksSha256'])
        m.retain(ur/'inputs.json',base/'accepted-uniform/inputs.json',manifest)
        for n in ('full.out','checks.json','uniform-response.pickle'):
            m.retain(base/'accepted-uniform'/n,base/'cases'/BASELINE/'continuum'/n,manifest)
    for n,v in pins.items():
        dst=base/'source'/n;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/n,dst);f.require(f.digest(dst)==v,'frozen uniform source')
    f.save(base/'inputs.json',manifest)
    return manifest


def root_equal(a,c):
    # BOUND_CARRIERS is the native mapping export. Preserve both original tuple
    # orders; compare its unique actual keys/values and every other root field.
    copies=[]
    for value in (a,c):
        copied=copy.deepcopy(value);items=copied['packet']['BOUND_CARRIERS'];mapping=dict(items)
        f.require(len(mapping)==len(items),'unique native bound-carrier keys')
        copied['packet']['BOUND_CARRIERS']=mapping;copies.append(copied)
    return m.same(*copies)


def rational_pencil_pair(target,pair):
    # The two native forms may have differently grouped rational denominators.
    # Save full original expressions and raw difference before exact numerator
    # certificates. No branch, root equation or numerical tolerance is used.
    a,c=pair;f.require(a.shape==c.shape,'actual rational pencil dimensions')
    f.atomic_pickle(target/'uniform-original-pencil-raw.pickle',{'pair':pair,'liveStrings':(str(a),str(c)),'rawDifference':a-c})
    entries=[]
    for index,(left,right) in enumerate(zip(a,c)):
        together=sp.together(left-right);numerator,denominator=together.as_numer_denom();expanded=sp.expand(numerator)
        original_denominators=tuple((v.base,v.exp) for source in (left,right) for v in sp.preorder_traversal(source) if isinstance(v,sp.Pow) and v.exp.is_negative)
        entries.append({'index':index,'left':left,'right':right,'together':together,'numerator':numerator,'denominator':denominator,'expandedNumerator':expanded,'originalDenominators':original_denominators})
    f.atomic_pickle(target/'uniform-original-pencil-certificate.pickle',entries)
    f.require(all(v['expandedNumerator']==0 for v in entries),'exact rational pencil numerator identity')
    return sp.ImmutableMatrix(a.rows,a.cols,[v['expandedNumerator'] for v in entries])


def input_routes(base):
    values=f.unpickle(base/'remaining-case-currents.pickle');old=f.unpickle(base/'accepted-uniform/uniform-response.pickle');old_inputs=json.loads((base/'accepted-uniform/inputs.json').read_text())
    inventory={};controls=[];labels=tuple(values['cases'])
    for label in labels:
        for end in ENDS:
            target=base/'cases'/label/end.lower();inp=f.unpickle(target/'mode-inputs.pickle');modal,known=f.unpickle(target/'modal.pickle')
            owner=NEW if end=='RIGHT' and label.endswith('RHOBR_CONSTANT') else BASELINE
            ref=f.unpickle(base/'cases'/owner/end.lower()/'mode-inputs.pickle');reference,unused=f.unpickle(base/'cases'/owner/end.lower()/'modal.pickle')
            pairs={k:(inp[k],ref[k]) for k in INPUT_KEYS};f.atomic_pickle(target/'uniform-input-pairs.pickle',pairs)
            f.require((inp['case'],inp['end'])==(label,end),'actual case/background input address')
            f.require(all(root_equal(a,c) if key=='rootPacket' else m.same(a,c) for key,(a,c) in pairs.items()),'complete uniform family inputs')
            f.require(m.same(modal,reference),'complete reused current normalized modal packet')
            value=values['cases'][label][end]
            f.require(m.same(inp['pairing'],value['pairing']) and m.same(inp['acoustic'],value['acoustic']),'actual current source before reuse')
            root=f.unpickle(base/'native-roots'/(label+'__'+end+'.pickle'))
            f.require(m.same(inp['rootPacket'],root) and m.same(modal['NATIVE_RECORDS'],root['records']),'complete actual root/lift candidates')
            f.require(all(complex(r['OMEGA'])==complex(inp['frequency']) and r['NULLITY']==int(n['NULLITY']) for r,n in zip(modal['RECORDS'],root['records'])),'actual frequency and complete nullity')
            f.require(tuple(inp['fieldUnits'])==tuple(old['fieldUnits']),'full inherited field units')
            reuse=owner==BASELINE
            if reuse:
                original=next(Path(n) for n in old_inputs['inputPackets'] if n.endswith('/'+end.lower()+'/modal.pickle'))
                f.require(f.digest(original)==old_inputs['inputPackets'][str(original)]==f.digest(base/'accepted-modes'/end.lower()/'modal.pickle')==f.digest(target/'modal.pickle'),'same actual original full modal input')
                symbolic=f.unpickle(base/'accepted-uniform'/(end.lower()+'-symbolic.pickle'))
                pair=(inp['physical'],symbolic['freshPencil']);residual=rational_pencil_pair(target,pair)
                f.atomic_pickle(target/'uniform-original-pencil-pair.pickle',{'pair':pair,'residual':residual})
                f.require(len(old['backgrounds'][end]['modes'])==len(modal['RECORDS']),'original homogeneous full candidate census')
            inventory[label+'__'+end]={'owner':owner+'__'+end,'reusedOriginalUniformResponse':reuse,'candidates':len(modal['RECORDS']),'basisDirections':sum(r['NULLITY'] for r in modal['RECORDS']),'inputSha256':f.digest(target/'mode-inputs.pickle'),'modalSha256':f.digest(target/'modal.pickle')}
    # Mutate actual input pairs without altering any saved operand.
    ref=f.unpickle(base/'cases'/NEW/'right/mode-inputs.pickle')
    for key in ('physical','fieldUnits','rootPacket','pairing'):
        changed=copy.deepcopy(ref)
        if key=='physical':changed[key]=sp.ImmutableMatrix(changed[key])+sp.eye(5)
        elif key=='fieldUnits':changed[key]=list(changed[key]);changed[key][0]=(2,0,0)
        elif key=='rootPacket':changed[key]['records'][0]['NORMAL_LIFT_SIGN']=-changed[key]['records'][0]['NORMAL_LIFT_SIGN']
        else:
            original=changed[key]['SLAB_CURRENT_MATRIX'];changed[key]['SLAB_CURRENT_MATRIX']=2*original
        f.require(not (root_equal(changed[key],ref[key]) if key=='rootPacket' else m.same(changed[key],ref[key])),'actual wrong uniform source/address rejection');controls.append(key)
    f.require(len(inventory)==12 and sum(v['reusedOriginalUniformResponse'] for v in inventory.values())==10,'full actual uniform reuse census')
    result={'routes':inventory,'wrongInputControls':controls,'labels':labels};f.save(base/'input-routes.json',result);return result


def prepare(base,manifest):
    target=base/'new-uniform';target.mkdir()
    sources=f.unpickle(base/'remaining-case-end-sources.pickle');values=f.unpickle(base/'remaining-case-currents.pickle')
    pair,inp,value=m.context(base,NEW,'RIGHT',sources,values,manifest)
    inputs=f.unpickle(base/'cases'/NEW/'right/mode-inputs.pickle');modal,known=f.unpickle(base/'cases'/NEW/'right/modal.pickle')
    # Restore accepted symbolic/scalar preparation; no native prepare/construct call.
    builder=boundary.restore_builder(pair,value,modal,inputs)
    keys=('PENCIL_PLUS','NORMAL_PENCIL_PLUS','FREQUENCY_PENCIL_PLUS','CURRENT_SLAB','CURRENT_BULK','INFINITE_DEPTH_INTEGRAL')
    symbolic={**builder.symbolic_operands,**builder.scalar_operands};bound={k:symbolic[k].xreplace(inputs['binding']) for k in keys}
    saved={'variables':builder.variables,'bound':bound,'physical':inputs['physical'],'relation':inputs['relation'],'binding':inputs['binding'],'sourceModalSha256':f.digest(base/'cases'/NEW/'right/modal.pickle'),'dimensionState':dict(engine.PHYSICAL_METADATA.dimensions.__dict__)}
    f.atomic_pickle(target/'prepared-evaluation.pickle',saved)
    evaluate={k:builder.evaluate[k] for k in keys};fresh=sp.lambdify((pair.modes.k,pair.modes.q),inputs['physical'],'numpy',cse=True)
    curve=sp.lambdify((pair.modes.k,pair.modes.q),inputs['relation'].xreplace(inputs['binding']),'numpy',cse=True)
    records=[];controls=[]
    for old in modal['RECORDS']:
        i=old['INDEX'];k,q,pq=(complex(old[n]) for n in ('K','Q','PHYSICAL_Q'));w=complex(old['OMEGA']);point=(w,w,k.conjugate(),k,pq.conjugate(),pq)
        actual={n:np.asarray(evaluate[n](*point),complex) for n in keys if n!='INFINITE_DEPTH_INTEGRAL'}
        matrix=np.asarray(fresh(k,q),complex);right=old['FORMS']['RIGHT'];left=old['FORMS']['LEFT'];n=old['NULLITY'];scale=1+b.norm(matrix)
        residuals={name:(v-old['OPERANDS'][name])/(1+b.norm(old['OPERANDS'][name])) for name,v in actual.items()}
        residuals.update(freshPencil=(matrix-actual['PENCIL_PLUS'])/scale,rightKernel=matrix@right/scale,leftKernel=matrix.conj().T@left/scale,
                         rightGram=right.conj().T@right-np.eye(n),leftGram=left.conj().T@left-np.eye(n))
        info={name:old[name] for name in ('INDEX','ROOT_DISK_INDEX','NORMAL_LIFT_SIGN','K','Q','PHYSICAL_Q','OMEGA','NULLITY','SHEET_MEMBERSHIP','EXACT_REAL_NORMAL','BULK_DECAY_DISK_CERTIFIED','PHYSICAL_RIGHT_CURRENT_NORMALIZATION_DEFINED')}
        info['CLASSIFIER_STATUS']=modal['NATIVE_RECORDS'][i]['CLASSIFIER_STATUS']
        record={'info':info,'right':right,'rawRight':right,'left':left,'singularValues':old['SINGULAR_VALUES'],'pencil':matrix,
                'normalPairing':old['FORMS']['N_NORMAL'],'frequencyPairing':old['FORMS']['N_FREQUENCY'],'normalRank':old['NORMAL_PAIRING_RANK'],'frequencyRank':old['FREQUENCY_PAIRING_RANK'],
                'currentOperands':{key:actual[key] for key in ('CURRENT_SLAB','CURRENT_BULK')},'residuals':residuals,'waveResidual':complex(curve(k,q)),'currentDefined':bool(old['PHYSICAL_RIGHT_CURRENT_NORMALIZATION_DEFINED'])}
        residuals['normalPairing']=left.conj().T@actual['NORMAL_PENCIL_PLUS']@right-record['normalPairing']
        residuals['frequencyPairing']=left.conj().T@actual['FREQUENCY_PENCIL_PLUS']@right-record['frequencyPairing']
        if old['BULK_DECAY_DISK_CERTIFIED']:
            depth=complex(evaluate['INFINITE_DEPTH_INTEGRAL'](*point));current=actual['CURRENT_SLAB']+depth*actual['CURRENT_BULK'];gram=right.conj().T@current@right
            record.update(depthIntegral=depth,currentMatrix=current,currentGram=gram,currentEigenvalues=old['FORMS']['CURRENT_EIGENVALUES'])
            residuals.update(depthIntegral=np.asarray([depth-old['INFINITE_DEPTH_INTEGRAL']]),currentGram=gram-old['FORMS']['CURRENT_INFINITE'],currentHermitian=(gram-gram.conj().T)/(1+b.norm(gram)))
            if record['currentDefined']:
                record.update(fluxRight=old['FORMS']['FLUX_RIGHT'],fieldToFlux=old['FORMS']['FIELD_TO_FLUX_MAP'],signedCurrent=old['FORMS']['SIGNED_CURRENT'])
                residuals['fluxNormalization']=record['fluxRight'].conj().T@current@record['fluxRight']-record['signedCurrent']
                mutation=record['fluxRight'].conj().T@(2*current)@record['fluxRight']-record['signedCurrent'];f.require(b.norm(mutation)>0,'actual current normalization mutation');controls.append(b.norm(mutation))
        f.atomic_pickle(target/('right-mode-'+str(i)+'.pickle'),record)
        f.require(b.norm(residuals)<1e-8 and abs(record['waveResidual'])<1e-7*(1+abs(q)**2+abs(k)**2),'saved full basis and actual source replay')
        f.require(m.same(right,old['FORMS']['RIGHT']) and m.same(left,old['FORMS']['LEFT']),'no new modal frame solve')
        records.append(record)
    f.require(len(records)==18 and controls,'all candidate dispositions and responding physical-current mutation')
    f.atomic_pickle(target/'prepared-modes.pickle',records)
    result={'candidates':len(records),'basisDirections':sum(r['info']['NULLITY'] for r in records),'maximumResidual':max(b.norm(r['residuals']) for r in records),'currentMutations':controls,'newModeConstructions':0,'newCurrentClosures':0,'newMatchingSolves':0}
    f.save(target/'preparation-checks.json',result);return result


def construct(base,manifest,routes):
    target=base/'new-uniform';prepared=f.unpickle(target/'prepared-evaluation.pickle');records=f.unpickle(target/'prepared-modes.pickle')
    evaluate={k:sp.lambdify(prepared['variables'],v,'numpy',cse=True) for k,v in prepared['bound'].items()}
    result=native.match(target,'RIGHT',records,evaluate)
    old=f.unpickle(base/'accepted-uniform/uniform-response.pickle');cases={}
    for label in routes['labels']:
        if label==BASELINE:continue
        backgrounds={end:old['backgrounds'][end] if routes['routes'][label+'__'+end]['reusedOriginalUniformResponse'] else {'modes':records,'response':result} for end in ENDS}
        packet={**old,'backgrounds':backgrounds,'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets'],'scope':manifest['scope']+' This case: '+label}
        directory=base/'cases'/label/'continuum';directory.mkdir(parents=True,exist_ok=True);f.atomic_pickle(directory/'uniform-response.pickle',packet);cases[label]=packet
    return cases


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True);ap.add_argument('--focused',action='store_true');ap.add_argument('--resume-focused',type=Path);args=ap.parse_args()
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);start=time.monotonic();base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    manifest=load(base,args.resume_focused)
    if args.resume_focused:routes=json.loads((base/'input-routes.json').read_text());prepared=json.loads((base/'new-uniform/preparation-checks.json').read_text())
    else:routes=input_routes(base);prepared=prepare(base,manifest)
    # Compile all actual output adapters before a matching solve.
    for label in routes['labels']:
        if label!=BASELINE:emitter(label)
    emissions={};aggregate=None;new=0
    if not args.focused:
        cases=construct(base,manifest,routes);new=1
        sources=f.unpickle(base/'remaining-case-end-sources.pickle');values=f.unpickle(base/'remaining-case-currents.pickle')
        pair,inp,value=m.context(base,NEW,'RIGHT',sources,values,manifest)
        state=f.unpickle(base/'new-uniform/prepared-evaluation.pickle')['dimensionState'];engine.PHYSICAL_METADATA.dimensions.__dict__.update(state)
        for label,result in cases.items():
            target=base/'cases'/label/'continuum';emissions[label]=emitter(label)(target,result,pair.r,manifest['sourceFiles'],manifest['inputPackets'],f.digest(target/'uniform-response.pickle'));f.save(target/'emission-checks.json',emissions[label])
        with (base/'original-combined.out').open('xb') as stream:
            for label in routes['labels']:stream.write((base/'cases'/label/'continuum/full.out').read_bytes())
        aggregate=output.aggregate(base,routes['labels']);f.save(base/'aggregation-checks.json',aggregate)
        f.atomic_pickle(base/'remaining-case-uniform.pickle',{'routes':routes,'cases':cases,'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets'],'scope':manifest['scope']})
    for n,v in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/n)==f.digest(base/'source'/n)==v,'current/frozen uniform source pre/post')
    for n,v in manifest['inputPackets'].items():f.require(f.digest(Path(n))==v,'original uniform input pre/post')
    for n,v in manifest['copiedInputs'].items():f.require(f.digest(base/n)==v,'immutable uniform input copies')
    artifacts={str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*') if p.is_file() and 'source' not in p.relative_to(base).parts and p!=base/'inputs.json' and p!=base/'checks.json'}
    checks={**manifest,'status':'PASSED_FOUR_CASE_UNIFORM_INPUTS' if args.focused else 'COMPLETED_FOUR_CASE_UNIFORM_RESPONSES','mode':'focused' if args.focused else 'construct','routes':routes,'prepared':prepared,'emissions':emissions,'aggregate':aggregate,'newMatchingSolves':new,'newModeConstructions':0,'newCurrentClosures':0,'artifacts':artifacts,'wallSeconds':time.monotonic()-start}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
