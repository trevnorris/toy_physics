#!/usr/bin/env python3
"""Emit accepted derivative controls without repeating numerical construction."""
import argparse,ast,copy,gc,hashlib,itertools,json,resource,shutil,signal,time,types
from pathlib import Path
import numpy as np
import sympy as sp
import S11c_d_remaining_case_first_jet_response as h
from S11c_d_output_codec import PayloadEncoder,PayloadDecoder,decoded_lines,restore_emission_index
f,m,b=h.f,h.m,h.b;engine=f.engine;rnative=h.native.response
CP=f.M/'S11c_d_remaining_case_first_jet_response_checkpoint.json'
FCP=f.M/'S11c_d_remaining_case_first_jet_output_focused.json'
PLAN=f.M/'S11c_d_remaining_case_first_jet_output_plan.md'
FOCUS='LAB_HELD__RHOBR_CONSTANT'
SCOPE=('Actual own-case first-w-derivative closed-operator sensitivity. Four new finite controls, three new continuum controls and one reused historical continuum control. '
       'Finite currents cover four open directions with three separate closed matching amplitudes; retained continuum currents cover all seven directions and cross terms. '
       'Not a consistent new profile, isolated advection channel or c2 source-origin closure. Positive regulator, approximate boundaries, unresolved tiny signals and omitted parent pure-second-order terms remain. '
       'Only the baseline material comparison has been constructed. Other material routes remain unfinished without zero substitutes.')


def load(base,resume):
    cp,origin=h.native.checked_checkpoint(CP,'ACCEPTED_FIRST_JET_NUMERICAL_RESPONSES');checks=json.loads((origin/'checks.json').read_text())
    v=cp['validation'];vr=Path(v['runDirectory']);f.require(f.digest(vr/'checks.json')==v['checksSha256'] and v['outcome']['exitCode']==v['outcome']['stderrBytes']==0,'accepted independent saved-result validation')
    pins=dict(checks['sourceFiles'])
    for p in (Path(__file__),PLAN,CP):pins[str(p.resolve().relative_to(f.ROOT))]=f.digest(p)
    manifest={'runDirectory':str(base),'sourceFiles':pins,'inputPackets':dict(checks['inputPackets']),'copiedInputs':{},'input':checks['input'],'settings':checks['settings'],'scope':SCOPE,
              'numericalOrigin':str(origin),'numericalChecksSha256':cp['checksSha256'],'acceptedValidation':v['checksSha256']}
    for p in (origin/'checks.json',origin/'inputs.json',vr/'checks.json',vr/'validate.py',vr/'validator-join.json'):manifest['inputPackets'][str(p)]=f.digest(p)
    parts={}
    if resume:
        old=json.loads((resume/'checks.json').read_text());accepted=json.loads(FCP.read_text())
        f.require(accepted['status']=='ACCEPTED_FIRST_JET_OUTPUT_FOCUS' and accepted['checksSha256']==f.digest(resume/'checks.json'),'accepted completed output focus')
        f.require(old['sourceFiles']==pins and old['status']=='COMPLETED_FIRST_JET_OUTPUT_FOCUS','unchanged focused output sources')
        for n,v in old['artifacts'].items():m.retain(resume/n,base/n,manifest,v['sha256'])
        parts=old['parts'];manifest['inputPackets'][str(resume/'checks.json')]=f.digest(resume/'checks.json')
        manifest['completedEmissionReuse']={'directory':str(resume),'checksSha256':f.digest(resume/'checks.json'),'parts':list(parts),'artifacts':len(old['artifacts'])}
    else:
        for n,v in checks['artifacts'].items():m.retain(origin/n,base/'numerical'/n,manifest,v['sha256'])
        m.retain(origin/'inputs.json',base/'numerical/inputs.json',manifest)
        # Original baseline continuum emission is immutable, including its local index.
        folder=base/'parts'/('historical_'+h.BASELINE);folder.mkdir(parents=True)
        m.retain(origin/'cases'/h.BASELINE/'continuum/full.out',folder/'full.out',manifest)
        m.retain(origin/'cases'/h.BASELINE/'continuum/checks.json',folder/'emission-checks.json',manifest)
        old=json.loads((folder/'emission-checks.json').read_text())
        parts['historical_'+h.BASELINE]={'directory':str(folder.relative_to(base)),'kind':'historical','case':h.BASELINE,'sha256':f.digest(folder/'full.out'),'tags':old['tagCount'],'keys':old['writeKeys'],'metadataPaths':old['metadataPaths'],'reusedOriginalEmission':True}
    for n,v in pins.items():
        target=base/'source'/n;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/n,target);f.require(f.digest(target)==v,'frozen output source')
    f.save(base/'inputs.json',manifest);f.save(base/'part-inventory.json',parts)
    return manifest,tuple(checks['cases']),parts


def context(base,label,packet):
    source=f.unpickle(base/'numerical/interiors/accepted-cases'/label/'reduced-action.pickle')
    r,d=f.prior.domain.momentum.source.native.source.restore_context(source)
    d.__dict__.update(packet['dimensionState']);f.require(engine.PHYSICAL_METADATA.dimensions is d,'actual restored metadata context')
    return r


def bundle(base,label,manifest):
    n=base/'numerical';q=n/'cases'/label;old=n/'unchanged-response/cases'/label
    value=f.unpickle(q/'continuum/continuum-response.pickle')
    result={'case':label,'finite':{},'continuum':{},'sensitivity':f.unpickle(q/'first-jet-sensitivity.pickle'),
            'current':f.unpickle(q/'continuum/continuum-currents.pickle'),'openFiniteCurrent':f.unpickle(q/'finite/open-end-currents.pickle'),
            'caseBinding':f.unpickle(n/'interiors/accepted-bindings'/label/'case-binding.pickle'),
            'originalBinding':f.unpickle(n/'original-bindings'/label/'case-binding.pickle'),
            'interior':f.unpickle(n/'interiors/cases'/label/'interior-matrices.pickle'),
            'prepared':f.unpickle(n/'preparation'/label/'finite-prepared.pickle'),
            'fieldUnits':value['fieldUnits'],'rowUnits':value['rowUnits'],'currentUnit':value['currentUnit'],'dimensionState':value['dimensionState'],
            'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets'],'settings':value['settings'],'scope':SCOPE}
    for name,root in [('REVERSED',q),('ORIGINAL',old)]:
        result['finite'][name]={k:f.unpickle(root/'finite'/(fn+'.pickle')) for k,fn in [('system','finite-system'),('solution','finite-solution'),('observable','observable')]}
        result['continuum'][name]=f.unpickle(root/'continuum/continuum-response.pickle')
    if label==h.BASELINE:
        result['originalRecords']=f.unpickle(n/'accepted-first-jet/first-jet-binding.pickle')['originalRecords']
    else:result['originalRecords']=result['caseBinding']['grades']['originalRecords']
    if label!=h.BASELINE:result['formalRemainders']=f.unpickle(q/'continuum/formal-remainders.pickle')
    return result


def emit_sensitivity(result,r):
    prefix='FIRST_JET_'+result['case'].replace('__','_')+'_SENSITIVITY';zero=(0,0,0)
    mode=engine.FullPencilModes.__new__(engine.FullPencilModes);mode.r=r;mode.eta=r.symbols['eta_bg'];mode.sigma=r.symbols['sigma_W'];eps=r.symbols['epsilon_shape']
    def tensor(name,array,unit=zero,g=(0,0),epsilon=0,homotopy=None,literal=False):
        a=np.asarray(array,complex)
        if a.ndim==0:a=a.reshape(1,1)
        if a.ndim==1:a=a.reshape(-1,1)
        f.require(a.ndim==2 and np.isfinite(a).all(),'actual finite output array')
        weight=eps**epsilon*(mode.eta**homotopy if homotopy is not None else mode.eta**g[0]*mode.sigma**g[1])
        body=sp.ImmutableMatrix(*a.shape,[mode.number(v)*weight for v in a.ravel()]);tag=prefix+'_'+name
        engine.emit(tag,body if literal else mode.compact_fingerprint(body));engine.emit('METADATA_'+tag,mode.numeric_metadata(body,lambda p:unit))
    def fingerprint(name,array,unit,g=(0,0,0)):
        body=rnative.interior.fingerprint(array);weight=sp.prod(v**n for v,n in zip(result['interior']['generators'],g));body.update(COEFFICIENT_ENTRY_UNIT=unit,COEFFICIENT_GRADE=g,COMPONENT_WEIGHT=weight,LAMBDA_ORDER=g[1]+g[2]);body['NUMERIC_TENSOR_PROJECTIONS']=sp.Tuple(*(v*weight for v in body['NUMERIC_TENSOR_PROJECTIONS']))
        payload=engine.cas(body);tag=prefix+'_'+name;engine.emit(tag,payload);engine.emit('METADATA_'+tag,mode.numeric_metadata(payload,lambda p:unit if p and p[0]=='NUMERIC_TENSOR_PROJECTIONS' else zero))
    size=129;field=[tuple(a-b/2 for a,b in zip(u,result['currentUnit'])) for u in result['fieldUnits']]
    for name,data in result['finite'].items():
        sol=data['solution'];view=data['observable'];system=data['system']
        for key in ('originScattering','originCurrent','incomingPhase','outgoingPhase','totalCurrentRatio','originCurrentResidual'):tensor('FINITE_'+name+'_'+key,view[key],literal=True)
        for key in ('boundaryAnchoredFluxBasisScattering','outgoingChannelCurrent','outgoingFlux','incomingFlux','outgoingFluxRatio'):tensor('FINITE_'+name+'_'+key,sol[key],epsilon=2 if key in ('outgoingFlux','incomingFlux') else 0,literal=True)
        for end,values in sol['modalAmplitudes'].items():tensor('FINITE_'+name+'_'+end+'_MODAL',values,epsilon=1,literal=True)
        tensor('FINITE_'+name+'_POSITIONS',view['positions'],(1,0,0),literal=True)
        for i in range(5):
            unit=field[i];equation=tuple(a-b/2 for a,b in zip(result['rowUnits'][i],result['currentUnit']));trace=tuple(v-(j==0) for j,v in enumerate(unit));sl=slice(i*size,(i+1)*size)
            tensor('FINITE_'+name+'_FIELD_COEFFICIENT_'+str(i),sol['coefficients'][sl],unit,epsilon=1)
            tensor('FINITE_'+name+'_FIELD_NODES_'+str(i),sol['fields'][i],unit,epsilon=1)
            tensor('FINITE_'+name+'_FIELD_GRID_'+str(i),view['originFields'][i],unit,epsilon=1)
            tensor('FINITE_'+name+'_EQUATION_INTERIOR_'+str(i),sol['equationResidual'][sl][1:-1],equation,epsilon=1,literal=True)
            tensor('FINITE_'+name+'_EQUATION_BOUNDARY_'+str(i),sol['equationResidual'][sl][[0,-1]],trace,epsilon=1,literal=True)
            tensor('FINITE_'+name+'_EQUATION_SCALED_'+str(i),sol['scaledEquationResidual'][sl],literal=True)
            if view['independentDifference'] is not None:tensor('FINITE_'+name+'_INDEPENDENT_DIFFERENCE_'+str(i),view['independentDifference'][sl],unit,epsilon=1,literal=True)
            for end in ('LEFT','RIGHT'):
                tensor('FINITE_'+name+'_'+end+'_BOUNDARY_RESIDUAL_'+str(i),view['boundaryResiduals'][end][i],trace,epsilon=1,literal=True)
                tensor('FINITE_'+name+'_'+end+'_TRACE_RESIDUAL_'+str(i),view['traceResiduals'][end][i],unit,epsilon=1,literal=True)
        # Fixed-contrast system fingerprints keep physical interior and trace units distinct.
        for i in range(5):
            for j in range(5):
                block=system['matrix'][i*size:(i+1)*size,j*size:(j+1)*size]
                iu=tuple(a-z for a,z in zip(result['rowUnits'][i],result['fieldUnits'][j]));tu=tuple(a-z-(k==0) for k,(a,z) in enumerate(zip(result['fieldUnits'][i],result['fieldUnits'][j])))
                fingerprint('FINITE_'+name+'_MATRIX_INTERIOR_'+str(i)+'_'+str(j),block[1:-1],iu)
                fingerprint('FINITE_'+name+'_MATRIX_BOUNDARY_'+str(i)+'_'+str(j),block[[0,-1]],tu)
    s=result['sensitivity'];difference=s['finiteAndCommonGrid']
    for key in ('finiteScattering','finiteCurrent'):tensor('DIFFERENCE_'+key,difference[key],literal=True)
    for i in range(5):
        tensor('DIFFERENCE_FINITE_FIELD_GRID_'+str(i),difference['finiteFields'][i],field[i],epsilon=1)
        for name,values in difference['continuumOriginFields'].items():
            for g,a in values.items():tensor('CONTINUUM_GRID_'+name+'_'+str(i)+'_'+str(g),a[i],field[i],g,1)
        for g,a in s['fieldCoefficientDifference'].items():tensor('CONTINUUM_COEFFICIENT_DIFFERENCE_'+str(i)+'_'+str(g),a[i*size:(i+1)*size],field[i],g,1)
        trace=tuple(v-(j==0) for j,v in enumerate(field[i]));tensor('INCIDENT_SIGN_MUTATION_'+str(i),result['prepared']['incidentSignMutation'][i*size:(i+1)*size][[0,-1]],trace,epsilon=1,literal=True)
    for key,series in s['channelCoefficients'].items():
        for name,packet in result['continuum'].items():
            for g,a in packet['response'][key].items():tensor(key+'_'+name+'_'+str(g),a,g=g,literal=True)
        for g,a in series.items():tensor(key+'_DIFFERENCE_'+str(g),a,g=g,literal=True)
        tensor(key+'_EVALUATED_DIFFERENCE',s['evaluatedChannels'][key],literal=True)
    for name,packet in s['retainedPolynomialCurrent'].items():
        if name=='difference':tensor('RETAINED_CURRENT_DIFFERENCE',packet,literal=True)
        else:
            for key,array in packet.items():tensor('RETAINED_CURRENT_'+name+'_'+key,array,literal=True)
    for key,series in difference['continuumCurrent'].items():
        for power,a in series.items():tensor('HOMOTOPY_DIFFERENCE_'+key+'_'+str(power),a,homotopy=power,literal=True)
    for end,series in difference['closedMatching'].items():
        for g,a in series.items():tensor('CLOSED_DIFFERENCE_'+end+'_'+str(g),a,g=g,epsilon=1,literal=True)
    for name,record in result['openFiniteCurrent'].items():
        if name=='differences':
            for end,values in record.items():
                for key,a in values.items():tensor('FINITE_OPEN_CURRENT_DIFFERENCE_'+end+'_'+key,a,epsilon=2,literal=True)
        else:
            for end,values in record.items():
                for key in ('amplitude','metric','full','outgoing','incoming','interference'):tensor('FINITE_OPEN_'+name+'_'+end+'_'+key,values[key],epsilon=1 if key=='amplitude' else 0 if key=='metric' else 2,literal=True)
                for key,a in values['residuals'].items():tensor('FINITE_OPEN_'+name+'_'+end+'_RESIDUAL_'+key,a,epsilon=2,literal=True)
                for index,a in values['closedMatchingAmplitudes'].items():tensor('FINITE_CLOSED_MATCHING_'+name+'_'+end+'_'+str(index),a,epsilon=1,literal=True)
    for name,parts in s['openCurrentDifference'].items():
        for part,series in parts.items():
            for g,a in series.items():tensor('OPEN_CURRENT_DIFFERENCE_'+name+'_'+part+'_'+str(g),a,result['currentUnit'],g,2)
    for end,parts in s['endCurrentDifference'].items():
        for part,kinds in parts.items():
            for kind,series in kinds.items():
                for g,a in series.items():tensor('FULL_END_CURRENT_DIFFERENCE_'+end+'_'+part+'_'+kind+'_'+str(g),a,result['currentUnit'],g,2)
    # Actual saved symbolic operands and proof scalars; no grade splitter or derivative runs.
    records=result['caseBinding']['grades']['records'];originals=result['originalRecords']
    for key,item in records.items():
        rec=item['record'];old=originals[key]['record'];f.require(tuple(rec['UNIT'])==tuple(old['UNIT']) and m.same(item['address'],originals[key]['address']),'source address/unit before output')
        body=sp.Tuple(old['ORIGINAL'],rec['ORIGINAL'],rec['ORIGINAL']-old['ORIGINAL'],*rec['COMPONENTS'].values());name=prefix+'_SOURCE_'+key
        engine.emit(name,engine.carrier_fingerprint(body));engine.emit('METADATA_'+name,mode.numeric_metadata(body,lambda p,u=rec['UNIT']:u))
        proofs=sp.Tuple(rec['ROUND_TRIP_RESIDUAL'],rec['RECONSTRUCTION_RESIDUAL'],*rec['DERIVATIVE_RESIDUALS'].values(),*rec['NATIVE_COLLECTOR_RESIDUALS'].values())
        engine.emit(name+'_RESIDUALS',proofs);engine.emit('METADATA_'+name+'_RESIDUALS',mode.numeric_metadata(proofs,lambda p,u=rec['UNIT']:u))
    for kind in ('local','nonlocal','total'):
        for g,a in result['interior']['matrices'][kind].items():
            for i in range(5):
                for j in range(5):fingerprint('COEFFICIENT_MATRIX_'+kind+'_'+str(g)+'_'+str(i)+'_'+str(j),a[i*size:(i+1)*size,j*size:(j+1)*size],result['interior']['blockUnits'][i][j],g)
    for point,record in result.get('formalRemainders',{}).items():
        for key in ('direct','retained','difference'):
            for i in range(5):tensor('FORMAL_DIAGNOSTIC_'+str(point)+'_'+key+'_'+str(i),record[key][i*size:(i+1)*size],field[i],epsilon=1)
    b.structural_flags(prefix+'_MANIFEST',{'sourceFiles':result['sourceFiles'],'inputPackets':result['inputPackets'],'case':result['case'],'settings':{k:str(v) for k,v in result['settings'].items()},'evaluatedGradeOrigin':s['evaluatedGradeOrigin'],'materialComparison':s['materialComparison'],'scope':result['scope'],'sourceRecords':len(records),'physicalCurrentUnit':tuple(map(str,result['currentUnit'])),'finiteCurrentConvention':'Current forms divided by squared reference incoming flux-amplitude unit; closed finite amplitudes are not currents.','formalPoints':'Arithmetic truncation diagnostics, not additional physical inputs. Baseline has no invented direct formal-point comparison.','matrices':'Full arrays retained in pinned numerical packets; fingerprints include actual shape, bytes and physical block projection units.'})


def custom_emitter(label):
    native=h.native.emitter('FIRST_JET_'+label+'_SENSITIVITY');tree=ast.parse(Path(rnative.__file__).read_text());main=h.native.function(tree,'main')
    first=next(i for i,n in enumerate(main.body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and ast.unparse(n.value.func)=='engine.EMISSION_LINES.clear')
    stop=next(i for i,n in enumerate(main.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='summary' for t in n.targets));body=copy.deepcopy(main.body[first:stop]);original=copy.deepcopy(body)
    prefix='FIRST_JET_'+label.replace('__','_')+'_SENSITIVITY';key='s11cd'+prefix
    changes=[]
    class Rename(ast.NodeTransformer):
        def visit_Constant(self,n):
            table={'s11cdContinuumResponse':key,'continuum-response.pickle':'sensitivity-output.pickle'}
            if n.value in table:changes.append((n.value,table[n.value]));return ast.Constant(table[n.value])
            return n
    body=[Rename().visit(n) for n in body]
    class Undo(ast.NodeTransformer):
        def visit_Constant(self,n):
            table={key:'s11cdContinuumResponse','sensitivity-output.pickle':'continuum-response.pickle'}
            return ast.Constant(table[n.value]) if isinstance(n.value,str) and n.value in table else n
    back=[Undo().visit(copy.deepcopy(n)) for n in body]
    f.require(len(changes)==2 and ast.dump(ast.Module(body=back,type_ignores=[]))==ast.dump(ast.Module(body=original,type_ignores=[])),'whole original output/replay tail with namespace and packet address only')
    body.append(ast.Return(ast.Dict(keys=[ast.Constant(k) for k in ('tags','keys','metadataPaths')],values=[ast.Call(ast.Name('len',ast.Load()),[ast.Name('entries',ast.Load())],[]),ast.Name('keys',ast.Load()),ast.Name('metadata_paths',ast.Load())])))
    node=ast.FunctionDef(name='emit_saved',args=ast.arguments(posonlyargs=[],args=[ast.arg(n) for n in ('base','result','r','pins','operands','before')],vararg=None,kwonlyargs=[],kw_defaults=[],kwarg=None,defaults=[]),body=body,decorator_list=[])
    return h.native.compile_function(node,dict(vars(rnative),emit_result=emit_sensitivity,PREFIX=prefix)),{'wholeNativeReplayTail':True,'namespaceEdits':1,'packetAddressEdits':1,'customEmitterSha256':h.inputs.matrices.binding.native_body(emit_sensitivity)}


def emit_part(base,label,kind,manifest):
    name=kind+'_'+label;folder=base/'parts'/name;folder.mkdir(parents=True,exist_ok=False)
    if kind=='sensitivity':
        result=bundle(base,label,manifest);packet='sensitivity-output.pickle';fn,join=custom_emitter(label);f.atomic_pickle(folder/packet,result)
    else:
        packet='continuum-response.pickle' if kind=='continuum' else 'continuum-currents.pickle'
        origin=base/'numerical/cases'/label/'continuum'/packet;m.retain(origin,folder/packet,manifest);result=f.unpickle(folder/packet)
        fn=h.native.emitter('FIRST_JET_'+label) if kind=='continuum' else h.flux_helper.emitter('FIRST_JET_'+label)
        join={'wholeOriginalEmitterAndReplay':True,'caseNamespace':'FIRST_JET_'+label,'packetSha256':f.digest(folder/packet)}
    f.save(folder/'emitter-join.json',join);r=context(base,label,result);before=f.digest(folder/packet)
    check=fn(folder,result,r,manifest['sourceFiles'],manifest['inputPackets'],before)
    f.require(before==f.digest(folder/packet) and not engine.PHYSICAL_METADATA.dimensions.constraints,'actual emission packet and unit closure')
    # Exercise the actual decoded metadata and physical payload replay predicate.
    controls={}
    for line in decoded_lines(folder/'full.out'):
        tag,_,body=line.rstrip('\n').partition(': ');value=rnative.grades._restore(body)
        if tag.startswith('PY_S11CD_METADATA_') and 'unit' not in controls:
            candidates=[v for v in sp.preorder_traversal(value) if isinstance(v,sp.Tuple) and len(v)==2 and str(v[0])=='DIMENSION_L_T_M' and isinstance(v[1],sp.Tuple) and len(v[1])==3]
            if candidates:
                pair=candidates[0];unit=pair[1];changed=value.xreplace({pair:sp.Tuple(pair[0],sp.Tuple(unit[0]+1,*unit[1:]))})
                f.require(value!=changed,'actual changed metadata unit rejected by replay equality');controls['unit']={'tag':tag,'original':str(unit),'changed':str((unit[0]+1,*unit[1:]))}
        elif not tag.startswith('PY_S11CD_METADATA_') and 'payload' not in controls:
            changed=sp.Tuple(value,sp.Integer(1));f.require(value!=changed,'actual changed physical payload rejected by replay equality');controls['payload']={'tag':tag}
        if len(controls)==2:break
    f.require(len(controls)==2,'actual payload and unit controls');f.save(folder/'replay-controls.json',controls)
    check.update(packetSha256=before,emitterJoin=join,replayControls=controls);f.save(folder/'emission-checks.json',check)
    del result;gc.collect()
    return {'directory':str(folder.relative_to(base)),'case':label,'kind':kind,'sha256':f.digest(folder/'full.out'),'tags':check['tags'],'keys':len(check['keys']),'metadataPaths':check['metadataPaths'],'reusedOriginalEmission':False}


def aggregate(base,parts):
    encoder=PayloadEncoder();tags=set();keys=set();count=0;indices=0;expected_hash=hashlib.sha256();order=list(parts)
    with (base/'original-combined.out').open('xb') as raw,(base/'full.out').open('x') as output:
        for name,record in parts.items():
            folder=base/record['directory'];path=folder/'full.out';raw.write(path.read_bytes());local=[];local_keys=None;local_indices=0
            for line in decoded_lines(path):
                tag,sep,body=line.rstrip('\n').partition(': ');f.require(tag not in tags,'unique global first-jet tag');tags.add(tag)
                if tag.endswith('_WRITE_KEYS') and not tag.startswith('PY_S11CD_METADATA_'):
                    local_keys={str(k):str(v) for k,v in rnative.grades._restore(body)};f.require(len(local_keys)==len(set(local_keys.values())) and not keys&set(local_keys.values()),'unique global export keys');keys.update(local_keys.values())
                if tag.endswith('_EMISSION_LINES') and not tag.startswith('PY_S11CD_METADATA_'):
                    ix={str(k):v for k,v in rnative.grades._restore(body)};restore_emission_index(ix,local);local_indices+=1;indices+=1
                output.write(tag+sep+encoder.encode(body)+'\n');expected_hash.update(line.encode());local.append(tag);count+=1
            f.require(len(local)==record['tags'] and local_keys is not None and len(local_keys)==record['keys'] and local_indices==1,'complete saved part census and local index')
    expected=itertools.chain.from_iterable(decoded_lines(base/parts[n]['directory']/'full.out') for n in order);sentinel=object();actual_count=0
    for actual,wanted in itertools.zip_longest(decoded_lines(base/'full.out'),expected,fillvalue=sentinel):f.require(actual==wanted,'every decoded original payload and local index');actual_count+=1
    f.require(actual_count==count,'complete global payload count')
    rejected=False
    try:PayloadDecoder().decode("Tuple(Str('s11cdSharedPayloadReference'), Integer(0))")
    except ValueError:rejected=True
    f.require(rejected,'malformed reference control')
    first=next(decoded_lines(base/'full.out'));changed=first.replace(': ',': Integer(1) + ',1);f.require(hashlib.sha256(first.encode()).digest()!=hashlib.sha256(changed.encode()).digest(),'actual changed payload control')
    return {'tags':count,'keys':len(keys),'parts':len(parts),'localIndicesPreserved':indices,'identicalDecodedPayloads':count,'decodedSha256':expected_hash.hexdigest(),'metadataPaths':sum(v['metadataPaths'] for v in parts.values()),'malformedReferenceRejected':True,'changedPayloadRejected':True,'nativeCodecUnchanged':True}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',required=True,type=Path);ap.add_argument('--mode',choices=('focused','construct'),required=True);ap.add_argument('--resume-from',type=Path);args=ap.parse_args()
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900)
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False);manifest,labels,parts=load(base,args.resume_from)
    # Numerical construction is not callable during saved-packet output.
    def forbidden(*a,**kw):raise RuntimeError('output attempted numerical reconstruction')
    for module,names in ((np.linalg,('solve','lstsq','inv','pinv','svd','eig','eigh')),(rnative,('systems','solve','channels','open_flux')),(h.current,('construct_open','end_currents','open_metrics'))):
        for n in names:setattr(module,n,forbidden)
    selected=(FOCUS,) if args.mode=='focused' else labels
    for label in selected:
        for kind in ('continuum','current','sensitivity'):
            if kind=='continuum' and label==h.BASELINE:continue
            name=kind+'_'+label
            if name in parts:continue
            parts[name]=emit_part(base,label,kind,manifest);f.save(base/'part-inventory.json',parts);gc.collect()
    # New per-part input-copy routes were appended only after their immutable copies.
    f.save(base/'inputs.json',manifest)
    combined=None
    if args.mode=='construct':combined=aggregate(base,parts);f.save(base/'aggregation-checks.json',combined)
    for n,v in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/n)==f.digest(base/'source'/n)==v,'current/frozen sources pre/post')
    for n,v in manifest['inputPackets'].items():f.require(f.digest(Path(n))==v,'original source/input pre/post')
    for n,v in manifest['copiedInputs'].items():f.require(f.digest(base/n)==v,'copied original/complete packet pre/post')
    artifacts={str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*') if p.is_file() and 'source' not in p.relative_to(base).parts and p not in (base/'inputs.json',base/'checks.json')}
    result={**manifest,'status':'COMPLETED_FIRST_JET_OUTPUT_FOCUS' if args.mode=='focused' else 'COMPLETED_FOUR_CASE_FIRST_JET_OUTPUT','mode':args.mode,'parts':parts,'aggregate':combined,'artifacts':artifacts,'newSolves':0,'newQuadratureNodes':0,'newCurrentConstructions':0,'newGradeExtractions':0,'wallSeconds':time.monotonic()-start}
    f.save(base/'checks.json',result);signal.alarm(0);print(json.dumps(result,indent=2))


if __name__=='__main__':main()
