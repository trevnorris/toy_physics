#!/usr/bin/env python3
"""Saved material-coordinate transcripts, with separate emission/replay phases."""
import argparse,ast,copy,gc,json,resource,shutil,signal,time,types
from pathlib import Path
import numpy as np
import sympy as sp
import S11c_d_remaining_case_coordinate_response as h
import S11c_d_remaining_case_first_jet_output as original
import S11c_d_remaining_case_first_jet_output_finish as finished
import S11c_d_remaining_case_first_jet_output_recover as recovered
f,m,b=h.f,h.m,h.b;engine=f.engine;rnative=h.native.response
CP=f.M/'S11c_d_remaining_case_coordinate_response_checkpoint.json'
FCP=f.M/'S11c_d_remaining_case_coordinate_output_focused.json'
PLAN=f.M/'S11c_d_remaining_case_coordinate_output_plan.md'
WIRING=f.M/'S11c_d_remaining_case_coordinate_output_wiring.json'
ADAPTER=f.M/'S11c_d_remaining_case_coordinate_output_emitter_adapter.json'
SCOPE=('Own-case material-coordinate controls constructed in common Eulerian coordinates before boundary replacement and response solves. '
 'Four genuinely new finite controls, three new continuum controls and the reused historical baseline continuum. '
 'Finite current is open-only 4 by 4 with three separate closed matching amplitudes; continuum retains seven directions and closed/cross terms. '
 'Coordinate consistency does not resolve tiny reflection/loss. Positive regulator, approximate boundaries and omitted parent pure-second-order scope remain. '
 'No post-hoc S conjugation, unlike-density subtraction, isolated advection or c2 closure claim.')


def load(base):
    wiring=json.loads(WIRING.read_text());f.require(wiring['helperSha256']==f.digest(Path(__file__)) and wiring['emitterAdapterSha256']==f.digest(ADAPTER),'exact reviewed output helper/wiring')
    cp,origin=h.native.checked_checkpoint(CP,'ACCEPTED_CASE_MATERIAL_NUMERICAL_RESPONSES');checks=json.loads((origin/'checks.json').read_text())
    v=cp['validation'];vr=Path(v['runDirectory']);inv=json.loads((vr/'validate.invocation.json').read_text());guard=json.loads((vr/'resource-guard/outcome.json').read_text())
    f.require(inv['exitCode']==guard['exitCode']==guard['childOutcome']['exitCode']==0 and inv['stderrBytes']==guard['stderrBytes']==0 and guard['limitsVerified'] and guard['childOutcome']['guardReason'] is None,'clean independent numerical validation')
    f.require(f.digest(vr/'checks.json')==v['checksSha256'] and (vr/'checks.json').read_bytes()==(vr/'validate.stdout').read_bytes(),'actual validated final checks/stdout')
    pins=dict(checks['sourceFiles'])
    for p in (Path(__file__),PLAN,WIRING,ADAPTER,CP,Path(original.__file__),Path(finished.__file__),Path(recovered.__file__)):
        n=str(p.resolve().relative_to(f.ROOT));f.require(n not in pins or pins[n]==f.digest(p),'unchanged inherited output helper');pins[n]=f.digest(p)
    manifest={'runDirectory':str(base),'sourceFiles':pins,'inputPackets':dict(checks['inputPackets']),'copiedInputs':{},'input':checks['input'],'settings':checks['settings'],'scope':SCOPE,
        'numericalOrigin':str(origin),'numericalChecksSha256':cp['checksSha256'],'acceptedValidation':v['checksSha256'],'materialUnitSupplement':checks['materialUnitSupplement']}
    for p in (origin/'checks.json',origin/'inputs.json',vr/'checks.json',vr/'validate.py',vr/'validator-join.json',vr/'current-validator-join.json'):manifest['inputPackets'][str(p)]=f.digest(p)
    for n,v in checks['artifacts'].items():m.retain(origin/n,base/'numerical'/n,manifest,v['sha256'])
    m.retain(origin/'inputs.json',base/'numerical/inputs.json',manifest)
    folder=base/'parts'/('historical_'+h.BASELINE);folder.mkdir(parents=True)
    for name,target in [('full.out','full.out'),('checks.json','emission-checks.json')]:m.retain(origin/'cases'/h.BASELINE/'continuum'/name,folder/target,manifest)
    old=json.loads((folder/'emission-checks.json').read_text());parts={'historical_'+h.BASELINE:{'directory':str(folder.relative_to(base)),'kind':'historical','case':h.BASELINE,'sha256':f.digest(folder/'full.out'),'tags':old['tagCount'],'keys':old['writeKeys'],'metadataPaths':old['metadataPaths'],'reusedOriginalEmission':True}}
    for n,sha in pins.items():
        p=base/'source'/n;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/n,p);f.require(f.digest(p)==sha,'frozen output source')
    f.save(base/'inputs.json',manifest);f.save(base/'part-inventory.json',parts)
    return manifest,tuple(checks['cases']),parts


def context(base,label,packet,folder):
    source_path=base/'numerical/boundary-inputs/bindings/sources/accepted-cases'/label/'reduced-action.pickle'
    source=f.unpickle(source_path);r,d=f.prior.domain.momentum.source.native.source.restore_context(source)
    d.__dict__.update(packet['dimensionState']);f.require(engine.PHYSICAL_METADATA.dimensions is d,'restored actual metadata context')
    supplement=base/'numerical/boundary-unit-supplement/unit-overlay';pairs=f.unpickle(supplement/'unit-source-pairs.pickle');view=f.unpickle(supplement/'dimension-state-view.pickle')
    raw=dict(d.known);joined=pairs['joinedKnown']
    for atom,unit in joined.items():
        f.require(tuple(view['known'][atom])==tuple(unit)==(-1,0,0),'exact accepted current-unit supplement')
        if atom in raw:f.require(m.same(raw[atom],unit) or m.same(raw[atom],pairs['rawKnown'][atom]),'raw current declaration joins original or accepted view')
    d.known=dict(raw);d.known.update(joined)
    proof={'source':str(source_path),'sourceSha256':f.digest(source_path),'packetDimensionState':packet['dimensionState'],'rawKnown':{k:raw.get(k) for k in joined},'joinedKnown':joined,'sourcePairsSha256':f.digest(supplement/'unit-source-pairs.pickle'),'savedViewSha256':f.digest(supplement/'dimension-state-view.pickle'),'originalPacketsUnchanged':True}
    path=folder/'context-unit-join.pickle'
    if path.exists():f.require(m.same(f.unpickle(path),proof),'reuse exact saved context/unit evidence')
    else:f.atomic_pickle(path,proof)
    return r


def bundle(base,label,manifest):
    n=base/'numerical';q=n/'cases'/label;old=n/'unchanged-response/cases'/label
    value=f.unpickle(q/'continuum/continuum-response.pickle');co=f.unpickle(n/'interiors/cases'/label/'interior-matrices.pickle')
    ends=f.unpickle(n/'boundary-inputs/material-cases'/label/'case-material-boundary.pickle');families={}
    for end,owner in ends['materialFamilyRoutes'].items():
        folder=n/'boundary-inputs/material-families'/owner;inv=json.loads((folder/'checks.json').read_text());source_end=inv['sourceEnd']
        families[end]={'owner':owner,'finite':f.unpickle(folder/'finite/finite-material-boundary.pickle'),
            'continuum':f.unpickle(folder/'continuum'/(source_end.lower()+'-material-boundary.pickle')),
            'units':f.unpickle(folder/'array-units.pickle'),'sourceCurrentUnits':f.unpickle(folder/'source-current-units.pickle'),
            'input':f.unpickle(folder/'input.pickle'),'checks':inv}
        f.require(m.same(families[end]['finite']['commonEulerian'],ends['finite'][end]) and m.same(families[end]['continuum']['commonEulerian'],ends['ends'][end]),'actual own constructed boundary maps')
    result={'case':label,'finite':{},'continuum':{},'comparison':f.unpickle(q/'coordinate-comparisons.pickle'),
        'currentDifferences':f.unpickle(q/'coordinate-current-comparisons.pickle'),'current':f.unpickle(q/'continuum/continuum-currents.pickle'),
        'openFiniteCurrent':f.unpickle(q/'finite/open-end-currents.pickle'),'commonCoordinateJoins':f.unpickle(q/'common-coordinate-joins.pickle'),
        'families':families,'chartState':f.unpickle(n/'boundary-inputs/bindings/chart-state.pickle'),'interior':co,'forcingControls':f.unpickle(n/'preparation'/label/'forcing-controls.pickle'),
        'gradeOriginNumerical':{str(g):float(co['gradeOrigin'][g]) for g in co['generators'][1:]},
        'fieldUnits':value['fieldUnits'],'rowUnits':value['rowUnits'],'currentUnit':value['currentUnit'],'dimensionState':value['dimensionState'],
        'sourceFiles':manifest['sourceFiles'],'inputPackets':manifest['inputPackets'],'settings':value['settings'],
        'materialUnitSupplement':manifest['materialUnitSupplement'],'scope':SCOPE}
    for name,root in [('MATERIAL',q),('EULERIAN',old)]:
        result['finite'][name]={k:f.unpickle(root/'finite'/(fn+'.pickle')) for k,fn in [('system','finite-system'),('solution','finite-solution'),('observable','observable')]}
        result['continuum'][name]=f.unpickle(root/'continuum/continuum-response.pickle')
    if label!=h.BASELINE:result['formalRemainders']=f.unpickle(q/'continuum/formal-remainders.pickle')
    result['sourceCensus']=json.loads((q/'response-disposition.json').read_text())
    return result


def boundary_arrays(result):
    """Route every accepted map-unit record to its actual saved map array."""
    for end,family in result['families'].items():
        for address,record in family['units'].items():
            kind,*path=address
            if kind=='finite':
                key='current' if path[0]=='openCurrent' else path[0]
                for frame in ('material','commonEulerian'):
                    a=family['finite'][frame][key];yield end+'_'+frame+'_'+str(address),a,record['units'],(0,0),0
                if key in family['finite']['differences']:yield end+'_RESIDUAL_'+str(address),family['finite']['differences'][key],record['units'],(0,0),0
            elif path[0]!='cluster':
                if path[0]=='current':key,name,g=path;arrays={frame:family['continuum'][frame]['currents'][name][g] for frame in ('material','commonEulerian')};epsilon=2
                else:
                    key,g=path;arrays={'commonEulerian':family['continuum']['commonEulerian'][key][g]};epsilon=0
                    if key in family['continuum']['material']['boundary']:arrays['material']=family['continuum']['material']['boundary'][key][g]
                for frame,a in arrays.items():yield end+'_'+frame+'_'+str(address),a,record['units'],g,epsilon
            else:
                # Full cluster subspaces are already pinned; export all actual
                # common-frame arrays with their accepted per-entry unit maps.
                _,index,key,*degree=path;cluster=next(c for c in family['continuum']['commonEulerian']['clusters'] if c['info']['INDEX']==index)
                a=cluster['diagnostics']['BASE_EQUATION'] if key=='baseEquation' else cluster[key][degree[0]]
                yield end+'_COMMON_'+str(address),a,record['units'],degree[0] if degree else (0,0),0


def emit_boundary_maps(result,tensor):
    for name,array,units,grade,epsilon in boundary_arrays(result):
        tensor('BOUNDARY_'+name,array,lambda i,j,u=units:u[i][j],grade,epsilon)
    for end,family in result['families'].items():
        raw=result['commonCoordinateJoins'][end]
        tensor('COMMON_'+end+'_NORMAL_MOMENTUM_RESIDUAL',raw['normalMomentum'],(-1,0,0),literal=True)
        for name,key in [('rightBasis','right'),('incomingValues','incomingValues'),('traceMap','traceMap'),('finiteCurrent','openCurrent')]:
            units=family['units'][('finite',key)]['units'];tensor('COMMON_'+end+'_'+name+'_RESIDUAL',raw[name],lambda i,j,u=units:u[i][j],literal=True)
        for direction,values in family['finite']['phases'].items():
            for name,array in values.items():tensor('PHASE_'+end+'_'+direction+'_'+name,array,literal=True)
        # These physical omissions have distinct dimensions; retain each matrix
        # row/column unit from its actual accepted derivative/current operands.
        units=family['units'][('finite','derivative')]['units']
        tensor('CONTROL_'+end+'_NORMAL_COVECTOR',family['finite']['controls']['omittedDerivativeCovector'],lambda i,j,u=units:u[i][j],literal=True)
        for direction,a in family['finite']['controls']['omittedShearPhase'].items():tensor('CONTROL_'+end+'_SHEAR_PHASE_'+direction,a,literal=True)
        for name in ('omittedMeasure','currentCoefficient'):tensor('CONTROL_'+end+'_'+name,family['finite']['controls'][name],literal=True)


def emit_coordinate(result,r):
    prefix='MATERIAL_COORDINATE_'+result['case'].replace('__','_')+'_CONTROL';zero=(0,0,0)
    mode=engine.FullPencilModes.__new__(engine.FullPencilModes);mode.r=r;mode.eta=r.symbols['eta_bg'];mode.sigma=r.symbols['sigma_W'];eps=r.symbols['epsilon_shape']
    def tensor(name,array,unit=zero,g=(0,0),epsilon=0,homotopy=None,literal=False):
        a=np.asarray(array,complex)
        if a.ndim==0:a=a.reshape(1,1)
        if a.ndim==1:a=a.reshape(-1,1)
        f.require(a.ndim==2 and np.isfinite(a).all(),'actual finite output array')
        weight=eps**epsilon*(mode.eta**homotopy if homotopy is not None else mode.eta**g[0]*mode.sigma**g[1])
        body=sp.ImmutableMatrix(*a.shape,[mode.number(v)*weight for v in a.ravel()]);tag=prefix+'_'+name
        engine.emit(tag,body if literal else mode.compact_fingerprint(body));engine.emit('METADATA_'+tag,mode.numeric_metadata(body,lambda p:unit(p[0]//a.shape[1],p[0]%a.shape[1]) if callable(unit) else unit))
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
            unit=field[i];equation=tuple(a-b/2 for a,b in zip(result['rowUnits'][i],result['currentUnit']));trace=tuple(v-(1 if j==0 else 0) for j,v in enumerate(unit));sl=slice(i*size,(i+1)*size)
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
                iu=tuple(a-z for a,z in zip(result['rowUnits'][i],result['fieldUnits'][j]));tu=tuple(a-z-(1 if k==0 else 0) for k,(a,z) in enumerate(zip(result['fieldUnits'][i],result['fieldUnits'][j])))
                fingerprint('FINITE_'+name+'_MATRIX_INTERIOR_'+str(i)+'_'+str(j),block[1:-1],iu)
                fingerprint('FINITE_'+name+'_MATRIX_BOUNDARY_'+str(i)+'_'+str(j),block[[0,-1]],tu)
    s=result['comparison'];difference=s
    for key in ('finiteScattering','finiteCurrent'):tensor('DIFFERENCE_'+key,difference[key],literal=True)
    for i in range(5):
        tensor('DIFFERENCE_FINITE_FIELD_GRID_'+str(i),difference['finiteFields'][i],field[i],epsilon=1)
        for name,values in difference['continuumFields'].items():
            for g,a in values.items():tensor('CONTINUUM_GRID_'+name+'_'+str(i)+'_'+str(g),a[i],field[i],g,1)
        for g,a in s['continuumFieldCoefficients'].items():tensor('CONTINUUM_COEFFICIENT_DIFFERENCE_'+str(i)+'_'+str(g),a[i*size:(i+1)*size],field[i],g,1)
        trace=tuple(v-(1 if j==0 else 0) for j,v in enumerate(field[i]));tensor('INCIDENT_SIGN_MUTATION_'+str(i),result['forcingControls']['incidentSign'][i*size:(i+1)*size][[0,-1]],trace,epsilon=1,literal=True)
    for key,series in s['continuumChannels'].items():
        for name,packet in result['continuum'].items():
            for g,a in packet['response'][key].items():tensor(key+'_'+name+'_'+str(g),a,g=g,literal=True)
        for g,a in series.items():tensor(key+'_DIFFERENCE_'+str(g),a,g=g,literal=True)
        tensor(key+'_EVALUATED_DIFFERENCE',s['evaluatedChannels'][key],literal=True)
    for name,packet in s['retainedPolynomialCurrents'].items():
        if name=='difference':tensor('RETAINED_CURRENT_DIFFERENCE',packet,literal=True)
        else:
            for key,array in packet.items():tensor('RETAINED_CURRENT_'+name+'_'+key,array,literal=True)
    for end,series in difference['closedMatching'].items():
        for g,a in series.items():tensor('CLOSED_DIFFERENCE_'+end+'_'+str(g),a,g=g,epsilon=1,literal=True)
    for i in range(5):
        trace=tuple(v-(1 if j==0 else 0) for j,v in enumerate(field[i]))
        tensor('OMITTED_INCOMING_COVECTOR_'+str(i),result['forcingControls']['omittedIncomingCovector'][i*size:(i+1)*size][[0,-1]],trace,epsilon=1,literal=True)
    for name,record in result['openFiniteCurrent'].items():
        if name=='differences':
            for end,values in record.items():
                for key,a in values.items():tensor('FINITE_OPEN_CURRENT_DIFFERENCE_'+end+'_'+key,a,epsilon=2,literal=True)
        else:
            for end,values in record.items():
                for key in ('amplitude','metric','full','outgoing','incoming','interference'):tensor('FINITE_OPEN_'+name+'_'+end+'_'+key,values[key],epsilon=1 if key=='amplitude' else 0 if key=='metric' else 2,literal=True)
                for key,a in values['residuals'].items():tensor('FINITE_OPEN_'+name+'_'+end+'_RESIDUAL_'+key,a,epsilon=2,literal=True)
                for index,a in values['closedMatchingAmplitudes'].items():tensor('FINITE_CLOSED_MATCHING_'+name+'_'+end+'_'+str(index),a,epsilon=1,literal=True)
    for name,parts in result['currentDifferences']['open'].items():
        for part,series in parts.items():
            for g,a in series.items():tensor('OPEN_CURRENT_DIFFERENCE_'+name+'_'+part+'_'+str(g),a,result['currentUnit'],g,2)
    for end,parts in result['currentDifferences']['ends'].items():
        for part,kinds in parts.items():
            for kind,series in kinds.items():
                for g,a in series.items():tensor('FULL_END_CURRENT_DIFFERENCE_'+end+'_'+part+'_'+kind+'_'+str(g),a,result['currentUnit'],g,2)
    emit_boundary_maps(result,tensor)
    for kind in ('local','nonlocal','total'):
        for g,a in result['interior']['matrices'][kind].items():
            for i in range(5):
                for j in range(5):fingerprint('COEFFICIENT_MATRIX_'+kind+'_'+str(g)+'_'+str(i)+'_'+str(j),a[i*size:(i+1)*size,j*size:(j+1)*size],result['interior']['blockUnits'][i][j],g)
    for point,record in result.get('formalRemainders',{}).items():
        for key in ('direct','retained','difference'):
            for i in range(5):tensor('FORMAL_DIAGNOSTIC_'+str(point)+'_'+key+'_'+str(i),record[key][i*size:(i+1)*size],field[i],epsilon=1)
    b.structural_flags(prefix+'_MANIFEST',{'sourceFiles':result['sourceFiles'],'inputPackets':result['inputPackets'],'case':result['case'],'settings':{k:str(v) for k,v in result['settings'].items()},'evaluatedGradeOrigin':result['gradeOriginNumerical'],'scope':result['scope'],'sourceRows':result['sourceCensus']['rows'],'nativeTerms':result['sourceCensus']['terms'],'sourceAmplitudes':result['sourceCensus']['sources'],'physicalCurrentUnit':tuple(map(str,result['currentUnit'])),'materialUnitSupplementRequired':True,'finiteCurrentConvention':'Open-current forms divided by squared reference incoming flux-amplitude unit; closed matching amplitudes are separate.','formalPoints':'Saved arithmetic truncation diagnostics; none invented for historical baseline.','coordinateRoute':'Actual material then common Eulerian maps precede response construction; every full input and source proof is pinned.','matrices':'Full arrays remain in pinned numerical packets; physical entry units and shape/byte fingerprints retained.'})


def coordinate_drivers(label,folder):
    main=h.native.function(ast.parse(Path(rnative.__file__).read_text()),'main')
    first=next(i for i,n in enumerate(main.body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and ast.unparse(n.value.func)=='engine.EMISSION_LINES.clear')
    stop=next(i for i,n in enumerate(main.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='entries' for t in n.targets))
    prefix='MATERIAL_COORDINATE_'+label.replace('__','_')+'_CONTROL';key='s11cd'+prefix
    original_body=copy.deepcopy(main.body[first:stop]);body=copy.deepcopy(original_body)
    class Rename(ast.NodeTransformer):
        def __init__(self,table):self.table=table;self.hits=[]
        def visit_Constant(self,n):
            if isinstance(n.value,str) and n.value in self.table:self.hits.append(n.value);return ast.copy_location(ast.Constant(self.table[n.value]),n)
            return n
    ren=Rename({'s11cdContinuumResponse':key});body=[ren.visit(n) for n in body]
    back=[Rename({key:'s11cdContinuumResponse'}).visit(copy.deepcopy(n)) for n in body]
    f.require(len(ren.hits)==1 and ast.dump(ast.Module(body=back,type_ignores=[]))==ast.dump(ast.Module(body=original_body,type_ignores=[])),'whole native emission prefix namespace join')
    body.extend(ast.parse('return dict(keys=keys,index=index)').body)
    node=ast.FunctionDef(name='emit_only',args=ast.arguments(posonlyargs=[],args=[ast.arg(n) for n in ('base','result','r')],vararg=None,kwonlyargs=[],kw_defaults=[],kwarg=None,defaults=[]),body=body,decorator_list=[])
    emit=h.native.compile_function(node,dict(vars(rnative),emit_result=emit_coordinate,PREFIX=prefix))
    # Existing saved-stream replay suffix: only its custom tag and packet names
    # change. All native metadata/key/index/hash guards and observers remain.
    text=Path(finished.__file__).read_text();fn=h.native.function(ast.parse(text),'validation_tail');source=copy.deepcopy(fn)
    table={'FIRST_JET_':'MATERIAL_COORDINATE_','_SENSITIVITY':'_CONTROL','sensitivity-output.pickle':'coordinate-output.pickle'}
    rename=Rename(table);adapted=rename.visit(copy.deepcopy(source));back=Rename({v:k for k,v in table.items()}).visit(copy.deepcopy(adapted))
    f.require(set(rename.hits)==set(table) and ast.dump(back)==ast.dump(source),'whole saved-stream validation factory reverse AST')
    proxy=types.SimpleNamespace(**vars(original));proxy.emit_sensitivity=emit_coordinate;proxy.__file__=__file__
    factory=h.native.compile_function(adapted,dict(vars(finished),h=proxy));replay,join=factory(label,folder)
    join.update(wholeNativeEmissionPrefix=True,emissionNamespaceEdits=1,wholeValidationFactoryReverseAST=True,validationFactoryEdits=rename.hits,originalFactorySha256=f.digest(Path(finished.__file__)))
    return emit,replay,join


def native_drivers(label,kind,folder):
    module=rnative if kind=='continuum' else h.current
    prefix='MATERIAL_COORDINATE_'+label.replace('__','_')+'_'+kind.upper();key='s11cd'+prefix
    packet='continuum-response.pickle' if kind=='continuum' else 'continuum-currents.pickle'
    oldkey='s11cdContinuumResponse' if kind=='continuum' else 's11cdContinuumCurrents'
    main=h.native.function(ast.parse(Path(module.__file__).read_text()),'main')
    first=next(i for i,n in enumerate(main.body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and ast.unparse(n.value.func)=='engine.EMISSION_LINES.clear')
    middle=next(i for i,n in enumerate(main.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='entries' for t in n.targets))
    stop=next(i for i,n in enumerate(main.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='summary' for t in n.targets))
    original_emit=copy.deepcopy(main.body[first:middle]);body=copy.deepcopy(original_emit)
    class Name(ast.NodeTransformer):
        def __init__(self,a,z):self.a=a;self.z=z;self.count=0
        def visit_Constant(self,n):
            if n.value==self.a:self.count+=1;return ast.copy_location(ast.Constant(self.z),n)
            return n
    rename=Name(oldkey,key);body=[rename.visit(n) for n in body];back=[Name(key,oldkey).visit(copy.deepcopy(n)) for n in body]
    f.require(rename.count==1 and ast.dump(ast.Module(body=back,type_ignores=[]))==ast.dump(ast.Module(body=original_emit,type_ignores=[])),'whole native emission prefix with one key namespace')
    namespace=dict(vars(module),PREFIX=prefix);namespace['emit_result']=types.FunctionType(module.emit_result.__code__,namespace,module.emit_result.__name__)
    args=lambda names:ast.arguments(posonlyargs=[],args=[ast.arg(n) for n in names],vararg=None,kwonlyargs=[],kw_defaults=[],kwarg=None,defaults=[])
    body.extend(ast.parse('return dict(keys=keys,index=index)').body)
    emit=h.native.compile_function(ast.FunctionDef(name='emit_only',args=args(('base','result','r')),body=body,decorator_list=[]),namespace)
    original_replay=copy.deepcopy(main.body[middle:stop]);replay_body=copy.deepcopy(original_replay)
    replay_body.insert(2,ast.parse('keys,index=saved_keys_index(entries,PREFIX)').body[0])
    back=copy.deepcopy(replay_body);del back[2]
    f.require(ast.dump(ast.Module(body=back,type_ignores=[]))==ast.dump(ast.Module(body=original_replay,type_ignores=[])),'whole native replay suffix with saved keys/index only')
    metric='metadata_paths' if kind=='continuum' else 'paths'
    replay_body.extend(ast.parse('return dict(tags=len(entries),keys=keys,metadataPaths='+metric+')').body)
    replay=h.native.compile_function(ast.FunctionDef(name='replay_saved',args=args(('base','result','r','pins','operands','before')),body=replay_body,decorator_list=[]),dict(namespace,saved_keys_index=finished.saved_keys_index))
    return emit,replay,{'wholeNativeEmissionPrefix':True,'wholeNativeReplaySuffix':True,'emitterBytecodeUnchanged':True,'keyNamespaceEdits':1,'savedKeysAndIndexOnly':True,'prefix':prefix,'packet':packet,'sourceSha256':f.digest(Path(module.__file__))}


def prohibit():
    h.prohibit_completed_work()
    def forbidden(*a,**kw):raise RuntimeError('scientific reconstruction prohibited during saved material output')
    for module in (np.linalg,h.native.la):
        for name in ('solve','inv','pinv','lstsq','svd','eig','eigh','eigvals','eigvalsh','matrix_rank','lu_factor','lu_solve','solve_sylvester'):
            if hasattr(module,name):setattr(module,name,forbidden)
    for module,names in ((h,('load','construct','comparison','wiring')),(rnative,('systems','solve','channels','open_flux')),
        (h.current,('selectors','open_metrics','construct_open','end_currents','bulk_domains')),(h.inputs,('load','prepare','baseline_views'))):
        for name in names:setattr(module,name,forbidden)


def sample_emitter(base,label,kind,result,folder):
    prefix='MATERIAL_COORDINATE_'+label
    emit=(emit_coordinate if kind=='coordinate' else types.FunctionType((rnative.emit_result if kind=='continuum' else h.current.emit_result).__code__,
        dict(vars(rnative if kind=='continuum' else h.current),PREFIX=prefix.replace('__','_')+'_'+kind.upper()),(rnative.emit_result if kind=='continuum' else h.current.emit_result).__name__))
    r=context(base,label,result,folder);captured=[];native_emit=engine.emit
    class Captured(Exception):pass
    def first(name,value):
        captured.append((name,engine.cas(value)))
        if len(captured)==2:raise Captured()
    engine.emit=first
    try:
        try:emit(result,r)
        except Captured:pass
    finally:engine.emit=native_emit
    f.atomic_pickle(folder/'first-payload-and-metadata.pickle',captured)
    f.require(len(captured)==2 and captured[1][0].startswith('METADATA_'),'actual native first payload and metadata')
    value=captured[1][1];units=[v for v in sp.preorder_traversal(value) if isinstance(v,sp.Tuple) and len(v)==2 and str(v[0])=='DIMENSION_L_T_M' and isinstance(v[1],sp.Tuple) and len(v[1])==3]
    f.require(units and all(not x.free_symbols for pair in units for x in pair[1]),'resolved actual output unit sample')
    pair=units[0];changed=value.xreplace({pair:sp.Tuple(pair[0],sp.Tuple(pair[1][0]+1,*pair[1][1:]))});mutated=sp.Tuple(captured[0][1],sp.Integer(1))
    f.atomic_pickle(folder/'sample-mutations.pickle',{'originalMetadata':value,'changedMetadata':changed,'originalPayload':captured[0][1],'changedPayload':mutated})
    f.require(value!=changed and captured[0][1]!=mutated and not engine.PHYSICAL_METADATA.dimensions.constraints,'actual unit/payload sample controls')
    return {'kind':kind,'firstTags':[x[0] for x in captured],'changedUnitRejected':True,'changedPayloadRejected':True,'instrumentOnly':True,'physicalTranscriptEmitted':False}


def focused(base,manifest,labels):
    results={}
    for label in labels:
        folder=base/'bundles'/label;folder.mkdir(parents=True)
        value=bundle(base,label,manifest);packet=folder/'coordinate-output.pickle';f.atomic_pickle(packet,value)
        count=0
        for name,a,units,g,eps in boundary_arrays(value):
            array=np.asarray(a);f.require(array.ndim==2 and array.shape==(len(units),len(units[0])) and np.isfinite(array).all(),'all saved boundary output arrays and entry units')
            f.require(all(len(unit)==3 and all(not sp.sympify(v).free_symbols for v in unit) for row in units for unit in row),'full boundary output unit resolution');count+=1
        for record in value['finite'].values():
            f.require(record['system']['matrix'].shape==(645,645) and record['system']['rhs'].shape==(645,4) and record['solution']['coefficients'].shape==(645,4),'full actual finite output systems')
        for atom in value['chartState']['values']['g']['coordinates']:
            f.require(tuple(value['dimensionState']['known'][atom])==(1,0,0),'material and Eulerian coordinate units agree')
        _,_,join=coordinate_drivers(label,folder);f.save(folder/'emitter-wiring.json',join)
        samples={}
        for kind in ('coordinate','current','continuum'):
            if kind=='continuum' and label==h.BASELINE:continue
            sample=folder/('sample-'+kind);sample.mkdir()
            if kind!='coordinate':
                _,_,native_join=native_drivers(label,kind,sample);f.save(sample/'native-driver-join.json',native_join)
            result=value if kind=='coordinate' else value['current'] if kind=='current' else value['continuum']['MATERIAL']
            samples[kind]=sample_emitter(base,label,kind,result,sample)
        results[label]={'bundleSha256':f.digest(packet),'bundleBytes':packet.stat().st_size,'boundaryArrayViews':count,'samples':samples,'fullFiniteSystems':2,'independentGrades':True,'sourceCensus':{k:value['sourceCensus'][k] for k in ('rows','terms','sources')},'newScientificWork':False}
        f.save(base/'output-case-inventory.json',results);del value;gc.collect()
    return results


def hashes(base,manifest):
    for n,sha in manifest['sourceFiles'].items():f.require(f.digest(f.ROOT/n)==f.digest(base/'source'/n)==sha,'current/frozen output source')
    for n,sha in manifest['inputPackets'].items():f.require(f.digest(Path(n))==sha,'original input pre/post')
    for n,sha in manifest['copiedInputs'].items():f.require(f.digest(base/n)==sha,'completed copy pre/post')


def resume(base,focus):
    cp=json.loads(FCP.read_text());old=json.loads((focus/'checks.json').read_text())
    f.require(cp['status']=='ACCEPTED_CASE_MATERIAL_OUTPUT_INPUTS' and cp['runDirectory']==str(focus) and cp['checksSha256']==f.digest(focus/'checks.json'),'accepted output input focus')
    f.require(old['status']=='COMPLETED_CASE_MATERIAL_OUTPUT_INPUTS','complete focus')
    manifest={k:old[k] for k in ('sourceFiles','inputPackets','input','settings','scope','numericalOrigin','numericalChecksSha256','acceptedValidation','materialUnitSupplement')}
    manifest.update(runDirectory=str(base),copiedInputs={},completedOutputInputReuse={'directory':str(focus),'checksSha256':f.digest(focus/'checks.json'),'artifacts':len(old['artifacts'])})
    for n,v in old['artifacts'].items():
        target='accepted-focus-part-inventory.json' if n=='part-inventory.json' else n
        m.retain(focus/n,base/target,manifest,v['sha256'])
    manifest['inputPackets'][str(focus/'checks.json')]=f.digest(focus/'checks.json');manifest['inputPackets'][str(FCP)]=f.digest(FCP)
    for n,sha in manifest['sourceFiles'].items():
        p=base/'source'/n;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(focus/'source'/n,p);f.require(f.digest(p)==sha,'exact focused source copy')
    parts=old['parts'];f.save(base/'part-inventory.json',parts);f.save(base/'inputs.json',manifest)
    return manifest,tuple(old['cases']),parts


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',type=Path,required=True);ap.add_argument('--stage',choices=('focused','prepare','native','emit','replay','aggregate'),required=True)
    ap.add_argument('--resume-from',type=Path);ap.add_argument('--case');ap.add_argument('--kind',choices=('continuum','current','coordinate'));ap.add_argument('--phase-directory',type=Path);args=ap.parse_args()
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);started=time.monotonic();base=args.run_directory.resolve();base.relative_to(f.STORE)
    prohibit();cases={};combined=None
    if args.stage in ('focused','prepare'):
        base.mkdir(parents=True,exist_ok=False)
        manifest,labels,parts=load(base) if args.stage=='focused' else resume(base,args.resume_from.resolve())
        if args.stage=='focused':cases=focused(base,manifest,labels)
    else:
        manifest=json.loads((base/'inputs.json').read_text());parts=json.loads((base/'part-inventory.json').read_text());labels=tuple(json.loads(CP.read_text())['checks']['cases'])
        if args.stage=='aggregate':
            wanted={'historical_'+h.BASELINE}|{kind+'_'+label for label in labels for kind in ('continuum','current','coordinate') if not(kind=='continuum' and label==h.BASELINE)}
            f.require(set(parts)==wanted and len(parts)==12,'all actual material control output parts accepted')
            for record in parts.values():f.require(f.digest(base/record['directory']/'full.out')==record['sha256'],'completed stream identity')
            combined=original.aggregate(base,parts);f.save(base/'aggregation-checks.json',combined)
        else:
            f.require(args.case in labels,'actual case address');name=args.kind+'_'+args.case;folder=base/'parts'/name;f.require(name not in parts,'no repeated completed output part')
            packet='coordinate-output.pickle' if args.kind=='coordinate' else 'continuum-response.pickle' if args.kind=='continuum' else 'continuum-currents.pickle'
            if args.stage in ('native','emit'):
                folder.mkdir(parents=True,exist_ok=False)
                origin=base/'bundles'/args.case/packet if args.kind=='coordinate' else base/'numerical/cases'/args.case/'continuum'/packet
                m.retain(origin,folder/packet,manifest);f.save(base/'inputs.json',manifest)
            value=f.unpickle(folder/packet);before=f.digest(folder/packet);r=context(base,args.case,value,folder)
            if args.stage=='native':
                f.require(args.kind in ('continuum','current') and not(args.kind=='continuum' and args.case==h.BASELINE),'only missing native output')
                fn=h.native.emitter('MATERIAL_COORDINATE_'+args.case) if args.kind=='continuum' else h.flux.emitter('MATERIAL_COORDINATE_'+args.case)
                join={'wholeOriginalEmitterReplayTail':True,'namespace':'MATERIAL_COORDINATE_'+args.case};f.save(folder/'emitter-join.json',join)
                check=fn(folder,value,r,manifest['sourceFiles'],manifest['inputPackets'],before)
            else:
                emit,replay,join=coordinate_drivers(args.case,folder) if args.kind=='coordinate' else native_drivers(args.case,args.kind,folder)
                if args.stage=='emit':
                    f.save(folder/'emitter-join.json',join);state=emit(folder,value,r);f.atomic_pickle(folder/'emission-state.pickle',state)
                    f.require(f.digest(folder/packet)==before and not engine.PHYSICAL_METADATA.dimensions.constraints,'saved emission packet/unit closure')
                    f.save(folder/'emission-complete.json',{'status':'COMPLETED_EMISSION_AWAITING_REPLAY','packetSha256':before,'transcriptSha256':f.digest(folder/'full.out'),'tags':len(engine.EMISSION_LINES),'keys':len(state['keys'])})
                else:
                    done=json.loads((folder/'emission-complete.json').read_text());transcript=f.digest(folder/'full.out');f.require(done['packetSha256']==before and done['transcriptSha256']==transcript,'completed stream reused without emission')
                    f.save(folder/'validation-join.json',join);check=replay(folder,value,r,manifest['sourceFiles'],manifest['inputPackets'],before)
                    f.require(check['tags']==done['tags'] and len(check['keys'])==done['keys'] and f.digest(folder/'full.out')==transcript,'unchanged full saved stream')
            if args.stage in ('native','replay'):
                parts[name]=recovered.finish_tail()(base,args.case,args.kind,manifest,folder,packet,value,join,before,check)
            f.save(base/'inputs.json',manifest);f.save(base/'part-inventory.json',parts);del value;gc.collect()
    hashes(base,manifest)
    artifacts={str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*') if p.is_file() and 'source' not in p.relative_to(base).parts and p not in (base/'inputs.json',base/'checks.json')}
    status='COMPLETED_CASE_MATERIAL_OUTPUT_INPUTS' if args.stage=='focused' else 'COMPLETED_FOUR_CASE_MATERIAL_OUTPUT' if args.stage=='aggregate' else 'COMPLETED_MATERIAL_OUTPUT_PHASE'
    result={**manifest,'status':status,'stage':args.stage,'parts':parts,'cases':cases,'aggregate':combined,'artifacts':artifacts,'newSolves':0,'newQuadratureNodes':0,'newCurrentConstructions':0,'newGradeExtractions':0,'newIndividualEmissions':int(args.stage in ('native','emit')),'wallSeconds':time.monotonic()-started}
    target=base/'checks.json' if args.stage in ('focused','aggregate') else args.phase_directory/'checks.json'
    f.save(target,result);signal.alarm(0);print(json.dumps(result,indent=2))


if __name__=='__main__':main()
