#!/usr/bin/env python3
"""Complete finite material-coordinate quadrature and continuum response."""
import argparse,ast,contextlib,copy,inspect,json,resource,shutil,signal,textwrap,time,types
from pathlib import Path
import numpy as np
import sympy as sp
import S11c_d_coordinate_response as c

f=c.f;engine=c.engine;boundary=c.boundary;response=c.response;interior=c.interior;grades=c.grades;G=c.G
PLAN=f.M/'S11c_d_coordinate_response_production_plan.md'
FOCUSED=f.M/'S11c_d_coordinate_response_focused.json'
PREFIX='COORDINATE_RESPONSE_COMPLETE_LAB_HELD_RHO4_CONSTANT'


def native_accumulator():
    original=ast.parse(textwrap.dedent(inspect.getsource(f.BasisMomentum.matrix_group))).body[0];new=copy.deepcopy(original)
    expected=ast.parse("(2 * setting['momentumBound']) ** len(variables)",mode='eval').body;replacement=ast.parse('self.material_box_volume(variables)',mode='eval').body;found=[]
    class Change(ast.NodeTransformer):
        def visit_BinOp(self,n):
            if ast.dump(n)==ast.dump(expected):found.append(1);return copy.deepcopy(replacement)
            return self.generic_visit(n)
    new=Change().visit(new);f.require(len(found)==1,'single native finite-box measure guard')
    restored=copy.deepcopy(new);count=[]
    class Reverse(ast.NodeTransformer):
        def visit_Call(self,n):
            if ast.dump(n)==ast.dump(replacement):count.append(1);return copy.deepcopy(expected)
            return self.generic_visit(n)
    restored=Reverse().visit(restored);f.require(len(count)==1 and ast.dump(restored)==ast.dump(original),'whole native accumulator reverse AST join')
    namespace=dict(vars(f));exec(compile(ast.fix_missing_locations(ast.Module(body=[new],type_ignores=[])),__file__,'exec'),namespace)
    return namespace['matrix_group'],{'wholeAccumulatorReverseAstJoin':True,'changedExpectedVolumeExpressions':1,'originalExpression':ast.unparse(expected),'materialExpression':ast.unparse(replacement),'nativeHelperSha256':f.digest(Path(f.__file__))}


class MaterialProduction(c.MaterialMomentum):
    def material_box_volume(self,variables):
        layouts=[r['limits'] for r in self.rows if tuple(v[0] for v in r['limits'])==variables]
        f.require(layouts and all(limits==layouts[0] for limits in layouts),'same actual material ordered limits')
        widths=[float(upper-lower) for variable,lower,upper in layouts[0]]
        f.require(all(np.isfinite(w) and w>0 for w in widths),'finite positive material integration widths')
        return float(np.prod(widths))


MaterialProduction.matrix_group,ACCUMULATOR_JOIN=native_accumulator()


def load(base):
    accepted=json.loads(FOCUSED.read_text());f.require(accepted['status']=='ACCEPTED_FOCUSED_INSTRUMENT','accepted material instrument')
    old=Path(accepted['runDirectory']);checks=json.loads((old/'checks.json').read_text());inputs=json.loads((old/'inputs.json').read_text())
    f.require(f.digest(old/'checks.json')==accepted['checksSha256'] and checks==accepted['checks'],'accepted focused checks identity')
    outcome=json.loads((old.parent/'coordinate_focused.invocation.json').read_text())
    f.require(outcome['exitCode']==0 and outcome['stderrBytes']==0 and not (old.parent/'coordinate_focused.stderr').stat().st_size and json.loads((old.parent/'coordinate_focused.stdout').read_text())==checks,'completed focused final guards')
    data=c.load(base);f.require(data['pins']==checks['sourceFiles'] and data['operands']==checks['inputPackets'],'exact focused consumed-source and input joins')
    for n,h in checks['sourceFiles'].items():f.require(f.digest(old/'source'/n)==h==f.digest(f.ROOT/n),('unchanged focused source',n))
    for n,h in checks['inputPackets'].items():f.require(f.digest(Path(n))==h,('unchanged accepted input',n))
    for n,v in checks['artifacts'].items():f.require(f.digest(old/n)==v['sha256'],('focused saved artifact',n))
    copies={}
    wanted=['material-binding.pickle','material-local.pickle','material-boundary.pickle',*[n for n in checks['artifacts'] if 'material-current-tables' in n]]
    for n in wanted:
        shutil.copyfile(old/n,base/n);f.require(f.digest(base/n)==checks['artifacts'][n]['sha256'],'byte-identical focused operand reuse');copies[n]=checks['artifacts'][n]
    for p in (Path(__file__),PLAN,FOCUSED):data['pins'][str(p.resolve().relative_to(f.ROOT))]=f.digest(p)
    data['operands'].update({str(old/n):v['sha256'] for n,v in checks['artifacts'].items()});data['operands'][str(old/'checks.json')]=f.digest(old/'checks.json')
    for n in data['pins']:
        target=base/'source'/n;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/n,target)
    manifest={'sourceFiles':data['pins'],'inputPackets':data['operands'],'settings':data['system']['settings'],'input':data['adapter'].input.specification,'focusedReuse':{'directory':str(old),'copies':copies,'checksSha256':f.digest(old/'checks.json')},'accumulatorJoin':ACCUMULATOR_JOIN,'scope':'New full material row quadratures and mixed continuum solve; reused exact focused bindings/local/end operands.'}
    f.save(base/'inputs.json',manifest)
    return data,f.unpickle(base/'material-binding.pickle'),f.unpickle(base/'material-local.pickle'),f.unpickle(base/'material-boundary.pickle'),manifest


def mass_controls(worker,binding):
    results=[]
    for variables in sorted({tuple(l[0] for l in row['limits']) for row in binding['rows']},key=len):
        volume=worker.material_box_volume(variables);physical=(2*binding['settings']['momentumBound'])**len(variables)
        determinant=float(worker.chart.d)**len(variables);residual=volume/determinant-physical
        f.require(abs(residual)<1e-12*(1+physical),'actual material/physical box measure identity')
        wrong=volume-physical;f.require(wrong!=0,'untransformed box measure rejected')
        # A constant shift cancels from a width. Check actual endpoints separately.
        row=next(r for r in binding['rows'] if tuple(l[0] for l in r['limits'])==variables)
        original=worker.data['old']['domain-binding.pickle']['bound']['rows'][row['index']]
        endpoints=[]
        for (_,a,b),(_,ea,eb) in zip(row['limits'],original['limits']):
            endpoints.extend([sp.expand(a-worker.chart.d*ea-worker.chart.kappa),sp.expand(b-worker.chart.d*eb-worker.chart.kappa)])
        f.require(all(v==0 for v in endpoints),'literal affine shifted endpoint joins')
        results.append({'variables':variables,'materialVolume':volume,'physicalVolume':physical,'momentumDeterminant':determinant,'residual':residual,'wrongVolumeResidual':wrong,'endpointResiduals':endpoints})
    return results


def compute(base,data,binding,local,end_route,chart):
    system=data['system'];setting=system['settings'];r=data['r'];size=len(system['nodes']);original_bound=data['old']['domain-binding.pickle']['bound']
    worker=MaterialProduction(binding,data,chart);worker.prepare_basis(binding['jets'],system['nodes'],setting['sourceBound'],size,setting['sourceNodes'])
    controls=mass_controls(worker,binding);f.atomic_pickle(base/'material-measure-controls.pickle',controls)
    width=float(original_bound['abel']['width'].subs(r.regulator,setting['regulator']));groups=[];rows={};row_comparisons={}
    old_rows={i:a for group in data['old']['finite-solution.pickle']['groups'] for i,a in zip(group['rowIndices'],group['matrices'])}
    for variables in sorted({tuple(l[0] for l in row['limits']) for row in binding['rows']},key=len):
        group=worker.matrix_group(variables,setting,original_bound['pairs'],width,system['nodes'],base);groups.append(group)
        for index,array in zip(group['rowIndices'],group['matrices']):
            rows[index]=array;old=old_rows[index];delta=array-old;row_comparisons[index]={'material':array,'eulerian':old,'residual':delta,'scaledReferenceFrame':boundary.norm(delta)/(1+boundary.norm(old))}
        f.save(base/'layout-inventory.json',{str(len(v['variables'])):{'path':f'layout-{len(v["variables"])}.pickle','sha256':f.digest(base/f'layout-{len(v["variables"])}.pickle'),'nodes':v['nodes'],'batches':v['batches']} for v in groups})
    f.require(set(rows)==set(range(80)),'all original material row matrices completed');f.atomic_pickle(base/'material-row-comparisons.pickle',row_comparisons)
    f.require(max(v['scaledReferenceFrame'] for v in row_comparisons.values())<1e-9,'complete material/Eulerian quadrature row agreement')
    shape=system['unreplacedOperator'].shape;nonlocal_={g:np.zeros(shape,complex) for g in binding['support']};term_count=0
    for term in data['packet']['termJoins']:
        i,j,index=term['row'],term['column'],term['integralIndex'];address=('cell',i,j,term['term']);native=original_bound['rows'][index]
        f.require(all(binding['jets'][v['sourceIndex']]['column']==j for v in native['factors']),'complete native source/field/term join')
        for g,a in binding['coefficients'][address].items():
            values=np.broadcast_to(np.asarray(sp.lambdify((r.z,r.regulator),a,'numpy',cse=True)(system['nodes']/float(chart.d),setting['regulator']),complex),system['nodes'].shape)
            f.require(np.isfinite(values).all(),'finite material cell coefficient');nonlocal_[g][i*size:(i+1)*size,j*size:(j+1)*size]+=values[:,None]*rows[index]
        term_count+=1
    f.require(term_count==160,'every native nonlocal cell term');matrices={'local':local['matrices'],'nonlocal':nonlocal_,'total':{g:local['matrices'][g]+a for g,a in nonlocal_.items()}}
    comparisons={kind+'_'+''.join(map(str,g)):interior.differences(value,data['coefficients']['matrices'][kind][g],size) for kind,series in matrices.items() for g,value in series.items()}
    f.atomic_pickle(base/'material-coefficient-matrices.pickle',{'matrices':matrices,'comparisons':comparisons})
    f.require(max(v['maximumScaledReferenceFrame'] for v in comparisons.values())<1e-9,'complete common-frame coefficient operator agreement')
    ends=end_route['ends'];a,b=response.systems({'matrices':matrices},ends,system);old_a,old_b=response.systems(data['coefficients'],data['ends'],system)
    system_comparisons={g:{'matrix':interior.differences(a[g],old_a[g],size),'rhsResidual':b[g]-old_b[g]} for g in G}
    f.atomic_pickle(base/'material-coefficient-systems.pickle',{'matrices':a,'rhs':b,'eulerianMatrices':old_a,'eulerianRhs':old_b,'comparisons':system_comparisons})
    solved=response.solve(a,b);f.atomic_pickle(base/'material-coefficient-solutions.pickle',solved)
    channel=response.channels(solved,ends,system,data['reference']);f.atomic_pickle(base/'material-channel-response.pickle',channel)
    f.require(boundary.norm(channel['residuals'])<1e-8,'actual boundary/phase/current normalization residuals')
    ratio=data['baseline']['ratio'];flux=response.open_flux(channel,ratio)
    channel_comparisons={key:boundary.subtract(channel[key],data['baseline']['response'][key]) for key in ('openBoundaryScattering','openOriginScattering','fieldOriginScattering','fluxOriginScattering','incomingCurrentOrigin','outgoingCurrentOrigin')}
    field_comparisons=boundary.subtract(solved['coefficients'],data['baseline']['solve']['coefficients']);flux_comparisons={n:a-data['baseline']['flux']['openOutgoingFractionCoefficients'][n] for n,a in flux['openOutgoingFractionCoefficients'].items()}
    for key,series in channel_comparisons.items():f.require(boundary.norm(series)/(1+boundary.norm(data['baseline']['response'][key]))<1e-8,('coordinate channel comparison',key))
    f.require(boundary.norm(field_comparisons)/(1+boundary.norm(data['baseline']['solve']['coefficients']))<1e-8,'coordinate full field comparison')
    coordinate={'channels':channel_comparisons,'fields':field_comparisons,'flux':flux_comparisons};f.atomic_pickle(base/'coordinate-comparisons.pickle',coordinate)
    physical_eta=float(data['adapter'].input.origin[r.symbols['eta_bg']]);physical_sigma=float(data['adapter'].input.origin[r.symbols['sigma_W']]);evaluated={k:boundary.evaluate(v,physical_eta,physical_sigma) for k,v in channel_comparisons.items()}
    out={'rows':rows,'rowComparisons':row_comparisons,'groups':groups,'matrices':matrices,'operatorComparisons':comparisons,'systemComparisons':system_comparisons,'solve':solved,'response':channel,'flux':flux,'coordinateComparisons':coordinate,'evaluatedCoordinateDifferences':evaluated,'size':size,'fieldUnits':ends['fieldUnits'],'rowUnits':ends['rowUnits'],'currentUnit':ends['currentUnit'],'blockUnits':data['coefficients']['blockUnits'],'ratio':ratio,'settings':setting,'generators':data['packet']['generators'],'scope':'Complete finite material route in the same transported trial space, common Eulerian operator/boundary/current data constructed before solving. One affine chart; approximate modal boundaries and positive regulator. No general nonaffine or c2 N3/N4/N6 closure.','sourceFiles':data['pins'],'inputPackets':data['operands'],'endRoute':end_route,'measureControls':controls,'geometry':chart.g,'binding':binding,'baseline':data['baseline']}
    return out


def renamed(function,prefix):
    namespace=dict(function.__globals__);namespace['PREFIX']=prefix
    return types.FunctionType(function.__code__,namespace,function.__name__,function.__defaults__,function.__closure__)


def emit_result(result,r):
    # Existing complete response emitter; only its tag namespace is fresh.
    renamed(response.emit_result,PREFIX+'_RESPONSE')(result,r)
    modes=engine.FullPencilModes.__new__(engine.FullPencilModes);modes.r=r;modes.eta=r.symbols['eta_bg'];modes.sigma=r.symbols['sigma_W'];eps=r.symbols['epsilon_shape'];zero=(0,0,0)
    def tensor(name,array,unit=zero,g=(0,0),epsilon=0,literal=False):
        value=np.asarray(array,complex)
        if value.ndim==1:value=value.reshape(-1,1)
        weight=eps**epsilon*modes.eta**g[0]*modes.sigma**g[1];body=sp.ImmutableMatrix(*value.shape,[modes.number(v)*weight for v in value.ravel()])
        engine.emit(PREFIX+'_'+name,body if literal else modes.compact_fingerprint(body));engine.emit('METADATA_'+PREFIX+'_'+name,modes.numeric_metadata(body,lambda p:unit(p[0]//value.shape[1],p[0]%value.shape[1]) if callable(unit) else unit))
    def array_fingerprint(name,array,unit,g=(0,0,0)):
        body=interior.fingerprint(array);body.update(COEFFICIENT_ENTRY_UNIT=unit,COEFFICIENT_GRADE=g,COMPONENT_WEIGHT=sp.prod(x**n for x,n in zip(result['generators'],g)),LAMBDA_ORDER=g[1]+g[2])
        body['NUMERIC_TENSOR_PROJECTIONS']=sp.Tuple(*(v*body['COMPONENT_WEIGHT'] for v in body['NUMERIC_TENSOR_PROJECTIONS']))
        payload=engine.cas(body);engine.emit(PREFIX+'_'+name,payload);engine.emit('METADATA_'+PREFIX+'_'+name,modes.numeric_metadata(payload,lambda p:unit if p and p[0]=='NUMERIC_TENSOR_PROJECTIONS' else zero))
    size=result['size'];fields=result['fieldUnits'];current=result['currentUnit']
    for kind,series in result['matrices'].items():
        for g,value in series.items():
            for i in range(5):
                for j in range(5):array_fingerprint('MATRIX_'+kind+'_'+''.join(map(str,g))+f'_ROW_{i}_COL_{j}',value[i*size:(i+1)*size,j*size:(j+1)*size],result['blockUnits'][i][j],g)
    for index,value in result['rows'].items():
        row=result['binding']['rows'][index];columns={result['binding']['jets'][v['sourceIndex']]['column'] for v in row['factors']};f.require(len(columns)==1,'one complete source field per row');column=next(iter(columns));unit=tuple(x-y for x,y in zip(row['unit'],fields[column]));array_fingerprint('INTEGRAL_ROW_'+str(index),value,unit)
        array_fingerprint('INTEGRAL_ROW_RESIDUAL_'+str(index),result['rowComparisons'][index]['residual'],unit)
    for name,record in result['operatorComparisons'].items():
        for v in record['blocks']:
            unit=result['blockUnits'][v['row']][v['column']];tensor('OPERATOR_RESIDUAL_'+name+f'_ROW_{v["row"]}_COL_{v["column"]}',[[v['absolute']]],unit,g=tuple(map(int,name.rsplit('_',1)[1][1:])),literal=True)
    for name,series in result['coordinateComparisons']['channels'].items():
        for label,values in [('EULERIAN',result['baseline']['response'][name]),('MATERIAL',result['response'][name]),('RESIDUAL',series)]:
            for g,value in values.items():tensor('COORDINATE_'+name+'_'+label+'_'+str(g),value,g=g,literal=label=='RESIDUAL')
    for g,value in result['coordinateComparisons']['fields'].items():
        for i in range(5):tensor('COORDINATE_FIELD_RESIDUAL_'+str(i)+'_'+str(g),value[i*size:(i+1)*size],tuple(x-y/2 for x,y in zip(fields[i],current)),g,1)
    for n,value in result['coordinateComparisons']['flux'].items():
        tensor('COORDINATE_FLUX_RESIDUAL_LAMBDA_'+str(n),value,g=(n,0),literal=True)
    for end,data in result['endRoute']['ends']['ends'].items():
        amplitude=[v['amplitudeUnit'] for v in data['clusters'] for _ in range(v['R'][(0,0)].shape[1])]
        for name,series in data['currents'].items():
            units=lambda i,j:tuple(v-a-b for v,a,b in zip(current,amplitude[i],amplitude[j]))
            for g,value in series.items():
                tensor(end+'_COMMON_CURRENT_'+name+'_'+str(g),value,units,g,2)
                tensor(end+'_MATERIAL_CURRENT_'+name+'_'+str(g),result['endRoute']['materialEnds'][end]['currents'][name][g],units,g,2)
                tensor(end+'_CURRENT_COORDINATE_RESIDUAL_'+name+'_'+str(g),result['endRoute']['proofs'][end]['differences']['currents'][name][g],units,g,2,literal=True)
                tensor(end+'_CURRENT_HERMITIAN_RESIDUAL_'+name+'_'+str(g),value-value.conj().T,units,g,2,literal=True)
        for name,rows,cols,extra in [('trace',fields,fields,(-1,0,0)),('outgoing',fields,amplitude[:5],zero),('outgoingDerivative',fields,amplitude[:5],(-1,0,0)),('incoming',fields,amplitude[5:],zero),('insertion',fields,amplitude[5:],(-1,0,0)),('outgoingInverse',amplitude[:5],fields,zero)]:
            units=lambda i,j,rows=rows,cols=cols,extra=extra:tuple(a-b+q for a,b,q in zip(rows[i],cols[j],extra))
            for g,value in data[name].items():
                tensor(end+'_COMMON_'+name+'_'+str(g),value,units,g)
                tensor(end+'_BOUNDARY_COORDINATE_RESIDUAL_'+name+'_'+str(g),result['endRoute']['proofs'][end]['differences'][name][g],units,g,literal=True)
    for name,packet in result['materialCurrentTables'].items():
        for key,material in packet['material'].items():
            a,b,jl,jr=key;weight=modes.eta**a*modes.sigma**b
            for label,value in [('EULERIAN',packet['original'][key]),('MATERIAL',material)]:
                body=sp.ImmutableMatrix(value)*weight;tag=PREFIX+'_CURRENT_TABLE_'+name.removesuffix('.pickle')+'_'+str(key)+'_'+label
                engine.emit(tag,engine.carrier_fingerprint(body))
                def table_unit(path):
                    i,j=divmod(path[0],5)
                    return tuple(v-u-w+(jl+jr if axis==0 else 0) for axis,(v,u,w) in enumerate(zip(current,fields[i],fields[j])))
                engine.emit('METADATA_'+tag,modes.numeric_metadata(body,table_unit))
    for rec in result['measureControls']:
        degree=len(rec['variables']);unit=(-degree,0,0)
        for name in ('materialVolume','physicalVolume','residual','wrongVolumeResidual'):tensor('BOX_'+str(degree)+'_'+name,[[rec[name]]],unit,literal=True)
    boundary.structural_flags(PREFIX+'_MANIFEST',{'sourceFiles':result['sourceFiles'],'inputPackets':result['inputPackets'],'accumulatorJoin':ACCUMULATOR_JOIN,'scope':result['scope'],'settings':{k:str(v) for k,v in result['settings'].items()},'chartJacobians':{k:str(result['geometry'][k]) for k in ('volumeJacobian','tangentialJacobian','normalJacobian')},'fullRows':len(result['rows']),'actualNewMomentumNodes':sum(v['nodes'] for v in result['groups']),'reusedFocusedOperands':'Byte-identical binding/local/end packets; Eulerian quadratures are comparison operands.'})


def main():
    p=argparse.ArgumentParser();p.add_argument('--run-directory',type=Path,required=True);args=p.parse_args();base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));started=time.monotonic()
    def timeout(*_):raise TimeoutError('material production budget; preserve completed layouts, matrices and responses')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(900)
    data,binding,local,ends,manifest=load(base);chart=c.Chart(data)
    result=compute(base,data,binding,local,ends,chart);result['materialCurrentTables']={n:f.unpickle(base/n) for n in manifest['focusedReuse']['copies'] if 'material-current-tables' in n};result['dimensionState']=dict(vars(engine.PHYSICAL_METADATA.dimensions));f.atomic_pickle(base/'coordinate-response.pickle',result);before=f.digest(base/'coordinate-response.pickle')
    engine.EMISSION_LINES.clear();engine.PAYLOAD_ENCODER=grades.PayloadEncoder()
    with (base/'full.out').open('x') as stream,contextlib.redirect_stdout(stream):
        emit_result(result,data['r']);keys={tag:'s11cdCoordinateResponse'+str(i) for i,tag in enumerate(engine.EMISSION_LINES) if not tag.startswith('PY_S11CD_METADATA_')};boundary.structural_flags(PREFIX+'_WRITE_KEYS',keys);index=engine.emission_index(engine.EMISSION_LINES);boundary.structural_flags(PREFIX+'_EMISSION_LINES',index)
    entries={}
    for line in grades.decoded_lines(base/'full.out'):
        tag,_,body=line.rstrip('\n').partition(': ');f.require(tag not in entries,'unique material response tag');entries[tag]=grades._restore(body)
    old_emit=engine.emit;seen=set()
    def replay(name,value):
        tag='PY_S11CD_'+name;f.require(tag not in seen and entries.get(tag)==engine.cas(value),('material full emission replay',tag));seen.add(tag)
    engine.emit=replay
    try:emit_result(result,data['r']);boundary.structural_flags(PREFIX+'_WRITE_KEYS',keys);boundary.structural_flags(PREFIX+'_EMISSION_LINES',index)
    finally:engine.emit=old_emit
    f.require(seen==set(entries) and len(keys)==len(set(keys.values())) and not set(keys.values())&set(engine.IMPORT_KEYS),'complete fresh material write-key/payload census')
    paths=0
    for tag,body in entries.items():
        if not tag.startswith('PY_S11CD_METADATA_'):continue
        structural=tag.endswith(('_MANIFEST','_WRITE_KEYS','_EMISSION_LINES'))
        for item in body:
            fields={str(k):v for k,v in (item[1] if structural else item)};f.require(len(fields['DIMENSION_L_T_M'])==3 and all(not v.free_symbols for v in fields['DIMENSION_L_T_M']),'material restored dimensions');f.require('MULTIGRADE' in fields and 'EPSILON_LAMBDA_SUPPORT' in fields,'material independent grades and lambda');paths+=1 if structural else len(fields['PATHS'])
    final='PY_S11CD_'+PREFIX+'_EMISSION_LINES';grades.restore_emission_index({str(k):v for k,v in entries[final]},list(entries)[:list(entries).index(final)])
    f.require(f.digest(base/'coordinate-response.pickle')==before and not engine.PHYSICAL_METADATA.dimensions.constraints,'material packet and dimension closure')
    f.require(all(f.digest(f.ROOT/n)==h for n,h in data['pins'].items()) and all(f.digest(Path(n))==h for n,h in data['operands'].items()),'material source/input post hashes')
    checks={'runDirectory':str(base),'sourceFiles':data['pins'],'inputPackets':data['operands'],'rows':len(result['rows']),'nativeTerms':160,'unknowns':5*result['size'],'incidentChannels':4,'rank':result['solve']['rank'],'condition':result['solve']['condition'],'newMomentumNodes':sum(v['nodes'] for v in result['groups']),'rowResidual':max(v['scaledReferenceFrame'] for v in result['rowComparisons'].values()),'operatorResidual':max(v['maximumScaledReferenceFrame'] for v in result['operatorComparisons'].values()),'scaledEquationResidual':boundary.norm(result['solve']['scaledResidual']),'independentCoefficientDifference':boundary.norm(result['solve']['independentDifference']),'coordinateCoefficientDifferences':{k:boundary.norm(v) for k,v in result['coordinateComparisons']['channels'].items()},'evaluatedCoordinateDifferences':{k:boundary.norm(v) for k,v in result['evaluatedCoordinateDifferences'].items()},'tagCount':len(entries),'writeKeys':len(keys),'metadataPaths':paths,'packetSha256BeforeEmission':before,'packetSha256AfterEmission':f.digest(base/'coordinate-response.pickle'),'artifacts':{str(p.relative_to(base)):{'bytes':p.stat().st_size,'sha256':f.digest(p)} for p in base.rglob('*') if p.suffix in ('.pickle','.out') and 'source' not in p.relative_to(base).parts},'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':result['scope']}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
