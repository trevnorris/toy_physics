#!/usr/bin/env python3
"""Literal first-w-derivative sensitivity from accepted reduced operands."""
import argparse,contextlib,copy,json,resource,shutil,signal,time,types
from pathlib import Path
import numpy as np
import sympy as sp
import S11c_d_coordinate_response as c

f=c.f;engine=c.engine;source=c.source;profile=c.profile;response=c.response
boundary=c.boundary;interior=c.interior;grades=c.grades;G=c.G
PLAN=f.M/'S11c_d_first_jet_plan.md';PREFIX='FIRST_JET_CONTROL_LAB_HELD_RHO4_CONSTANT'


def load(base):
    data=c.load(base)
    checkpoint=f.M/'S11c_d_coordinate_response_checkpoint.json';accepted=json.loads(checkpoint.read_text())
    f.require(accepted['status']=='PUBLISHED_ANNEX_VERIFIED','accepted unchanged material scattering route')
    published=f.ROOT/accepted['publication']['path']
    f.require(published.is_symlink() and f.digest(published)==accepted['publication']['sha256'],'actual material publication')
    # The small saved comparison is sufficient to join the unchanged material
    # route to its Eulerian reference; do not load its full quadrature packet.
    path=Path(accepted['runDirectory'])/'coordinate-comparisons.pickle'
    f.require(f.digest(path)==accepted['artifacts'][path.name]['sha256'],'accepted actual two-route comparison')
    for n,h in accepted['sourceFiles'].items():
        f.require(f.digest(f.ROOT/n)==h,('unchanged material control source',n))
    data['materialComparison']=f.unpickle(path);data['materialCheckpoint']=accepted
    data['operands'][str(path)]=f.digest(path)
    for p in (Path(__file__),PLAN,checkpoint):data['pins'][str(p.resolve().relative_to(f.ROOT))]=f.digest(p)
    for n in data['pins']:
        dest=base/'source'/n;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/n,dest)
    manifest={'sourceFiles':data['pins'],'inputPackets':data['operands'],'settings':data['system']['settings'],
              'input':data['adapter'].input.specification,'scope':'One-sided literal first-w-derivative sensitivity; no consistent new profile or channel isolation.'}
    f.save(base/'inputs.json',manifest);return data


def bind(base,data):
    r=data['r'];adapter=data['adapter'];packet=data['packet'];old=data['old'];native=old['domain-binding.pickle']['bound'];cuts=native['cutoffBindings']
    # The unchanged production binder is also checked against accepted rows,
    # jets and the independently assembled full local operator.
    baseline=profile.rebind(base/'baseline-binding.pickle',r,native,packet,adapter)
    baseline_checks=profile.baseline_check(baseline,old['source-binding.pickle'])
    size=len(data['system']['nodes']);local_check=np.zeros_like(data['system']['localMatrix'])
    for n,matrix in baseline['local'].items():
        for i in range(5):
            for j in range(5):local_check[i*size:(i+1)*size,j*size:(j+1)*size]+=interior.values(matrix[i,j],r,data['system'])[:,None]*data['system']['derivativeMatrices'][n]
    local_residual=local_check-data['system']['localMatrix'];f.require(boundary.norm(local_residual)==0,'independent baseline local matrix')
    changed_records={};addresses={};inventory={};target=base/'records';target.mkdir();counts=dict.fromkeys(('local','cell','factor','source'),0)
    mutation_controls=[];endpoint=[]
    for key,item in packet['records'].items():
        original=item['record'];saved=data['coordinate']['records'][key];mutation=source.first_jet_mutation(original['ORIGINAL'],r)
        f.require(mutation==saved['shape'] and saved['original']==original['ORIGINAL'],'accepted literal mutation operand identity')
        selected=mutation['mutated'];record=grades.split(selected,packet['generators'],original['UNIT']);grades.check(record)
        f.require(not record['OUTSIDE_RECTANGLE'],'retained grade rectangle after mutation')
        address=tuple(item['address']);changed=selected!=original['ORIGINAL'];counts[address[0]]+=int(changed)
        # Reverse the same literal derivative substitution, including repeated
        # occurrences. A zero-valued physical endpoint does not waive this check.
        restored=source.first_jet_mutation(selected,r)['mutated'];f.require(restored==original['ORIGINAL'],'literal reversal involution')
        if changed:
            omitted=selected-original['ORIGINAL'];f.require(omitted!=0,'omitted reversal responds symbolically');mutation_controls.append((key,omitted))
        for atom in mutation['occurrences']:
            f.require(isinstance(atom.expr,sp.Derivative) and atom.expr.expr.func==r.profiles['w'] and sum(n for _,n in atom.expr.variable_count)==1,'first w derivative only')
        changed_records[key]=dict(item,record=record);addresses[address]=record
        value={'address':address,'original':original,'mutation':mutation,'record':record,'involutionResidual':restored-original['ORIGINAL'],'unit':original['UNIT']}
        path=target/(key+'.pickle');f.atomic_pickle(path,value);inventory[key]={'path':str(path.relative_to(base)),'sha256':f.digest(path),'bytes':path.stat().st_size};f.save(base/'record-inventory.json',inventory)
    # Actual approved thickness profile and its first derivative at both ends.
    derivative=sp.diff(adapter.input.profiles['w'],r.xi)
    for direction in (-sp.oo,sp.oo):
        before=sp.limit(derivative,r.xi,direction);after=sp.limit(-derivative,r.xi,direction)
        value=sp.limit(adapter.input.profiles['w'],r.xi,direction)
        f.require(before==after==0,'unchanged actual derivative endpoint before mode/current reuse')
        endpoint.append({'direction':direction,'profileLimit':value,'derivativeLimit':before,'reversedDerivativeLimit':after,'residual':after-before})
    rows=[];profiles={};changed_factors=set()
    for oldrow in native['rows']:
        factors=[]
        for fi,oldfactor in enumerate(oldrow['factors']):
            rec=addresses['factor',oldrow['index'],fi];actual=engine.memo_xreplace(adapter.bind(rec['ORIGINAL']),cuts)
            factor=dict(oldfactor,coefficient=actual,symbolicCoefficient=rec['ORIGINAL']);factors.append(factor)
            if actual!=oldfactor['coefficient']:changed_factors.add(oldrow['index'])
            f.require(set(rec['COEFFICIENTS']) in (set(),{(0,0,0)}),'momentum factors grade free')
            for integral in rec['ORIGINAL'].atoms(sp.Integral):
                fresh=engine.memo_xreplace(adapter.bind(integral),cuts);unit=engine.PHYSICAL_METADATA.dimensions.measure(integral)
                if fresh in profiles:f.require(profiles[fresh]==unit,'nested profile unit identity')
                profiles[fresh]=unit
        rows.append(dict(oldrow,factors=factors))
    jets={};sources={};changed_sources=set()
    for si in range(35):
        rec=addresses['source',si];original=native['sources'][0,si];jet=f.source_jets(dict(original,symbolicAmplitude=rec['ORIGINAL']),adapter,r)
        prior=old['source-binding.pickle']['jets'][si]
        f.require(jet['column']==prior['column'] and jet['probe']==prior['probe'],'unchanged source field and position')
        f.require(set(rec['COEFFICIENTS']) in (set(),{(0,0,0)}),'source amplitudes grade free')
        jets[si]=jet
        if jet['coefficients']!=prior['coefficients']:changed_sources.add(si)
        for ti in (0,1):
            original=native['sources'][ti,si];width,k=native['tests'][ti];field=sp.exp(-(r.zp/width)**2+sp.I*k*r.zp)
            actual=sum(a*sp.diff(field,r.zp,n) for n,a in enumerate(jet['coefficients']))
            sources[ti,si]=dict(original,symbolicAmplitude=rec['ORIGINAL'],boundAmplitude=actual)
    changed_rows=[];reused_rows=[];row_joins=[]
    for actual,prior in zip(rows,native['rows']):
        f.require(actual['index']==prior['index'] and actual['limits']==prior['limits'] and actual['sourceLimit']==prior['sourceLimit'],'full row and ordered limit identity')
        sources_used={v['sourceIndex'] for v in actual['factors']}
        f.require(len({jets[i]['column'] for i in sources_used})==1,'complete row input field')
        changed=actual['index'] in changed_factors or bool(sources_used & changed_sources)
        for ti in (0,1):
            for si in sources_used:f.require(sources[ti,si]['frequency']==native['sources'][ti,si]['frequency'],'same actual source frequency')
        if changed:changed_rows.append(actual['index'])
        else:
            f.require(all(a['coefficient']==b['coefficient'] for a,b in zip(actual['factors'],prior['factors'])),'unchanged numerical row factors')
            f.require(all(jets[i]['coefficients']==old['source-binding.pickle']['jets'][i]['coefficients'] for i in sources_used),'unchanged numerical source derivatives')
            reused_rows.append(actual['index'])
        row_joins.append({'index':actual['index'],'sources':sorted(sources_used),'changed':changed,'limits':actual['limits'],'sourceLimit':actual['sourceLimit']})
    f.require(set(changed_rows).isdisjoint(reused_rows) and set(changed_rows+reused_rows)==set(range(80)),'complete computed row partition')
    f.require(set().union(*(factor['coefficient'].atoms(sp.Integral) for row in rows for factor in row['factors']))==set(profiles),'complete nested profile census')
    local_orders=sorted({a[1] for a in addresses if a[0]=='local'})
    local={n:sp.ImmutableMatrix(5,5,[adapter.bind(addresses['local',n,i,j]['ORIGINAL']) for i in range(5) for j in range(5)]) for n in local_orders}
    rebound=dict(native,rows=rows,sources=sources,profileUnits=profiles,profiles=[{'bound':v,'unit':u} for v,u in profiles.items()])
    result={'bound':rebound,'jets':jets,'local':local,'records':changed_records,'originalRecords':packet['records'],'termJoins':packet['termJoins'],'changedRows':changed_rows,'reusedRows':reused_rows,'rowJoins':row_joins,'changedSources':sorted(changed_sources),'changedFactorRows':sorted(changed_factors),'changedRecordCounts':counts,'endpointChecks':endpoint,'baselineChecks':baseline_checks,'baselineLocalResidual':local_residual,'mutationControls':mutation_controls,'settings':data['system']['settings'],'sourceFiles':data['pins'],'inputPackets':data['operands']}
    f.atomic_pickle(base/'first-jet-binding.pickle',result);return result


def compute(base,data,binding):
    r=data['r'];system=data['system'];settings=system['settings'];native=data['old']['domain-binding.pickle']['bound'];size=len(system['nodes'])
    f.require(settings==binding['settings'],'same actual finite quadrature/basis settings')
    old_rows={i:v for group in data['old']['finite-solution.pickle']['groups'] for i,v in zip(group['rowIndices'],group['matrices'])}
    rows={i:old_rows[i] for i in binding['reusedRows']};new_groups=[];work=base/'changed-layouts';work.mkdir()
    selected=[v for v in binding['bound']['rows'] if v['index'] in binding['changedRows']]
    worker=f.BasisMomentum(selected,binding['bound']['sources'],r);worker.prepare_basis(binding['jets'],system['nodes'],settings['sourceBound'],size,settings['sourceNodes'])
    width=float(native['abel']['width'].subs(r.regulator,settings['regulator']))
    for variables in sorted({tuple(v[0] for v in row['limits']) for row in selected},key=len):
        group=worker.matrix_group(variables,settings,native['pairs'],width,system['nodes'],work);new_groups.append(group)
        rows.update(zip(group['rowIndices'],group['matrices']))
        f.save(base/'layout-inventory.json',{str(len(v['variables'])):{'path':str((work/f'layout-{len(v["variables"])}.pickle').relative_to(base)),'nodes':v['nodes'],'sha256':f.digest(work/f'layout-{len(v["variables"])}.pickle')} for v in new_groups})
    f.require(set(rows)==set(range(80)),'all80 computed or source-identical reused rows')
    for i in binding['reusedRows']:f.require(np.array_equal(rows[i],old_rows[i]),'unchanged row arrays')
    f.atomic_pickle(base/'first-jet-rows.pickle',{'rows':rows,'newGroups':new_groups,'reusedRows':binding['reusedRows']})
    matrices=interior.assemble(binding['records'],binding['termJoins'],data['adapter'],r,system,rows,True)
    unsplit=interior.assemble(binding['records'],binding['termJoins'],data['adapter'],r,system,rows,False)
    recombination=interior.differences(interior.recombine(matrices['total'],data['adapter'].input.origin,data['packet']['generators']),unsplit['total'][(0,0,0)],size)
    f.require(recombination['maximumScaledReferenceFrame']<1e-10,'actual grade/unsplit operator recombination')
    f.atomic_pickle(base/'first-jet-matrices.pickle',{'matrices':matrices,'unsplit':unsplit,'recombination':recombination})
    a,b=response.systems({'matrices':matrices},data['ends'],system);f.atomic_pickle(base/'first-jet-systems.pickle',{'matrices':a,'rhs':b})
    solved=response.solve(a,b);f.atomic_pickle(base/'first-jet-solutions.pickle',solved)
    channel=response.channels(solved,data['ends'],system,data['reference']);f.atomic_pickle(base/'first-jet-channels.pickle',channel)
    f.require(boundary.norm(channel['residuals'])<1e-8,'actual boundary/current/phase normalization')
    flux=response.open_flux(channel,data['baseline']['ratio'])
    comparisons={k:boundary.subtract(channel[k],data['baseline']['response'][k]) for k in ('openBoundaryScattering','openOriginScattering','fieldOriginScattering','fluxOriginScattering','incomingCurrentOrigin','outgoingCurrentOrigin')}
    material={k:boundary.subtract(comparisons[k],data['materialComparison']['channels'][k]) for k in comparisons}
    fields=boundary.subtract(solved['coefficients'],data['baseline']['solve']['coefficients'])
    flux_difference={n:v-data['baseline']['flux']['openOutgoingFractionCoefficients'][n] for n,v in flux['openOutgoingFractionCoefficients'].items()}
    eta=float(data['adapter'].input.origin[r.symbols['eta_bg']]);sigma=float(data['adapter'].input.origin[r.symbols['sigma_W']])
    evaluated={k:boundary.evaluate(v,eta,sigma) for k,v in comparisons.items()}
    # Use the actual current matrices and incident denominators for the finite
    # evaluation of each retained amplitude/current polynomial as well.
    def current_fraction(ch):
        sm=boundary.evaluate(ch['openOriginScattering'],eta,sigma);ji=boundary.evaluate(ch['incomingCurrentOrigin'],eta,sigma);jo=boundary.evaluate(ch['outgoingCurrentOrigin'],eta,sigma)
        return np.diag(sm.conj().T@jo@sm)/np.diag(ji)
    fractions={'reversed':current_fraction(channel),'baseline':current_fraction(data['baseline']['response'])}
    fractions['difference']=fractions['reversed']-fractions['baseline']
    differences={'channels':comparisons,'againstUnchangedMaterial':material,'fields':fields,'homotopyFlux':flux_difference,'evaluated':evaluated,'retainedPolynomialCurrentFractions':fractions}
    f.atomic_pickle(base/'first-jet-comparisons.pickle',differences)
    result={'solve':solved,'response':channel,'flux':flux,'size':size,'fieldUnits':data['ends']['fieldUnits'],'rowUnits':data['ends']['rowUnits'],'currentUnit':data['ends']['currentUnit'],'blockUnits':data['coefficients']['blockUnits'],'ratio':data['baseline']['ratio'],'settings':settings,'scope':'Literal first-w-derivative reversal on the Eulerian route. Finite positive regulator and approximate modal boundaries; no isolated tilt channel, new physical profile, nonaffine covariance or c2 source-origin closure.','sourceFiles':data['pins'],'inputPackets':data['operands'],'baseline':data['baseline'],'binding':binding,'matrices':matrices,'rows':rows,'newGroups':new_groups,'recombination':recombination,'comparisons':differences,'generators':data['packet']['generators'],'densityAdvection':data['coordinate']['densityAdvection']}
    return result


def renamed(function,prefix):
    namespace=dict(function.__globals__);namespace['PREFIX']=prefix
    return types.FunctionType(function.__code__,namespace,function.__name__,function.__defaults__,function.__closure__)


def emit_result(result,r):
    renamed(response.emit_result,PREFIX+'_RESPONSE')(result,r)
    modes=engine.FullPencilModes.__new__(engine.FullPencilModes);modes.r=r;modes.eta=r.symbols['eta_bg'];modes.sigma=r.symbols['sigma_W'];zero=(0,0,0)
    def tensor(name,array,unit=zero,g=(0,0),epsilon=0,homotopy=None):
        a=np.asarray(array,complex)
        if a.ndim==1:a=a.reshape(-1,1)
        weight=r.symbols['epsilon_shape']**epsilon*(modes.eta**homotopy if homotopy is not None else modes.eta**g[0]*modes.sigma**g[1])
        body=sp.ImmutableMatrix(*a.shape,[modes.number(v)*weight for v in a.ravel()])
        engine.emit(PREFIX+'_'+name,body);engine.emit('METADATA_'+PREFIX+'_'+name,modes.numeric_metadata(body,lambda p:unit(p[0]//a.shape[1],p[0]%a.shape[1]) if callable(unit) else unit))
    for key,series in result['comparisons']['channels'].items():
        for label,values in [('BASELINE',result['baseline']['response'][key]),('REVERSED',result['response'][key]),('DIFFERENCE',series),('VERSUS_MATERIAL',result['comparisons']['againstUnchangedMaterial'][key])]:
            for g,value in values.items():tensor(key+'_'+label+'_'+str(g),value,g=g)
        tensor(key+'_EVALUATED_DIFFERENCE',result['comparisons']['evaluated'][key])
    for label,values in result['comparisons']['retainedPolynomialCurrentFractions'].items():tensor('RETAINED_POLYNOMIAL_CURRENT_'+label,values)
    size=result['size']
    for g,value in result['comparisons']['fields'].items():tensor('FIELD_DIFFERENCE_'+str(g),value,lambda i,j:tuple(x-y/2 for x,y in zip(result['fieldUnits'][i//size],result['currentUnit'])),g,1)
    for label,values in [('BASELINE',result['baseline']['flux']['openOutgoingFractionCoefficients']),('REVERSED',result['flux']['openOutgoingFractionCoefficients']),('DIFFERENCE',result['comparisons']['homotopyFlux'])]:
        for power,value in values.items():tensor('HOMOTOPY_FRACTION_'+label+'_'+str(power),value,homotopy=power)
    for kind,series in result['matrices'].items():
        if kind not in ('local','nonlocal','total'):continue
        for g,value in series.items():
            for i in range(5):
                for j in range(5):
                    unit=result['blockUnits'][i][j];body=interior.fingerprint(value[i*size:(i+1)*size,j*size:(j+1)*size]);weight=sp.prod(v**n for v,n in zip(result['generators'],g))
                    body.update(COEFFICIENT_ENTRY_UNIT=unit,COEFFICIENT_GRADE=g,COMPONENT_WEIGHT=weight,LAMBDA_ORDER=g[1]+g[2]);body['NUMERIC_TENSOR_PROJECTIONS']=sp.Tuple(*(v*weight for v in body['NUMERIC_TENSOR_PROJECTIONS']))
                    payload=engine.cas(body);name=PREFIX+'_MATRIX_'+kind+'_'+str(g)+f'_{i}_{j}';engine.emit(name,payload);engine.emit('METADATA_'+name,modes.numeric_metadata(payload,lambda p:unit if p and p[0]=='NUMERIC_TENSOR_PROJECTIONS' else zero))
    for v in result['recombination']['blocks']:tensor('RECOMBINATION_'+str(v['row'])+'_'+str(v['column']),[[v['absolute']]],result['blockUnits'][v['row']][v['column']])
    for g in result['newGroups']:
        for index,value in zip(g['rowIndices'],g['actionResidual']):
            row=result['binding']['bound']['rows'][index];column=result['binding']['jets'][row['factors'][0]['sourceIndex']]['column'];unit=tuple(a-b for a,b in zip(row['unit'],result['fieldUnits'][column]));tensor('ROW_ACTION_RESIDUAL_'+str(index),value,unit)
    for key,item in result['binding']['records'].items():
        rec=item['record'];original=result['binding']['originalRecords'][key]['record']['ORIGINAL'];body=sp.Tuple(original,rec['ORIGINAL'],rec['ORIGINAL']-original,*rec['COMPONENTS'].values());name=PREFIX+'_SOURCE_'+key
        engine.emit(name,engine.carrier_fingerprint(body));engine.emit('METADATA_'+name,modes.numeric_metadata(body,lambda p:rec['UNIT']))
        proofs=sp.Tuple(rec['ROUND_TRIP_RESIDUAL'],rec['RECONSTRUCTION_RESIDUAL'],*rec['DERIVATIVE_RESIDUALS'].values(),*rec['NATIVE_COLLECTOR_RESIDUALS'].values())
        engine.emit(name+'_RESIDUALS',proofs);engine.emit('METADATA_'+name+'_RESIDUALS',modes.numeric_metadata(proofs,lambda p:rec['UNIT']))
    boundary.structural_flags(PREFIX+'_MANIFEST',{'sourceFiles':result['sourceFiles'],'inputPackets':result['inputPackets'],'changedRows':result['binding']['changedRows'],'reusedRows':result['binding']['reusedRows'],'changedSources':result['binding']['changedSources'],'changedRecordCounts':result['binding']['changedRecordCounts'],'scope':result['scope'],'settings':{k:str(v) for k,v in result['settings'].items()},'currentEvaluation':'Ratio from retained polynomial amplitude/current at approved grades, not additional parent-theory accuracy.','advection':'Accepted constant-rho4 structural absence; no fictitious A-minus-A control.'})


def main():
    p=argparse.ArgumentParser();p.add_argument('--run-directory',type=Path,required=True);p.add_argument('--focused',action='store_true');p.add_argument('--resume-binding',type=Path);args=p.parse_args()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False);started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    def timeout(*_):raise TimeoutError('first-jet control budget; preserve completed operands')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(900);data=load(base)
    if args.resume_binding:
        old=args.resume_binding.resolve();cp=json.loads((old/'checks.json').read_text());outcome=json.loads((old.parent/'active.json').read_text())
        f.require(outcome['exitCode']==0 and outcome['stderrBytes']==0 and not Path(outcome['stderr']).stat().st_size and json.loads(Path(outcome['stdout']).read_text())==cp,'completed focused binding')
        f.require(cp['sourceFiles']==data['pins'] and cp['inputPackets']==data['operands'],'exact binding source/input joins')
        for n,v in cp['artifacts'].items():
            f.require(f.digest(old/n)==v['sha256'],'focused binding artifact');dest=base/n;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(old/n,dest);f.require(f.digest(dest)==v['sha256'],'byte-identical binding reuse')
        binding=f.unpickle(base/'first-jet-binding.pickle');f.save(base/'binding-reuse.json',{'directory':str(old),'checksSha256':f.digest(old/'checks.json'),'artifacts':cp['artifacts']})
    else:binding=bind(base,data)
    if args.focused:
        summary={'status':'FOCUSED_BINDING_COMPLETE','runDirectory':str(base),'sourceFiles':data['pins'],'inputPackets':data['operands'],'changedRecordCounts':binding['changedRecordCounts'],'changedRows':binding['changedRows'],'reusedRows':binding['reusedRows'],'changedSources':binding['changedSources'],'baselineCoefficientJoins':len(binding['baselineChecks']['coefficientResiduals']),'baselineLocalResidual':boundary.norm(binding['baselineLocalResidual']),'endpointChecks':[{k:str(v) for k,v in x.items()} for x in binding['endpointChecks']],'scriptSha256':f.digest(Path(__file__))}
    else:
        result=compute(base,data,binding);result['dimensionState']=dict(vars(engine.PHYSICAL_METADATA.dimensions));f.atomic_pickle(base/'first-jet-response.pickle',result);before=f.digest(base/'first-jet-response.pickle')
        engine.EMISSION_LINES.clear();engine.PAYLOAD_ENCODER=grades.PayloadEncoder()
        with (base/'full.out').open('x') as stream,contextlib.redirect_stdout(stream):
            emit_result(result,data['r']);keys={tag:'s11cdFirstJetControl'+str(i) for i,tag in enumerate(engine.EMISSION_LINES) if not tag.startswith('PY_S11CD_METADATA_')};boundary.structural_flags(PREFIX+'_WRITE_KEYS',keys);index=engine.emission_index(engine.EMISSION_LINES);boundary.structural_flags(PREFIX+'_EMISSION_LINES',index)
        entries={}
        for line in grades.decoded_lines(base/'full.out'):
            tag,_,body=line.rstrip('\n').partition(': ');f.require(tag not in entries,'unique first-jet tag');entries[tag]=grades._restore(body)
        old_emit=engine.emit;seen=set()
        def replay(name,value):
            tag='PY_S11CD_'+name;f.require(tag not in seen and entries.get(tag)==engine.cas(value),('first-jet full output replay',tag));seen.add(tag)
        engine.emit=replay
        try:emit_result(result,data['r']);boundary.structural_flags(PREFIX+'_WRITE_KEYS',keys);boundary.structural_flags(PREFIX+'_EMISSION_LINES',index)
        finally:engine.emit=old_emit
        f.require(seen==set(entries) and len(keys)==len(set(keys.values())) and not set(keys.values())&set(engine.IMPORT_KEYS),'full first-jet payload/key census')
        paths=0
        for tag,body in entries.items():
            if not tag.startswith('PY_S11CD_METADATA_'):continue
            structural=tag.endswith(('_MANIFEST','_WRITE_KEYS','_EMISSION_LINES'))
            for item in body:
                fields={str(k):v for k,v in (item[1] if structural else item)};f.require(len(fields['DIMENSION_L_T_M'])==3 and all(not v.free_symbols for v in fields['DIMENSION_L_T_M']),'restored units');f.require('MULTIGRADE' in fields and 'EPSILON_LAMBDA_SUPPORT' in fields,'independent grade/lambda');paths+=1 if structural else len(fields['PATHS'])
        final='PY_S11CD_'+PREFIX+'_EMISSION_LINES';grades.restore_emission_index({str(k):v for k,v in entries[final]},list(entries)[:list(entries).index(final)])
        f.require(before==f.digest(base/'first-jet-response.pickle') and not engine.PHYSICAL_METADATA.dimensions.constraints,'packet and dimension closure')
        summary={'status':'COMPLETED_FIRST_JET_RESPONSE','runDirectory':str(base),'sourceFiles':data['pins'],'inputPackets':data['operands'],'changedRows':binding['changedRows'],'reusedRows':binding['reusedRows'],'newMomentumNodes':sum(v['nodes'] for v in result['newGroups']),'rank':result['solve']['rank'],'condition':result['solve']['condition'],'scaledEquationResidual':boundary.norm(result['solve']['scaledResidual']),'independentCoefficientDifference':boundary.norm(result['solve']['independentDifference']),'operatorRecombination':result['recombination']['maximumScaledReferenceFrame'],'evaluatedChannelDifferences':{k:boundary.norm(v) for k,v in result['comparisons']['evaluated'].items()},'retainedPolynomialCurrentDifference':boundary.norm(result['comparisons']['retainedPolynomialCurrentFractions']['difference']),'tagCount':len(entries),'writeKeys':len(keys),'metadataPaths':paths,'packetSha256BeforeEmission':before,'packetSha256AfterEmission':f.digest(base/'first-jet-response.pickle'),'scope':result['scope']}
    f.require(all(f.digest(f.ROOT/n)==h for n,h in data['pins'].items()) and all(f.digest(Path(n))==h for n,h in data['operands'].items()),'source/input post hashes')
    summary.update(artifacts={str(p.relative_to(base)):{'bytes':p.stat().st_size,'sha256':f.digest(p)} for p in base.rglob('*') if p.suffix in ('.pickle','.out') and 'source' not in p.relative_to(base).parts},wallSeconds=time.monotonic()-started,peakRssKiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    f.save(base/'checks.json',summary);signal.alarm(0);print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
