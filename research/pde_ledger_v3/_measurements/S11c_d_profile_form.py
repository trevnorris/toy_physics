#!/usr/bin/env python3
"""Source-driven shape ablation using the accepted finite construction."""
import argparse,contextlib,copy,json,resource,shutil,signal,time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import sympy as sp
import S11c_d_continuum_response as response
import S11c_d_finite_scattering_domain as domain

f=response.f;engine=f.engine;boundary=response.boundary;interior=response.interior;grades=response.grades
PLAN=f.M/'S11c_d_profile_form_plan.md';PREFIX='PROFILE_FORM_LAB_HELD_RHO4_CONSTANT'


def shape_input(specification):
    result=copy.deepcopy(specification)
    result['profiles']={'w':'(1+tanh(xi)**3)/2','m':'(1-tanh(xi)**2)*(1+tanh(xi)/2)/3'}
    return result


def load(base):
    r,coefficients,ends,system,reference,pins,operands=response.load(base)
    result,cp,rp=f.accepted_packet(f.M/'S11c_d_continuum_response_checkpoint.json','continuum-response.pickle')
    grade_packet,gcp,gp=f.accepted_packet(f.M/'S11c_d_continuum_grade_checkpoint.json','continuum-grades.pickle')
    mcp=json.loads((f.M/'S11c_d_continuum_matrix_checkpoint.json').read_text());directory=Path(coefficients['finiteCase'])
    packets={}
    for name in ('domain-binding.pickle','source-binding.pickle','finite-solution.pickle'):
        p=directory/name;f.require(f.digest(p)==mcp['checks']['inputPackets'][str(p)],('accepted finite artifact',name));packets[name]=f.unpickle(p);operands[str(p)]=f.digest(p)
    operands.update({str(rp):f.digest(rp),str(gp):f.digest(gp)})
    specification=json.loads((f.M/'S11c_d_variable_profile_development_input.json').read_text());changed=shape_input(specification)
    baseline=engine.NumericalReducedAction(SimpleNamespace(r=r),{},specification)
    altered=engine.NumericalReducedAction(SimpleNamespace(r=r),{},changed)
    limit_residuals={key:altered.input.limits[key]-baseline.input.limits[key] for key in baseline.input.limits}
    f.require(all(v==0 for v in limit_residuals.values()),'same actual profile endpoints before channel reuse')
    for p in (Path(__file__),PLAN,Path(domain.__file__),f.M/'S11c_d_continuum_response_checkpoint.json',f.M/'S11c_d_continuum_grade_checkpoint.json'):
        pins[str(p.resolve().relative_to(f.ROOT))]=f.digest(p)
    for n,h in pins.items():
        f.require(f.digest(f.ROOT/n)==h,('unchanged consumed source',n));dest=base/'source'/n;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes((f.ROOT/n).read_bytes())
    f.save(base/'inputs.json',{'sourceFiles':pins,'operandHashes':operands,'input':changed,'baselineInput':specification,
        'fieldUnits':{str(i):list(map(str,v)) for i,v in enumerate(ends['fieldUnits'])},
        'equationUnits':[list(map(str,v)) for v in ends['rowUnits']],
        'operatorBlockUnits':[[list(map(str,v)) for v in row] for row in coefficients['blockUnits']],
        'endpointResiduals':{str(k):str(v) for k,v in limit_residuals.items()},'boundaryReuse':'Exact unchanged endpoints, parameters, anchoring and density representative.'})
    return r,coefficients,ends,system,reference,result,grade_packet,packets,baseline,altered,pins,operands


def rebind(base,r,bound,packet,adapter):
    records={tuple(v['address']):v['record'] for v in packet['records'].values()};cuts=bound['cutoffBindings'];rows=[];profiles={};pairs=[]
    for row in bound['rows']:
        factors=[]
        for index,old in enumerate(row['factors']):
            source=records[('factor',row['index'],index)]['ORIGINAL'];f.require(source==old['symbolicCoefficient'],'native factor identity')
            actual=engine.memo_xreplace(adapter.bind(source),cuts)
            factors.append(dict(old,coefficient=actual));pairs.append(('factor',row['index'],index,old['coefficient'],actual))
            for integral in source.atoms(sp.Integral):
                fresh=engine.memo_xreplace(adapter.bind(integral),cuts);unit=engine.PHYSICAL_METADATA.dimensions.measure(integral)
                if fresh in profiles:f.require(profiles[fresh]==unit,'profile unit consistency')
                profiles[fresh]=unit
        rows.append(dict(row,factors=factors))
    jets={i:f.source_jets(bound['sources'][(0,i)],adapter,r) for i in range(35)}
    sources={}
    for key,old in bound['sources'].items():
        ti,si=key;width,k=bound['tests'][ti];field=sp.exp(-(r.zp/width)**2+sp.I*k*r.zp)
        actual=sum(a*sp.diff(field,r.zp,n) for n,a in enumerate(jets[si]['coefficients']))
        sources[key]=dict(old,boundAmplitude=actual)
        f.require(old['symbolicAmplitude']==records[('source',si)]['ORIGINAL'],'native source identity')
    orders=sorted({address[1] for address in records if address[0]=='local'})
    local={n:sp.ImmutableMatrix(5,5,[adapter.bind(records[('local',n,i,j)]['ORIGINAL']) for i in range(5) for j in range(5)]) for n in orders}
    cells=[]
    for cell in bound['cells']:
        ti,i,j=cell['test'],cell['row'],cell['column'];width,k=bound['tests'][ti];field=sp.exp(-(r.z/width)**2+sp.I*k*r.z)
        terms=[]
        for index,(row,old) in enumerate(cell['terms']):
            actual=adapter.bind(records[('cell',i,j,index)]['ORIGINAL']);terms.append((row,actual));pairs.append(('cell',i,j,index,old,actual))
        cells.append(dict(cell,terms=terms,local=sum(a[i,j]*sp.diff(field,r.z,n) for n,a in local.items())))
    actual_profiles=set().union(*(a['coefficient'].atoms(sp.Integral) for row in rows for a in row['factors']))
    f.require(actual_profiles==set(profiles),'complete new nested profile census')
    rebound=dict(bound,rows=rows,sources=sources,cells=cells,profileUnits=profiles,
        profiles=[{'bound':p,'unit':u} for p,u in profiles.items()])
    record={'bound':rebound,'local':local,'jets':jets,'pairs':pairs,'input':adapter.input.specification if hasattr(adapter.input,'specification') else None}
    f.atomic_pickle(base,record)
    return record


def baseline_check(record,accepted):
    residuals=[]
    for pair in record['pairs']:
        before,after=pair[-2:];residual=after-before
        f.require(residual==0,('exact baseline coefficient replay',pair[:-2]));residuals.append(residual)
    for si,jet in record['jets'].items():
        original=accepted['jets'][si]
        f.require(jet['column']==original['column'] and jet['probe']==original['probe'],'source field/probe identity')
        f.require(jet['coefficients']==original['coefficients'],'baseline source-derivative replay')
    return {'coefficientResiduals':residuals,'sourceCount':len(record['jets']),'profileCount':len(record['bound']['profileUnits'])}


def moments(r,baseline,altered):
    t=sp.Symbol('s11cdProfileTanhCoordinate',real=True);xi=r.xi;rows={}
    profiles={name+'_'+p:value for name,adapter in [('baseline',baseline),('altered',altered)] for p,value in adapter.input.profiles.items()}
    profiles['thickness_bump']=(1-sp.tanh(xi)**2)/2
    jacobian=sp.diff(sp.tanh(xi),xi).subs(sp.tanh(xi),t)
    for name,value in profiles.items():
        derivative=sp.diff(value,xi);transformed=sp.cancel(derivative.subs(sp.tanh(xi),t)/jacobian)
        integrand_residual=sp.cancel(transformed*jacobian-derivative.subs(sp.tanh(xi),t))
        lower=sp.limit(sp.tanh(xi),xi,-sp.oo);upper=sp.limit(sp.tanh(xi),xi,sp.oo)
        integral=sp.Integral(transformed,(t,lower,upper));zero=integral.doit()
        limits=(sp.limit(value,xi,-sp.oo),sp.limit(value,xi,sp.oo));jump=limits[1]-limits[0]
        rows[name]={'profile':value,'derivative':derivative,'originalIntegral':sp.Integral(derivative,(xi,-sp.oo,sp.oo)),
            'coordinate':t,'jacobian':jacobian,'transformedIntegral':integral,'integrandResidual':integrand_residual,
            'zeroTransfer':zero,'endLimits':limits,'jump':jump,'residual':sp.simplify(zero-jump)}
    return rows


def compute(base,data,binding):
    r,coefficients,ends,oldsystem,reference,baseline,packet,old,oldadapter,adapter,pins,operands=data
    construct,_,_,joins=domain.adapters();work=base/'ablation';work.mkdir()
    args=(r,binding['bound'],binding['local'],binding['jets'],oldsystem['channels'],pins,operands)
    f.atomic_pickle(work/'source-binding.pickle',{'jets':binding['jets'],'channels':oldsystem['channels']})
    f.save(work/'inputs.json',json.loads((base/'inputs.json').read_text()))
    finite,settings=construct(work,args,129,16,4,256,512,interval=64.,regulator=0.1)
    system=f.unpickle(work/'finite-system.pickle');f.require(system['settings']==oldsystem['settings'],'same numerical settings')
    independent=np.linalg.solve(system['matrix'],system['rhs']);independent_difference=independent-finite['coefficients']
    f.require(boundary.norm(finite['scaledEquationResidual'])<1e-9 and boundary.norm(independent_difference)<1e-8 and finite['rank']==645,'finite equation and independent solve')
    row_matrices={index:array for group in finite['groups'] for index,array in zip(group['rowIndices'],group['matrices'])}
    matrices=interior.assemble(packet['records'],packet['termJoins'],adapter,r,system,row_matrices,True)
    recombined=interior.recombine(matrices['total'],adapter.input.origin,packet['generators'])
    recombination=interior.differences(recombined,system['unreplacedOperator'],129)
    f.require(recombination['maximumScaledReferenceFrame']<1e-10,'actual continuum/full operator recombination')
    f.atomic_pickle(base/'coefficient-matrices.pickle',{'matrices':matrices,'recombination':recombination})
    a,b=response.systems({'matrices':matrices},ends,system);solved=response.solve(a,b)
    f.atomic_pickle(base/'coefficient-response.pickle',solved)
    channels=response.channels(solved,ends,system,reference);f.atomic_pickle(base/'channel-response.pickle',channels)
    finite_phase_in,finite_phase_out=domain.phases(system,finite)
    old_finite=old['finite-solution.pickle'];old_in,old_out=domain.phases(oldsystem,old_finite)
    original=old_out[:,None]*old_finite['boundaryAnchoredFluxBasisScattering']*old_in[None,:]
    altered=finite_phase_out[:,None]*finite['boundaryAnchoredFluxBasisScattering']*finite_phase_in[None,:]
    comparisons={'finiteScattering':altered-original,'finiteCurrent':finite['outgoingFluxRatio']-old_finite['outgoingFluxRatio'],
        'finiteFields':finite['fields']*finite_phase_in[None,None,:]-old_finite['fields']*old_in[None,None,:],
        'continuum':{key:boundary.subtract(channels[key],baseline['response'][key]) for key in ('fieldOriginScattering','fluxOriginScattering')},
        'closedMatching':{end:boundary.subtract(channels['outgoingFieldAmplitudes'][end],baseline['response']['outgoingFieldAmplitudes'][end]) for end in ('LEFT','RIGHT')}}
    result={'finite':finite,'finiteOriginScattering':altered,'baselineFiniteOriginScattering':original,'solve':solved,'response':channels,
        'baselineResponse':baseline['response'],'comparisons':comparisons,'adapterJoins':joins,'settings':settings,
        'independentFiniteDifference':independent_difference,'recombination':recombination}
    f.atomic_pickle(base/'profile-response.pickle',result);return result


def emit_result(result,r):
    modes=engine.FullPencilModes.__new__(engine.FullPencilModes);modes.r=r;modes.eta=r.symbols['eta_bg'];modes.sigma=r.symbols['sigma_W'];eps=r.symbols['epsilon_shape']
    def tensor(name,value,unit=(0,0,0),g=(0,0),epsilon=0,literal=False):
        a=np.asarray(value,complex)
        if a.ndim==1:a=a.reshape(-1,1)
        weight=eps**epsilon*modes.eta**g[0]*modes.sigma**g[1];body=sp.ImmutableMatrix(*a.shape,[modes.number(v)*weight for v in a.ravel()])
        engine.emit(PREFIX+'_'+name,body if literal else modes.compact_fingerprint(body));engine.emit('METADATA_'+PREFIX+'_'+name,modes.numeric_metadata(body,lambda p:unit))
    data=result['computed']
    for key in ('fieldOriginScattering','fluxOriginScattering'):
        for label,series in [('BASE',data['baselineResponse'][key]),('ALTERED',data['response'][key]),('RESIDUAL',data['comparisons']['continuum'][key])]:
            for g,a in series.items():tensor(key+'_'+label+'_'+str(g),a,g=g,literal=label=='RESIDUAL')
    for name,key in [('BASE','baselineFiniteOriginScattering'),('ALTERED','finiteOriginScattering')]:tensor('FINITE_'+name,data[key])
    tensor('FINITE_SCATTERING_DIFFERENCE',data['comparisons']['finiteScattering'],literal=True)
    tensor('FINITE_CURRENT_DIFFERENCE',data['comparisons']['finiteCurrent'],literal=True)
    for end,series in data['comparisons']['closedMatching'].items():
        for g,a in series.items():tensor(end+'_FIELD_MODAL_DIFFERENCE_'+str(g),a,g=g,epsilon=1,literal=True)
    for i in range(5):
        unit=tuple(a-b/2 for a,b in zip(result['fieldUnits'][i],result['currentUnit']))
        tensor('FINITE_INDEPENDENT_RESIDUAL_'+str(i),data['independentFiniteDifference'][i*129:(i+1)*129],unit,epsilon=1)
        for g,a in data['solve']['scaledResidual'].items():tensor('CONTINUUM_SCALED_RESIDUAL_'+str(i)+'_'+str(g),a[i*129:(i+1)*129],g=g,literal=True)
    for name,record in result['moments'].items():
        for key in ('profile','derivative','jacobian','originalIntegral','transformedIntegral','zeroTransfer','endLimits','jump','integrandResidual','residual'):
            body=engine.cas(record[key]);tag=PREFIX+'_'+name+'_'+key;engine.emit(tag,body);engine.emit('METADATA_'+tag,modes.numeric_metadata(body,lambda p:(0,0,0)))
    boundary.structural_flags(PREFIX+'_MANIFEST',{'sourceFiles':result['sourceFiles'],'inputPackets':result['inputPackets'],
        'baselineShapes':result['baselineShapes'],'alteredShapes':result['alteredShapes'],'scope':result['scope'],
        'numericalContrastStatus':'FINITE tags are numerical evaluations at the approved grades, not additional retained coefficients.',
        'bumpScope':'Zero-transfer and end-jump discrimination only; no bump scattering solve.'})


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True);parser.add_argument('--focused',action='store_true')
    parser.add_argument('--resume-binding',type=Path);args=parser.parse_args()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));started=time.monotonic()
    def timeout(*_):raise TimeoutError('profile-form budget; preserve all completed operands')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(900)
    data=load(base);r,coefficients,ends,system,reference,baseline,packet,old,original,altered,pins,operands=data
    old_bound=old['domain-binding.pickle']['bound']
    if args.resume_binding:
        previous=args.resume_binding.resolve();old_checks=json.loads((previous/'checks.json').read_text())
        f.require(old_checks['scriptSha256']==f.digest(Path(__file__)) and old_checks['sourceFiles']==pins and old_checks['inputPackets']==operands,'exact completed binding source/input join')
        for name,item in old_checks['artifacts'].items():
            f.require(f.digest(previous/name)==item['sha256'],'completed binding packet');shutil.copyfile(previous/name,base/name)
            f.require(f.digest(base/name)==item['sha256'],'byte-identical binding reuse')
        new_binding=f.unpickle(base/'altered-binding.pickle');saved=f.unpickle(base/'binding-checks.pickle')
        exact,local_difference,moment=saved['exact'],saved['localOperatorResidual'],saved['moments']
        f.save(base/'binding-reuse.json',{'directory':str(previous),'checksSha256':f.digest(previous/'checks.json'),'artifacts':old_checks['artifacts']})
    else:
        old_binding=rebind(base/'baseline-binding.pickle',r,old_bound,packet,original)
        # Join local coefficients through the accepted independently assembled
        # finite local operator, not a second invocation of this same binder.
        exact=baseline_check(old_binding,old['source-binding.pickle'])
        assembled=np.zeros_like(system['localMatrix'])
        for n,matrix in old_binding['local'].items():
            for i in range(5):
                for j in range(5):assembled[i*129:(i+1)*129,j*129:(j+1)*129]+=interior.values(matrix[i,j],r,system)[:,None]*system['derivativeMatrices'][n]
        local_difference=assembled-system['localMatrix'];f.require(boundary.norm(local_difference)==0,'accepted full local operator baseline identity')
        new_binding=rebind(base/'altered-binding.pickle',r,old_bound,packet,altered);moment=moments(r,original,altered)
        f.atomic_pickle(base/'binding-checks.pickle',{'exact':exact,'localOperatorResidual':local_difference,'moments':moment})
    if args.focused:
        report={'scriptSha256':f.digest(Path(__file__)),'runDirectory':str(base),'exactBaselineCoefficients':len(exact['coefficientResiduals']),
            'sourceAmplitudes':exact['sourceCount'],'nestedProfiles':exact['profileCount'],'fullLocalOperatorResidual':boundary.norm(local_difference),
            'moments':{n:{k:str(v[k]) for k in ('zeroTransfer','jump','residual')} for n,v in moment.items()},
            'artifacts':{p.name:{'bytes':p.stat().st_size,'sha256':f.digest(p)} for p in base.glob('*.pickle')},'sourceFiles':pins,'inputPackets':operands,
            'wallSeconds':time.monotonic()-started,'scope':'Complete source/form binding and zero-transfer checks; no new scattering solve.'}
        f.save(base/'checks.json',report);print(json.dumps(report,indent=2));return
    computed=compute(base,data,new_binding)
    result={'computed':computed,'moments':moment,'baselineShapes':original.input.profiles,'alteredShapes':altered.input.profiles,
        'fieldUnits':ends['fieldUnits'],'currentUnit':ends['currentUnit'],'sourceFiles':pins,'inputPackets':operands,
        'scope':'One finite and retained-continuum shape ablation with identical endpoints; approximate boundaries and positive regulator. Bump moment control only. No global profile claim or pole solve.'}
    f.atomic_pickle(base/'profile-form.pickle',result);before=f.digest(base/'profile-form.pickle')
    engine.EMISSION_LINES.clear();engine.PAYLOAD_ENCODER=grades.PayloadEncoder()
    with (base/'full.out').open('x') as stream,contextlib.redirect_stdout(stream):
        emit_result(result,r);keys={tag:'s11cdProfileForm'+str(i) for i,tag in enumerate(engine.EMISSION_LINES) if not tag.startswith('PY_S11CD_METADATA_')}
        boundary.structural_flags(PREFIX+'_WRITE_KEYS',keys);index=engine.emission_index(engine.EMISSION_LINES);boundary.structural_flags(PREFIX+'_EMISSION_LINES',index)
    entries={}
    for line in grades.decoded_lines(base/'full.out'):
        tag,_,body=line.rstrip('\n').partition(': ');f.require(tag not in entries,'unique ablation tag');entries[tag]=grades._restore(body)
    old_emit=engine.emit;seen=set()
    def replay(name,value):
        tag='PY_S11CD_'+name;f.require(tag not in seen and entries.get(tag)==engine.cas(value),('profile output replay',tag));seen.add(tag)
    engine.emit=replay
    try:
        emit_result(result,r);boundary.structural_flags(PREFIX+'_WRITE_KEYS',keys);boundary.structural_flags(PREFIX+'_EMISSION_LINES',index)
    finally:engine.emit=old_emit
    f.require(seen==set(entries) and len(keys)==len(set(keys.values())) and not set(keys.values())&set(engine.IMPORT_KEYS),'ablation complete key/payload census')
    paths=0
    for tag,body in entries.items():
        if not tag.startswith('PY_S11CD_METADATA_'):continue
        structural=tag.endswith(('_MANIFEST','_WRITE_KEYS','_EMISSION_LINES'))
        for item in body:
            fields={str(k):v for k,v in (item[1] if structural else item)}
            f.require(len(fields['DIMENSION_L_T_M'])==3 and all(not v.free_symbols for v in fields['DIMENSION_L_T_M']),'ablation dimensions')
            f.require('MULTIGRADE' in fields and 'EPSILON_LAMBDA_SUPPORT' in fields,'ablation grades');paths+=1 if structural else len(fields['PATHS'])
    final='PY_S11CD_'+PREFIX+'_EMISSION_LINES';grades.restore_emission_index({str(k):v for k,v in entries[final]},list(entries)[:list(entries).index(final)])
    f.require(before==f.digest(base/'profile-form.pickle') and not engine.PHYSICAL_METADATA.dimensions.constraints,'ablation packet/dimension closure')
    f.require(all(f.digest(f.ROOT/n)==h for n,h in pins.items()) and all(f.digest(Path(n))==h for n,h in operands.items()),'ablation source post hashes')
    summary={'runDirectory':str(base),'sourceFiles':pins,'inputPackets':operands,'finiteRank':computed['finite']['rank'],
        'finiteScaledResidual':boundary.norm(computed['finite']['scaledEquationResidual']),'independentFiniteDifference':boundary.norm(computed['independentFiniteDifference']),
        'continuumScaledResidual':boundary.norm(computed['solve']['scaledResidual']),'operatorRecombination':computed['recombination']['maximumScaledReferenceFrame'],
        'finiteAmplitudeChange':boundary.norm(computed['comparisons']['finiteScattering']),'finiteCurrentChange':boundary.norm(computed['comparisons']['finiteCurrent']),
        'continuumCoefficientChanges':{k:{str(g):boundary.norm(v) for g,v in series.items()} for k,series in computed['comparisons']['continuum'].items()},
        'tagCount':len(entries),'writeKeys':len(keys),'metadataPaths':paths,'packetSha256BeforeEmission':before,'packetSha256AfterEmission':f.digest(base/'profile-form.pickle'),
        'artifacts':{str(p.relative_to(base)):{'bytes':p.stat().st_size,'sha256':f.digest(p)} for p in base.rglob('*') if p.suffix in ('.pickle','.out') and 'source' not in p.relative_to(base).parts},
        'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':result['scope']}
    f.save(base/'checks.json',summary);signal.alarm(0);print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
