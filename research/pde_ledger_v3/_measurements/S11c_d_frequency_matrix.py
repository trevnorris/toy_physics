#!/usr/bin/env python3
"""Actual finite frequency matrices from analytic sources and continued ends."""
import argparse,json,resource,shutil,signal,time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import scipy.linalg as la
import sympy as sp
import S11c_d_frequency_end as end
f=end.f;engine=end.engine;boundary=end.boundary;interior=end.u.q.interior
PLAN=f.M/'S11c_d_frequency_matrix_plan.md'


def load(base):
    data=end.load(base);cp_path=f.M/'S11c_d_frequency_end_focused.json'
    cp=json.loads(cp_path.read_text());f.require(cp['status']=='ACCEPTED_FOCUSED_END_CONTINUATION','accepted complete end continuation')
    directory=Path(cp['runDirectory'])
    for name,h in cp['sourceFiles'].items():f.require(data['pins'][name]==h==f.digest(f.ROOT/name),'unchanged full end continuation helper')
    completed={}
    for name in ('end-seeds.pickle','focused.pickle'):
        p=directory/name;h=cp['artifacts'][name]['sha256'];f.require(f.digest(p)==h,'accepted full end packet');target=base/('accepted-'+name);shutil.copyfile(p,target);f.require(f.digest(target)==h,'byte-identical end reuse');completed[name]=f.unpickle(target);data['operands'][str(p)]=h
    data['endSeeds']=completed['end-seeds.pickle'];data['endFocused']=completed['focused.pickle']
    for p in (Path(__file__),PLAN,cp_path):data['pins'][str(p.resolve().relative_to(f.ROOT))]=f.digest(p)
    for name,h in data['pins'].items():
        target=base/'source'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/name,target)
    data['manifest'].update(sourceFiles=data['pins'],inputPackets=data['operands'],scope='First actual finite analytic-frequency matrix and four-column response at 1-0.01i. Finite rules and approximate modal ends; no pole search or physical complex-frequency current.')
    f.save(base/'inputs.json',data['manifest']);return data


def bind(base,data,frequency):
    r=data['r'];native=data['old']['domain-binding.pickle']['bound'];oldsrc=data['frequencyPacket']['sources'];chart=data['chart'];origin=data['adapter'].input.origin;w=chart['chart']['frequency']
    actual={};mapping={};joins=[]
    for key,record in chart['records'].items():
        old=oldsrc['records'][key];original=old['original'];value=record['analytic'].subs({w:frequency,**origin},simultaneous=True)
        f.require(not value.has(w,*origin),'complete physical frequency/background binding')
        if original in mapping:f.require(mapping[original]==value,'identical source expression has identical analytic binding')
        mapping[original]=value;actual[tuple(record['address'])]=value
        joins.append({'key':key,'address':record['address'],'original':original,'analytic':record['analytic'],'bound':value,'unit':record['unit'],'limits':record['limits']})
    adapter=SimpleNamespace(bind=lambda expression:mapping[expression]);jets={}
    for si in range(35):
        source=native['sources'][0,si];jet=f.source_jets(source,adapter,r);old=data['old']['source-binding.pickle']['jets'][si]
        f.require(jet['probe']==old['probe'] and jet['column']==old['column'],'actual source field and derivative identity');jets[si]=jet
    sources={}
    for key,old in native['sources'].items():
        ti,si=key;record=oldsrc['sourceFrequencies'][si]
        f.require(not record['dependsOnFrequency'] and record['frequency']==old['frequency'],'native real Fourier character remains frequency independent')
        width,k=native['tests'][ti];field=sp.exp(-(r.zp/width)**2+sp.I*k*r.zp);value=sum(a*sp.diff(field,r.zp,n) for n,a in enumerate(jets[si]['coefficients']))
        sources[key]=dict(old,boundAmplitude=value)
    rows=[];profiles=set();source_columns={}
    for old in native['rows']:
        factors=[dict(v,coefficient=actual['factor',old['index'],i]) for i,v in enumerate(old['factors'])]
        rows.append(dict(old,factors=factors));columns={jets[v['sourceIndex']]['column'] for v in factors};f.require(len(columns)==1,'complete native row/source field join');source_columns[old['index']]=next(iter(columns))
        for factor in factors:profiles.update(factor['coefficient'].atoms(sp.Integral))
    f.require(profiles==set(native['profileUnits']),'all six nested profile operands unchanged')
    cells=[]
    for old in native['cells']:
        terms=[]
        for t,(row,prior) in enumerate(old['terms']):
            f.require(old['column']==source_columns[row],'actual native term field join');terms.append((row,actual['cell',old['row'],old['column'],t]))
        cells.append(dict(old,terms=terms))
    orders=sorted({a[1] for a in actual if a[0]=='local'});local={n:sp.ImmutableMatrix(5,5,[actual['local',n,i,j] for i in range(5) for j in range(5)]) for n in orders}
    result={'frequency':frequency,'bound':dict(native,rows=rows,sources=sources,cells=cells),'local':local,'jets':jets,'sourceMap':mapping,'sourceJoins':joins,'settings':data['system']['settings'],'fieldUnits':data['ends']['fieldUnits'],'rowUnits':data['ends']['rowUnits'],'sourceColumns':source_columns,'scope':'Actual analytic source evaluations at approved background contrasts; the symbolic independent grades remain in the accepted chart.'}
    f.atomic_pickle(base/'frequency-binding.pickle',result);return result


def seed_prefix(base,data,binding):
    native=data['old']['domain-binding.pickle']['bound'];r=data['r'];system=data['system'];setting=system['settings'];size=len(system['nodes']);positions=np.array([-12.,0.,15.]);width=float(native['abel']['width'].subs(r.regulator,setting['regulator']))
    fresh=f.BasisMomentum(binding['bound']['rows'],binding['bound']['sources'],r);old=f.BasisMomentum(native['rows'],native['sources'],r)
    fresh.prepare_basis(binding['jets'],positions,setting['sourceBound'],size,setting['sourceNodes']);old.prepare_basis(data['old']['source-binding.pickle']['jets'],positions,setting['sourceBound'],size,setting['sourceNodes'])
    records=[]
    for variables in sorted({tuple(v[0] for v in row['limits']) for row in native['rows']},key=len):
        points,weights=next(old.batches(variables,setting,native['pairs'],width));points,weights=points[:32],weights[:32]
        environment={v:points[:,j] for j,v in enumerate(variables)};environment[r.regulator]=np.full(len(weights),setting['regulator'])
        caches=[({},{}),({},{})]
        for row in native['rows']:
            if tuple(v[0] for v in row['limits'])!=variables:continue
            values=[]
            for worker,(pc,sc) in zip((old,fresh),caches):
                current=worker.rows[row['index']];matrix=np.zeros((len(positions),size),complex)
                for factor in current['factors']:
                    c=worker.coefficient_value(factor['coefficient'],environment,positions,setting,pc);s=worker.source_basis(factor['sourceIndex'],environment,sc);matrix+=(c*weights[:,None]).T@s
                values.append(matrix)
            difference=values[1]-values[0];record={'row':row['index'],'variables':variables,'points':points,'weights':weights,'old':values[0],'analytic':values[1],'difference':difference,'scaledDifference':end.norm(difference)/(1+end.norm(values[0]))}
            f.atomic_pickle(base/f'seed-prefix-{row["index"]}.pickle',record);records.append(record)
    f.require(len(records)==80 and max(v['scaledDifference'] for v in records)<1e-10,'all80 seed analytic/native row contractions')
    return records


def end_maps(base,data,frequency):
    results={};w=complex(frequency)
    for label,seeddata in data['endSeeds'].items():
        clusters=[]
        for seed in seeddata['clusters']:
            if w==1.:state=seed['seedState'];origin='accepted seed'
            else:
                f.require(w==1.-.01j,'only the accepted focused endpoint in this first finite frequency case')
                record=next(v for v in data['endFocused'] if v['label']==label and v['index']==seed['index']);state=record['fine'];origin='accepted fine-step corrected endpoint'
            f.require(state['frequency']==w and end.norm(state['residual'])<2e-12,'actual saved frequency end point')
            clusters.append({'seed':seed,'state':state,'origin':origin})
        result=end.maps(clusters,-64. if label=='LEFT' else 64.);result['clusters']=clusters
        result['openIndices']=[i for i,item in enumerate(seeddata['acceptedChannel']['outgoing']) if item['kind']=='open']
        f.require(result['right'].shape==(5,5) and result['incomingRight'].shape==(5,2),'all five outgoing and two incoming directions')
        f.require(end.norm(result['inverseResidual'])<1e-10 and end.norm(result['traceResidual'])<1e-10,'actual full end maps');results[label]=result
    f.atomic_pickle(base/'frequency-end-maps.pickle',results);return results


def assemble(base,data,binding,row_matrices):
    adapter=SimpleNamespace(bind=lambda expression:binding['sourceMap'][expression]);system=data['system'];r=data['r'];size=len(system['nodes'])
    matrices=interior.assemble(data['packet']['records'],data['packet']['termJoins'],adapter,r,system,row_matrices,False)
    f.require(len(data['packet']['termJoins'])==160 and set(row_matrices)==set(range(80)),'full160 native cell terms and80 row matrices')
    original=matrices['total'][(0,0,0)];local=matrices['local'][(0,0,0)]
    # Independent direct accumulation uses the actual bound cell occurrences.
    direct=np.zeros_like(original)
    for order,coefficients in binding['local'].items():
        for i in range(5):
            for j in range(5):direct[i*size:(i+1)*size,j*size:(j+1)*size]+=interior.values(coefficients[i,j],r,system)[:,None]*system['derivativeMatrices'][order]
    for cell in binding['bound']['cells']:
        if cell['test']!=0:continue
        i,j=cell['row'],cell['column']
        for index,coefficient in cell['terms']:direct[i*size:(i+1)*size,j*size:(j+1)*size]+=interior.values(coefficient,r,system)[:,None]*row_matrices[index]
    difference=interior.differences(original,direct,size);f.require(difference['maximumScaledReferenceFrame']<1e-12,'direct native cell versus term-join matrix assembly')
    result={'local':local,'nonlocal':matrices['nonlocal'][(0,0,0)],'total':original,'direct':direct,'comparison':difference,'rowCount':80,'nativeTerms':160}
    f.atomic_pickle(base/'frequency-interior.pickle',result);return result


def system_and_solve(base,data,binding,interior_data,maps,fixed_scales=None):
    reference=data['system'];size=len(reference['nodes']);derivative=reference['derivativeMatrices'];matrix=interior_data['total'].copy();rhs=np.zeros((5*size,4),complex)
    Yunits=[tuple(data['ends']['rowUnits'][i]) if j not in (0,size-1) else tuple(v-(1 if d==0 else 0) for d,v in enumerate(data['ends']['fieldUnits'][i])) for i in range(5) for j in range(size)]
    for block,(label,index) in enumerate((('LEFT',0),('RIGHT',size-1))):
        current=maps[label]
        for i in range(5):
            row=i*size+index
            for j in range(5):matrix[row,j*size:(j+1)*size]=(derivative[1][index] if i==j else 0)-current['trace'][i,j]*derivative[0][index]
            rhs[row,2*block:2*block+2]=current['forcingAtCommonOrigin'][i]
    if fixed_scales is None:
        rows=np.max(abs(matrix),axis=1);columns=np.linalg.norm(matrix/rows[:,None],axis=0)
    else:rows,columns=fixed_scales
    f.require(np.all(rows>0) and np.all(columns>0),'positive fixed coordinate-frame scales')
    balanced=matrix/rows[:,None]/columns[None,:];u,s,vh=np.linalg.svd(balanced,full_matrices=False);threshold=np.finfo(float).eps*len(matrix)*s[0];rank=int(np.sum(s>threshold))
    system={'frequency':binding['frequency'],'matrix':matrix,'rhs':rhs,'fixedRowScale':rows,'fixedColumnScale':columns,'balanced':balanced,'singularValues':s,'rank':rank,'rankThreshold':threshold,'fieldUnits':data['ends']['fieldUnits'],'sourceRowUnits':Yunits,'incomingAmplitudeUnit':tuple(v/2 for v in data['ends']['currentUnit']),'spaceConvention':'X=C^645 Chebyshev field coefficients, Y=C^645 interior equation/endpoint trace entries in their inherited reference-unit frames; fixed seed diagonal scales and Euclidean norms, U=C^4 fixed seed current coordinates.','settings':reference['settings'],'scope':'Finite analytic retained-operator evaluation, not a continuum grade expansion or a pole search.'}
    f.atomic_pickle(base/'frequency-system.pickle',system);f.require(rank==len(matrix),'regular sampled finite frequency matrix')
    coefficients=la.lu_solve(la.lu_factor(balanced),rhs/rows[:,None])/columns[:,None]
    independent=(vh.conj().T@((u.conj().T@(rhs/rows[:,None]))/s[:,None]))/columns[:,None]
    residual=(matrix@coefficients-rhs)/rows[:,None];difference=coefficients-independent
    f.require(end.norm(residual)<1e-9 and end.norm(difference)/(1+end.norm(coefficients))<1e-8,'complete frequency solve and independent SVD')
    fields=np.stack([derivative[0]@coefficients[i*size:(i+1)*size] for i in range(5)]);modal={};channels=[];boundary_residuals={}
    for block,(label,index) in enumerate((('LEFT',0),('RIGHT',size-1))):
        current=maps[label];value=current['observationAtCommonOrigin']@fields[:,index];value[:,2*block:2*block+2]+=current['directIncomingSubtraction'];modal[label]=value;channels.append(value[current['openIndices']])
        direct_trace=np.stack([derivative[1][index]@coefficients[i*size:(i+1)*size] for i in range(5)])-current['trace']@fields[:,index]
        direct_trace[:,2*block:2*block+2]-=current['forcingAtCommonOrigin'];boundary_residuals[label]=direct_trace
    result={'coefficients':coefficients,'independentDifference':difference,'scaledEquationResidual':residual,'fields':fields,'modalAmplitudes':modal,'openOriginScattering':np.vstack(channels),'boundaryResiduals':boundary_residuals,'fixedFrameCondition':float(s[0]/s[-1]),'rank':rank,'scope':'Analytically continued amplitudes in real-seed current coordinates; no complex-frequency physical current interpretation.'}
    f.atomic_pickle(base/'frequency-solution.pickle',result);return system,result


def seed_case(base,data):
    binding=bind(base,data,sp.S.One);prefix=seed_prefix(base,data,binding)
    old=data['old']['finite-solution.pickle'];rows={i:v for g in old['groups'] for i,v in zip(g['rowIndices'],g['matrices'])}
    interior_data=assemble(base,data,binding,rows);comparison=interior.differences(interior_data['total'],data['system']['unreplacedOperator'],129)
    f.require(comparison['maximumScaledReferenceFrame']<1e-10,'approved finite operator seed replay from joined old quadrature')
    maps=end_maps(base,data,sp.S.One);system,solution=system_and_solve(base,data,binding,interior_data,maps)
    transitions=[];transition_residuals={}
    for label,current in maps.items():
        previous=data['system']['channels'][label];T,_,rank,_=np.linalg.lstsq(previous['incomingValues'],current['incomingRight']@current['incomingPhase'],rcond=None)
        residual=previous['incomingValues']@T-current['incomingRight']@current['incomingPhase'];f.require(rank==2 and end.norm(residual)<1e-10,'actual seed incoming coordinate and phase transition')
        transitions.append(T);transition_residuals[label]=residual
    transition=la.block_diag(*transitions);rhs_residual=system['rhs']-data['system']['rhs']@transition;coeff_residual=solution['coefficients']-old['coefficients']@transition
    f.require(end.norm(rhs_residual)<1e-10 and end.norm(coeff_residual)/(1+end.norm(solution['coefficients']))<1e-8,'seed response compared only after actual incoming coordinate join')
    result={'rowPrefixCount':len(prefix),'maximumPrefixScaled':max(v['scaledDifference'] for v in prefix),'operatorComparison':comparison,'incomingTransition':transition,'transitionResiduals':transition_residuals,'rhsResidual':rhs_residual,'coefficientResidual':coeff_residual,'sourceReuse':'All80 old quadrature matrices reused only at omega=1 after actual source certificates, row-prefix and field/limit joins.'}
    f.atomic_pickle(base/'seed-comparisons.pickle',result);return system,solution,result


def complex_case(base,data,fixed_scales):
    frequency=sp.S.One-sp.I/100;binding=bind(base,data,frequency);r=data['r'];system=data['system'];settings=system['settings'];native=binding['bound'];size=len(system['nodes'])
    worker=f.BasisMomentum(native['rows'],native['sources'],r);worker.prepare_basis(binding['jets'],system['nodes'],settings['sourceBound'],size,settings['sourceNodes'])
    width=float(native['abel']['width'].subs(r.regulator,settings['regulator']));groups=[];rows={}
    for variables in sorted({tuple(v[0] for v in row['limits']) for row in native['rows']},key=len):
        group=worker.matrix_group(variables,settings,native['pairs'],width,system['nodes'],base);groups.append(group);rows.update(zip(group['rowIndices'],group['matrices']))
        f.save(base/'layout-inventory.json',{str(len(v['variables'])):{'path':str(base/f'layout-{len(v["variables"])}.pickle'),'sha256':f.digest(base/f'layout-{len(v["variables"])}.pickle'),'nodes':v['nodes']} for v in groups})
    f.atomic_pickle(base/'frequency-rows.pickle',{'rows':rows,'groups':groups});interior_data=assemble(base,data,binding,rows);maps=end_maps(base,data,frequency);system,solution=system_and_solve(base,data,binding,interior_data,maps,fixed_scales)
    return system,solution,groups


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',required=True,type=Path);args=ap.parse_args();base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False);resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    def timeout(*_):raise TimeoutError('bounded first finite frequency matrix budget; preserve all saved operands')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(900);start=time.monotonic();data=load(base);seed=base/'seed';seed.mkdir();a0,x0,comparisons=seed_case(seed,data)
    fresh=base/'complex';fresh.mkdir();a1,x1,groups=complex_case(fresh,data,(a0['fixedRowScale'],a0['fixedColumnScale']))
    difference={'matrix':a1['matrix']-a0['matrix'],'forcing':a1['rhs']-a0['rhs'],'coefficients':x1['coefficients']-x0['coefficients'],'scattering':x1['openOriginScattering']-x0['openOriginScattering']}
    f.require(end.norm(difference['matrix'])>0 and end.norm(difference['forcing'])>0,'actual matrix and incoming map frequency sensitivity')
    f.atomic_pickle(base/'frequency-comparison.pickle',difference)
    f.require(all(f.digest(f.ROOT/n)==h for n,h in data['pins'].items()) and all(f.digest(Path(n))==h for n,h in data['operands'].items()),'post-construction source/input hashes')
    checks={'status':'COMPLETE_FIRST_FINITE_FREQUENCY_MATRIX','sourceFiles':data['pins'],'inputPackets':data['operands'],'unknowns':645,'incidentDirections':4,'nativeRows':80,'nativeTerms':160,'sourceAmplitudes':35,'seedPrefixRows':80,'maximumSeedPrefixScaled':comparisons['maximumPrefixScaled'],'seedOperatorScaled':comparisons['operatorComparison']['maximumScaledReferenceFrame'],'seedCoefficientDifference':end.norm(comparisons['coefficientResidual']),'seedRank':a0['rank'],'complexRank':a1['rank'],'seedCondition':x0['fixedFrameCondition'],'complexCondition':x1['fixedFrameCondition'],'complexScaledEquationResidual':end.norm(x1['scaledEquationResidual']),'complexIndependentDifference':end.norm(x1['independentDifference']),'complexFrequency':'1-I/100','newMomentumNodes':sum(v['nodes'] for v in groups),'matrixFrequencyChange':end.norm(difference['matrix']),'forcingFrequencyChange':end.norm(difference['forcing']),'amplitudeFrequencyChange':end.norm(difference['scattering']),'artifacts':{str(p.relative_to(base)):{'sha256':f.digest(p),'bytes':p.stat().st_size} for p in base.rglob('*.pickle') if 'source' not in p.relative_to(base).parts},'wallSeconds':time.monotonic()-start,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':data['manifest']['scope']}
    f.save(base/'checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
