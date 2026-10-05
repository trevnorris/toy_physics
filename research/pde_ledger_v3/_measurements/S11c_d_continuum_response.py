#!/usr/bin/env python3
"""Solve the finite independent-grade response with computed end maps."""
import argparse
import contextlib
import json
from pathlib import Path
import resource
import signal
import time

import numpy as np
import scipy.linalg as la
import sympy as sp
import S11c_d_continuum_boundary as boundary
import S11c_d_continuum_matrices as interior

f=boundary.f;engine=f.engine;grades=interior.grades;J=boundary.J;G=boundary.G
PLAN=f.M/'S11c_d_continuum_response_plan.md'
PREFIX='CONTINUUM_RESPONSE_LAB_HELD_RHO4_CONSTANT'


def adjoint(series):return {g:v.conj().T for g,v in series.items()}


def gram(series,metric):return J.multiply(J.multiply(adjoint(series),metric),series)


def diagonal_series(parts):return {g:la.block_diag(*(v[g] for v in parts)) for g in G}


def sqrt_series(series):
    values,vectors=np.linalg.eigh(series[(0,0)]);f.require(np.all(values>0),'computed positive open current')
    root={(0,0):(vectors*np.sqrt(values))@vectors.conj().T}
    for g in J.grades:
        known=J.multiply(root,root).get(g,np.zeros_like(root[(0,0)]))
        root[g]=la.solve_sylvester(root[(0,0)],root[(0,0)],series[g]-known)
    residual=boundary.subtract(J.multiply(root,root),series)
    f.require(boundary.norm(residual)<1e-9,'current square-root coefficient equations')
    return root,residual


def load(base):
    coefficients,icp,ip=f.accepted_packet(f.M/'S11c_d_continuum_matrix_checkpoint.json','continuum-matrices.pickle')
    ends,bcp,bp=f.accepted_packet(f.M/'S11c_d_continuum_boundary_checkpoint.json','continuum-boundary.pickle')
    pins={}
    for cp in (icp,bcp):
        for name,sha in cp['sourceFiles'].items():
            f.require(f.digest(f.ROOT/name)==sha,('accepted current source',name))
            if name in pins:f.require(pins[name]==sha,'shared actual source identity')
            pins[name]=sha
    directory=Path(coefficients['finiteCase']);system_path=directory/'finite-system.pickle'
    f.require(f.digest(system_path)==icp['checks']['inputPackets'][str(system_path)],'accepted finite basis/system source')
    system=f.unpickle(system_path)
    reference_path=next(Path(p) for p in ends['inputPackets'] if Path(p).name=='modal.pickle' and Path(p).parent.name=='reference')
    f.require(f.digest(reference_path)==ends['inputPackets'][str(reference_path)],'accepted reference normalization source')
    reference=f.unpickle(reference_path)[0]
    rp=next(Path(p) for p in coefficients['inputPackets'] if Path(p).name=='reduced-action.pickle')
    f.require(f.digest(rp)==coefficients['inputPackets'][str(rp)]==ends['inputPackets'][str(rp)],'common reduced physics source')
    r,d=f.prior.domain.momentum.source.native.source.restore_context(f.unpickle(rp));d.__dict__.update(ends['dimensionState'])
    f.require(coefficients['settings']==system['settings'] and coefficients['size']==len(system['nodes'])==129,'same accepted finite domain/basis/rules')
    f.require(system['nodes'][0]==-64 and system['nodes'][-1]==64,'actual endpoint positions')
    for p in (Path(__file__),PLAN,f.M/'S11c_d_continuum_matrix_checkpoint.json',f.M/'S11c_d_continuum_boundary_checkpoint.json'):
        pins[str(p.resolve().relative_to(f.ROOT))]=f.digest(p)
    operands={str(p):f.digest(p) for p in (ip,bp,system_path,reference_path,rp)}
    for n in pins:
        target=base/'source'/n;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes((f.ROOT/n).read_bytes())
    f.save(base/'inputs.json',{'sourceFiles':pins,'inputPackets':operands,'settings':coefficients['settings'],
        'fieldUnits':[list(map(str,u)) for u in ends['fieldUnits']],'currentUnit':list(map(str,ends['currentUnit'])),
        'scope':'Finite retained continuum coefficient response; approximate modal boundaries and positive regulator.'})
    return r,coefficients,ends,system,reference,pins,operands


def systems(coefficients,ends,system):
    size=len(system['nodes']);derivative=system['derivativeMatrices']
    matrix={g:coefficients['matrices']['total'][(0,*g)].copy() for g in G}
    rhs={g:np.zeros((5*size,4),complex) for g in G}
    for offset,(end,index) in enumerate((('LEFT',0),('RIGHT',size-1))):
        data=ends['ends'][end]
        for g in G:
            for i in range(5):
                row=i*size+index
                for j in range(5):
                    matrix[g][row,j*size:(j+1)*size]=(derivative[1][index] if g==(0,0) and i==j else 0)-data['trace'][g][i,j]*derivative[0][index]
                rhs[g][row,2*offset:2*offset+2]=data['insertion'][g][i]
    return matrix,rhs


def solve(matrix,rhs):
    base=matrix[(0,0)];rows=np.max(abs(base),axis=1);f.require(np.all(rows>0),'nonempty baseline equations')
    scaled=base/rows[:,None];columns=np.linalg.norm(scaled,axis=0);f.require(np.all(columns>0),'nonempty baseline unknowns')
    balanced=scaled/columns[None,:];lu=la.lu_factor(balanced)
    u,s,vh=np.linalg.svd(balanced,full_matrices=False);threshold=np.finfo(float).eps*max(balanced.shape)*s[0]
    rank=int(np.sum(s>threshold));f.require(rank==len(base),'complete baseline rank')
    solutions={};independent={};forcing={}
    for g in G:
        known=J.multiply(matrix,solutions).get(g,np.zeros_like(rhs[g]));forcing[g]=rhs[g]-known
        solutions[g]=la.lu_solve(lu,forcing[g]/rows[:,None])/columns[:,None]
        other_rhs=rhs[g]-J.multiply(matrix,independent).get(g,np.zeros_like(rhs[g]))
        independent[g]=(vh.conj().T@((u.conj().T@(other_rhs/rows[:,None]))/s[:,None]))/columns[:,None]
    residual=boundary.subtract(J.multiply(matrix,solutions),rhs)
    scaled_residual={g:a/rows[:,None] for g,a in residual.items()}
    difference=boundary.subtract(solutions,independent)
    f.require(boundary.norm(scaled_residual)<1e-9 and boundary.norm(difference)/(1+boundary.norm(solutions))<1e-8,'complete coefficient equation and independent solve')
    # Omit the actual mixed operator forcing while retaining both cross terms.
    omitted=forcing[(1,1)]+matrix[(1,1)]@solutions[(0,0)]
    mutation=la.lu_solve(lu,omitted/rows[:,None])/columns[:,None]-solutions[(1,1)]
    f.require(boundary.norm(mutation)>0,'actual mixed operator forcing sensitivity')
    return {'coefficients':solutions,'independentCoefficients':independent,'forcing':forcing,'residual':residual,
        'scaledResidual':scaled_residual,'independentDifference':difference,'mixedForcingMutation':mutation,
        'rowScale':rows,'columnScale':columns,'rank':rank,'singularValues':s,'condition':float(s[0]/s[-1])}


def channels(solution,ends,system,reference):
    size=len(system['nodes']);derivative=system['derivativeMatrices'][0]
    values={g:np.stack([derivative@a[i*size:(i+1)*size] for i in range(5)]) for g,a in solution['coefficients'].items()}
    modal={};open_rows=[];incoming_current=[];outgoing_current=[];in_phases=[];out_phases=[];in_field=[];out_field=[];labels={};boundary_checks={}
    for offset,(end,index) in enumerate((('LEFT',0),('RIGHT',size-1))):
        e=ends['ends'][end];trace={g:a[:,index].copy() for g,a in values.items()}
        for g in G:trace[g][:,2*offset:2*offset+2]-=e['incoming'][g]
        modal[end]=J.multiply(e['outgoingInverse'],trace)
        boundary_checks[end]=boundary.subtract(J.multiply(e['outgoing'],modal[end]),trace)
        out_indices=[];out_maps=[];in_maps=[];record_labels=[];start=0
        for cluster in e['clusters']:
            i=cluster['info']['INDEX'];n=cluster['R'][(0,0)].shape[1];out=cluster['info']['direction']=='outgoing';opened=cluster['info']['kind']=='open'
            source=reference['RECORDS'][i];native=reference['NATIVE_RECORDS'][i]
            mapping=source['FORMS']['FIELD_TO_FLUX_MAP'] if opened else np.eye(n)
            expected=source['FORMS']['RIGHT']@mapping
            f.require(np.allclose(expected,cluster['R'][(0,0)],atol=1e-12,rtol=1e-12),'reference field/flux mode coordinate map')
            if out:
                out_maps.append(mapping)
                if opened:out_indices.extend(range(start,start+n))
                record_labels.append({'record':i,'columns':list(range(start,start+n)),'kind':cluster['info']['kind'],
                    'classifier':str(native['CLASSIFIER_STATUS'])})
            else:in_maps.append(mapping)
            start+=n
        f.require(len(out_indices)==2 and start==7,'complete open and evanescent channel coordinates')
        open_rows.append({g:a[out_indices] for g,a in modal[end].items()});labels[end]=record_labels
        incoming_current.append({g:-e['orientation']*a[5:,5:] for g,a in e['currents']['total'].items()})
        outgoing_current.append({g:e['orientation']*a[np.ix_(out_indices,out_indices)] for g,a in e['currents']['total'].items()})
        in_phases.append(e['incomingOriginPhase']);out_phases.append(e['outgoingOriginPhase'])
        in_field.append(la.block_diag(*in_maps));out_field.append(la.block_diag(*out_maps))
    scattering={g:np.vstack([v[g] for v in open_rows]) for g in G}
    din=diagonal_series(in_phases);dout=diagonal_series(out_phases);dout_inv,dout_res=boundary.inverse(dout)
    origin=J.multiply(J.multiply(dout,scattering),din)
    bin_=gram(din,diagonal_series(incoming_current));bout=gram(dout_inv,diagonal_series(outgoing_current))
    hin,rin=sqrt_series(bin_);hout,rout=sqrt_series(bout);hin_inv,rhin=boundary.inverse(hin)
    normalized=J.multiply(J.multiply(hout,origin),hin_inv)
    field_in=la.block_diag(*in_field);field_in_inv=np.linalg.inv(field_in)
    field_modal={end:{g:out_field[i]@a@field_in_inv for g,a in J.multiply(modal[end],din).items()} for i,end in enumerate(('LEFT','RIGHT'))}
    open_field_maps=[]
    for i,end in enumerate(('LEFT','RIGHT')):
        indices=[j for v in labels[end] if v['kind']=='open' for j in v['columns']]
        open_field_maps.append(out_field[i][np.ix_(indices,indices)])
    field_out=la.block_diag(*open_field_maps)
    # Phase maps act in their computed reference-flux coordinates. Change the
    # coordinates afterward; the field map need not commute with the phase.
    field_scattering={g:field_out@a@field_in_inv for g,a in origin.items()}
    din_inv,din_res=boundary.inverse(din)
    phase_current_in=boundary.subtract(gram(din_inv,bin_),diagonal_series(incoming_current))
    phase_current_out=boundary.subtract(gram(dout,bout),diagonal_series(outgoing_current))
    normalized_current=boundary.subtract(gram(normalized,{(0,0):np.eye(4)}),gram(hin_inv,gram(origin,bout)))
    root_hermitian={name:boundary.subtract(value,adjoint(value)) for name,value in (('incoming',hin),('outgoing',hout))}
    f.require(boundary.norm(boundary_checks)<1e-8,'every outgoing trace reconstruction')
    return {'fields':values,'modalBoundary':modal,'openBoundaryScattering':scattering,'openOriginScattering':origin,
        'outgoingFieldAmplitudes':field_modal,'fieldOriginScattering':field_scattering,
        'fluxOriginScattering':normalized,'incomingCurrentOrigin':bin_,'outgoingCurrentOrigin':bout,
        'incomingCurrentRoot':hin,'outgoingCurrentRoot':hout,'incomingCurrentRootInverse':hin_inv,
        'fieldIncomingMap':field_in,'fieldOutgoingMaps':out_field,'fieldOpenOutgoingMap':field_out,
        'incomingOriginPhase':din,'outgoingOriginPhase':dout,'labels':labels,'residuals':{'boundary':boundary_checks,
        'incomingRoot':rin,'outgoingRoot':rout,'incomingRootInverse':rhin,'outgoingPhaseInverse':dout_res,'incomingPhaseInverse':din_res,
        'phaseCurrentIncoming':phase_current_in,'phaseCurrentOutgoing':phase_current_out,
        'normalizedCurrent':normalized_current,'rootHermitian':root_hermitian}}


def lambda_series(series,ratio):
    result={}
    for (a,b),v in series.items():result[a+b]=result.get(a+b,np.zeros_like(v))+float(ratio)**b*v
    return result


def product(a,b):
    result={}
    for i,x in a.items():
        for j,y in b.items():result[i+j]=result.get(i+j,np.zeros((x.shape[0],y.shape[1]),complex))+x@y
    return result


def open_flux(response,ratio):
    amplitude=lambda_series(response['openOriginScattering'],ratio);current=lambda_series(response['outgoingCurrentOrigin'],ratio)
    incoming=lambda_series(response['incomingCurrentOrigin'],ratio)
    outgoing=product(product({g:a.conj().T for g,a in amplitude.items()},current),amplitude)
    q={}
    for g in range(3):
        numerator=np.diag(outgoing[g]);denominator=np.diag(incoming[0])
        forcing=numerator-sum((np.diag(incoming[j])*q[g-j] for j in range(1,g+1) if j in incoming),np.zeros(4,complex))
        q[g]=forcing/denominator
    return {'amplitudeHomotopy':amplitude,'outgoingCurrentHomotopy':current,'incomingCurrentHomotopy':incoming,
            'openOutgoingFluxHomotopy':outgoing,'openOutgoingFractionCoefficients':q,
            'scope':'Retained-rectangle open-current baseline/interference/quadratic coefficients; no thickness/bulk-loss identification.'}


def emit_result(result,r):
    modes=engine.FullPencilModes.__new__(engine.FullPencilModes);modes.r=r;modes.eta=r.symbols['eta_bg'];modes.sigma=r.symbols['sigma_W']
    zero=(0,0,0);eps=r.symbols['epsilon_shape'];eta=modes.eta;sigma=modes.sigma
    def tensor(name,array,unit=zero,g=(0,0),epsilon=0,literal=False,homotopy=None):
        a=np.asarray(array,complex)
        if a.ndim==1:a=a.reshape(-1,1)
        weight=eps**epsilon*(eta**g[0]*sigma**g[1] if homotopy is None else eta**homotopy)
        body=sp.ImmutableMatrix(*a.shape,[engine.FullPencilModes.number(v)*weight for v in a.ravel()])
        engine.emit(PREFIX+'_'+name,body if literal else modes.compact_fingerprint(body))
        engine.emit('METADATA_'+PREFIX+'_'+name,modes.numeric_metadata(body,lambda p:unit))
    for key in ('openBoundaryScattering','openOriginScattering','fieldOriginScattering','fluxOriginScattering',
                'incomingCurrentOrigin','outgoingCurrentOrigin','incomingCurrentRoot','outgoingCurrentRoot','incomingCurrentRootInverse',
                'incomingOriginPhase','outgoingOriginPhase'):
        for g,a in result['response'][key].items():tensor(key+'_'+''.join(map(str,g)),a,g=g)
    for end,series in result['response']['outgoingFieldAmplitudes'].items():
        for g,a in series.items():tensor(end+'_FIELD_MODAL_'+''.join(map(str,g)),a,g=g,epsilon=1)
    size=result['size']
    for g,a in result['solve']['coefficients'].items():
        for i in range(5):
            unit=tuple(x-y/2 for x,y in zip(result['fieldUnits'][i],result['currentUnit']))
            tensor('FIELD_'+str(i)+'_'+''.join(map(str,g)),a[i*size:(i+1)*size],unit,g,1)
            tensor('SCALED_EQUATION_FRAME_RESIDUAL_'+str(i)+'_'+''.join(map(str,g)),result['solve']['scaledResidual'][g][i*size:(i+1)*size],g=g,literal=True)
            tensor('INDEPENDENT_FIELD_RESIDUAL_'+str(i)+'_'+''.join(map(str,g)),result['solve']['independentDifference'][g][i*size:(i+1)*size],unit,g,1,literal=True)
            equation_unit=tuple(x-y/2 for x,y in zip(result['rowUnits'][i],result['currentUnit']))
            boundary_unit=tuple(v-(1 if j==0 else 0) for j,v in enumerate(unit))
            residual=result['solve']['residual'][g][i*size:(i+1)*size]
            tensor('INTERIOR_EQUATION_RESIDUAL_'+str(i)+'_'+''.join(map(str,g)),residual[1:-1],equation_unit,g,1,literal=True)
            tensor('BOUNDARY_EQUATION_RESIDUAL_'+str(i)+'_'+''.join(map(str,g)),residual[[0,-1]],boundary_unit,g,1,literal=True)
    for i in range(5):
        unit=tuple(x-y/2 for x,y in zip(result['fieldUnits'][i],result['currentUnit']))
        tensor('MIXED_FORCING_MUTATION_FIELD_'+str(i),result['solve']['mixedForcingMutation'][i*size:(i+1)*size],unit,(1,1),1)
    for name in ('incomingRoot','outgoingRoot','incomingRootInverse','outgoingPhaseInverse','incomingPhaseInverse','phaseCurrentIncoming','phaseCurrentOutgoing','normalizedCurrent'):
        for g,a in result['response']['residuals'][name].items():tensor('RESIDUAL_'+name+'_'+''.join(map(str,g)),a,g=g,literal=True)
    for end,series in result['response']['residuals']['boundary'].items():
        for g,a in series.items():
            for i in range(5):
                unit=tuple(x-y/2 for x,y in zip(result['fieldUnits'][i],result['currentUnit']))
                tensor('RESIDUAL_'+end+'_TRACE_'+str(i)+'_'+''.join(map(str,g)),a[i:i+1],unit,g,1,literal=True)
    for name,series in result['response']['residuals']['rootHermitian'].items():
        for g,a in series.items():tensor('RESIDUAL_'+name+'_HERMITIAN_'+''.join(map(str,g)),a,g=g,literal=True)
    for power,a in result['flux']['openOutgoingFractionCoefficients'].items():tensor('OPEN_FRACTION_LAMBDA_'+str(power),a,literal=True,homotopy=power)
    for power,a in result['flux']['openOutgoingFluxHomotopy'].items():tensor('OPEN_FLUX_PER_INCIDENT_CURRENT_UNIT_LAMBDA_'+str(power),a,epsilon=2,homotopy=power)
    for power,a in result['flux']['incomingCurrentHomotopy'].items():tensor('INCOMING_FLUX_PER_CURRENT_UNIT_LAMBDA_'+str(power),a,epsilon=2,homotopy=power)
    boundary.structural_flags(PREFIX+'_MANIFEST',{'sourceFiles':result['sourceFiles'],'inputPackets':result['inputPackets'],
        'labels':result['response']['labels'],'settings':{k:str(v) for k,v in result['settings'].items()},'scope':result['scope'],
        'homotopy':{'lambda':'eta_bg','sigmaOverLambda':str(result['ratio']),
        'status':'After path substitution only; independent grades retained in channel coefficients.',
        'currentConvention':'Current tensors divided by the squared reference incoming flux-amplitude unit; multiply by its physical current unit.',
        'quadraticScope':'Retained rectangle only; omitted pure second-order amplitudes/current can interfere with a nonzero baseline.'}})


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True);args=parser.parse_args()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    def timeout(*_):raise TimeoutError('continuum response budget; preserve completed systems')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(900);started=time.monotonic()
    r,coefficients,ends,system,reference,pins,operands=load(base)
    matrix,rhs=systems(coefficients,ends,system);f.atomic_pickle(base/'coefficient-systems.pickle',{'matrices':matrix,'rhs':rhs})
    solved=solve(matrix,rhs);f.atomic_pickle(base/'coefficient-solutions.pickle',solved)
    response=channels(solved,ends,system,reference);f.atomic_pickle(base/'channel-response.pickle',response)
    f.require(boundary.norm(response['residuals'])<1e-8,'all boundary/phase/current/normalization coefficient residuals')
    specification=json.loads((f.M/'S11c_d_variable_profile_development_input.json').read_text());ratio=sp.Rational(specification['parameters']['W_0'])/sp.Rational(specification['parameters']['L_W'])
    flux=open_flux(response,ratio);remainder={}
    for eta,sigma in ((0.01,0.001),(0.005,0.0005),(0.0025,0.00025)):
        a=boundary.evaluate(matrix,eta,sigma);b=boundary.evaluate(rhs,eta,sigma);direct=la.solve(a,b)
        predicted=boundary.evaluate(solved['coefficients'],eta,sigma);difference=direct-predicted
        remainder[(eta,sigma)]={'direct':direct,'retained':predicted,'difference':difference,
            'maximumReferenceFrame':boundary.norm(difference),'directEquationResidual':a@direct-b}
    f.atomic_pickle(base/'formal-remainders.pickle',remainder)
    result={'solve':solved,'response':response,'flux':flux,'remainders':remainder,'size':len(system['nodes']),
        'fieldUnits':ends['fieldUnits'],'rowUnits':ends['rowUnits'],'currentUnit':ends['currentUnit'],'ratio':ratio,'settings':coefficients['settings'],
        'sourceFiles':pins,'inputPackets':operands,'dimensionState':dict(vars(engine.PHYSICAL_METADATA.dimensions)),
        'scope':'Finite continuum coefficient response with actual end/current variation; approximate modal boundaries, positive regulator and retained rectangle.'}
    f.atomic_pickle(base/'continuum-response.pickle',result);before=f.digest(base/'continuum-response.pickle')
    engine.EMISSION_LINES.clear();engine.PAYLOAD_ENCODER=grades.PayloadEncoder()
    with (base/'full.out').open('x') as stream,contextlib.redirect_stdout(stream):
        emit_result(result,r)
        keys={tag:'s11cdContinuumResponse'+str(i) for i,tag in enumerate(engine.EMISSION_LINES) if not tag.startswith('PY_S11CD_METADATA_')}
        boundary.structural_flags(PREFIX+'_WRITE_KEYS',keys);index=engine.emission_index(engine.EMISSION_LINES);boundary.structural_flags(PREFIX+'_EMISSION_LINES',index)
    entries={}
    for line in grades.decoded_lines(base/'full.out'):
        tag,_,body=line.rstrip('\n').partition(': ');f.require(tag not in entries,'unique response tags');entries[tag]=grades._restore(body)
    original=engine.emit;seen=set()
    def replay(name,value):
        tag='PY_S11CD_'+name;f.require(tag not in seen and entries.get(tag)==engine.cas(value),('response emission replay',tag));seen.add(tag)
    engine.emit=replay
    try:
        emit_result(result,r);boundary.structural_flags(PREFIX+'_WRITE_KEYS',keys);boundary.structural_flags(PREFIX+'_EMISSION_LINES',index)
    finally:engine.emit=original
    f.require(seen==set(entries) and len(keys)==len(set(keys.values())) and not set(keys.values())&set(engine.IMPORT_KEYS),'response payload/key census')
    metadata_paths=0
    for tag,body in entries.items():
        if not tag.startswith('PY_S11CD_METADATA_'):continue
        structural=tag.endswith(('_MANIFEST','_WRITE_KEYS','_EMISSION_LINES'))
        for item in body:
            fields={str(k):v for k,v in (item[1] if structural else item)}
            f.require(len(fields['DIMENSION_L_T_M'])==3 and all(not v.free_symbols for v in fields['DIMENSION_L_T_M']),'resolved response units')
            f.require('MULTIGRADE' in fields and 'EPSILON_LAMBDA_SUPPORT' in fields,'response grades')
            metadata_paths+=1 if structural else len(fields['PATHS'])
    final='PY_S11CD_'+PREFIX+'_EMISSION_LINES';grades.restore_emission_index({str(k):v for k,v in entries[final]},list(entries)[:list(entries).index(final)])
    f.require(before==f.digest(base/'continuum-response.pickle'),'response pre/post packet hash')
    f.require(all(f.digest(f.ROOT/n)==h for n,h in pins.items()) and all(f.digest(Path(n))==h for n,h in operands.items()),'unchanged response sources/inputs')
    summary={'runDirectory':str(base),'sourceFiles':pins,'inputPackets':operands,'unknowns':5*result['size'],'grades':G,'incidentChannels':4,
        'rank':solved['rank'],'condition':solved['condition'],'maximumScaledEquationResidual':boundary.norm(solved['scaledResidual']),
        'maximumIndependentCoefficientDifference':boundary.norm(solved['independentDifference']),
        'mixedForcingMutationMaximum':boundary.norm(solved['mixedForcingMutation']),
        'normalizationResidualMaximum':boundary.norm(response['residuals']),
        'formalRemainderMaxima':{str(g):a['maximumReferenceFrame'] for g,a in remainder.items()},
        'tagCount':len(entries),'writeKeys':len(keys),'metadataPaths':metadata_paths,
        'packetSha256BeforeEmission':before,'packetSha256AfterEmission':f.digest(base/'continuum-response.pickle'),
        'artifacts':{p.name:{'bytes':p.stat().st_size,'sha256':f.digest(p)} for p in base.glob('*') if p.suffix in ('.pickle','.out')},
        'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':result['scope']}
    f.save(base/'checks.json',summary);signal.alarm(0);print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
