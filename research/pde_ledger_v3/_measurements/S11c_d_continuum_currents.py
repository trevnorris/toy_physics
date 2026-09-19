#!/usr/bin/env python3
"""Contract saved continuum amplitudes with source-derived channel currents."""
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
import S11c_d_continuum_response as response

f=response.f;boundary=response.boundary;engine=response.engine;grades=response.grades;G=response.G;J=response.J
PLAN=f.M/'S11c_d_continuum_currents_plan.md'
PREFIX='CONTINUUM_CURRENTS_LAB_HELD_RHO4_CONSTANT'


def multiply(a,b):
    result={}
    for g,x in a.items():
        for h,y in b.items():
            k=tuple(i+j for i,j in zip(g,h));value=x@y
            result[k]=result.get(k,np.zeros_like(value))+value
    return result


def quadratic(a,b):return multiply(multiply(response.adjoint(a),b),a)


def subtract(a,b):
    template=next(iter(a.values())) if a else next(iter(b.values()))
    return {g:a.get(g,np.zeros_like(template))-b.get(g,np.zeros_like(template)) for g in set(a)|set(b)}


def quotient(numerator,denominator,degree=2):
    result={};residual={}
    for i in range(degree+1):
        known=sum((np.diag(denominator[j])*result[i-j] for j in denominator if 0<j<=i),np.zeros(denominator[0].shape[0],complex))
        coefficient=np.diag(numerator.get(i,np.zeros_like(denominator[0])))
        result[i]=(coefficient-known)/np.diag(denominator[0])
        residual[i]=sum((np.diag(denominator[j])*result[i-j] for j in denominator if j<=i),np.zeros_like(result[i]))-coefficient
    return result,residual


def direct_quadratic(amplitude,metric,computed):
    # Independent scalar contraction; all interference entries are retained.
    expected={g:np.zeros_like(v) for g,v in computed.items()}
    for g,x in amplitude.items():
        for h,b in metric.items():
            for k,y in amplitude.items():
                degree=tuple(i+j+l for i,j,l in zip(g,h,k))
                for a in range(x.shape[1]):
                    for z in range(y.shape[1]):
                        expected[degree][a,z]+=sum(x[i,a].conjugate()*b[i,j]*y[j,z]
                            for i in range(len(x)) for j in range(len(y)))
    return subtract(expected,computed)


def load(base):
    result,cp,rp=f.accepted_packet(f.M/'S11c_d_continuum_response_checkpoint.json','continuum-response.pickle')
    ends,bcp,bp=f.accepted_packet(f.M/'S11c_d_continuum_boundary_checkpoint.json','continuum-boundary.pickle')
    f.require(result['inputPackets'][str(bp)]==f.digest(bp),'response/end identity')
    reference_path=next(Path(n) for n in result['inputPackets'] if Path(n).name=='modal.pickle')
    f.require(f.digest(reference_path)==result['inputPackets'][str(reference_path)],'reference classifier identity')
    reference=f.unpickle(reference_path)[0]
    pencil={};operands={str(p):f.digest(p) for p in (rp,bp,reference_path)}
    for end in ('LEFT','RIGHT'):
        pencil[end],_,p=f.accepted_packet(f.M/'S11c_d_continuum_boundary_checkpoint.json',end.lower()+'-pencil.pickle')
        operands[str(p)]=f.digest(p)
    reduction=next(Path(n) for n in result['inputPackets'] if Path(n).name=='reduced-action.pickle')
    f.require(f.digest(reduction)==result['inputPackets'][str(reduction)],'accepted reduction context')
    r,d=f.prior.domain.momentum.source.native.source.restore_context(f.unpickle(reduction));d.__dict__.update(result['dimensionState'])
    operands[str(reduction)]=f.digest(reduction);pins=dict(cp['sourceFiles'])
    for n,h in pins.items():f.require(f.digest(f.ROOT/n)==h,('unchanged consumed source',n))
    for p in (Path(__file__),PLAN,f.M/'S11c_d_continuum_response_checkpoint.json',f.M/'S11c_d_continuum_boundary_checkpoint.json'):
        pins[str(p.resolve().relative_to(f.ROOT))]=f.digest(p)
    for n in pins:
        dest=base/'source'/n;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes((f.ROOT/n).read_bytes())
    f.save(base/'inputs.json',{'sourceFiles':pins,'inputPackets':operands,'settings':result['settings'],
        'scope':'Retained finite continuum channel and normal-current bookkeeping; bulk-depth kinematic domain is separate.'})
    return r,result,ends,reference,pencil,pins,operands


def selectors(result,ends,reference):
    labels=[];census={};outgoing=[]
    for end in ('LEFT','RIGHT'):
        v=ends['ends'][end];records=[];indices=[];offset=0
        f.require(len(v['census'])==18,'full original candidate dispositions')
        for cluster in v['clusters']:
            i=cluster['info']['INDEX'];n=cluster['R'][(0,0)].shape[1]
            status=str(reference['NATIVE_RECORDS'][i]['CLASSIFIER_STATUS'])
            f.require(status in ('THICKNESS_LIKE','TRANSVERSE_LIKE'),'resolved actual channel classifier')
            f.require(reference['RECORDS'][i]['NULLITY']==n,'full classifier subspace')
            records.append({'record':i,'columns':list(range(offset,offset+n)),'direction':cluster['info']['direction'],
                'kind':cluster['info']['kind'],'classifier':status})
            if cluster['info']['direction']=='outgoing' and cluster['info']['kind']=='open':
                labels.extend([{'end':end,'record':i,'direction':j,'classifier':status} for j in range(n)])
                indices.extend(range(offset,offset+n))
            offset+=n
        old=result['response']['labels'][end]
        f.require([{k:v for k,v in row.items() if k!='direction'} for row in records if row['direction']=='outgoing']==old,'response/classifier coordinate join')
        f.require(offset==7 and len(indices)==2,'complete outgoing and incident basis census')
        census[end]=records;outgoing.append(indices)
    identity=np.eye(len(labels));projectors={name:np.diag([float(v['classifier']==status) for v in labels])
        for name,status in [('transverse','TRANSVERSE_LIKE'),('thickness','THICKNESS_LIKE')]}
    projectors['allOpen']=identity
    residual=projectors['transverse']+projectors['thickness']-identity
    f.require(boundary.norm(residual)==0,'classified complete open subspace')
    return {'labels':labels,'census':census,'indices':outgoing,'projectors':projectors,'partitionResidual':residual,
        'ranks':{n:int(np.linalg.matrix_rank(v)) for n,v in projectors.items()}}


def open_metrics(result,ends,selection):
    data=result['response'];din=data['incomingOriginPhase'];dout_inverse,residual=boundary.inverse(data['outgoingOriginPhase'])
    fi=np.linalg.inv(data['fieldIncomingMap']);fo=np.linalg.inv(data['fieldOpenOutgoingMap']);parts={}
    for part in ('slab','bulk','total'):
        outgoing=[];incoming=[]
        for e,indices in zip(('LEFT','RIGHT'),selection['indices']):
            v=ends['ends'][e];s=v['orientation']
            outgoing.append({g:s*a[np.ix_(indices,indices)] for g,a in v['currents'][part].items()})
            incoming.append({g:-s*a[5:,5:] for g,a in v['currents'][part].items()})
        raw_out=response.gram(dout_inverse,response.diagonal_series(outgoing));raw_in=response.gram(din,response.diagonal_series(incoming))
        parts[part]={'outgoing':{g:fo.conj().T@a@fo for g,a in raw_out.items()},
                     'incoming':{g:fi.conj().T@a@fi for g,a in raw_in.items()}}
    joins={key:subtract({g:parts['slab'][key][g]+parts['bulk'][key][g] for g in G},parts['total'][key]) for key in ('outgoing','incoming')}
    f.require(boundary.norm(joins)<1e-10,'actual slab plus bulk open-current source join')
    return parts,{'parts':joins,'phaseInverse':residual}


def amplitude_bookkeeping(series):
    baseline={g:(a.copy() if g==(0,0) else np.zeros_like(a)) for g,a in series.items()}
    zero={g:(a.copy() if g==(1,0) else np.zeros_like(a)) for g,a in series.items()}
    first={g:(a.copy() if g[1] else np.zeros_like(a)) for g,a in series.items()}
    delta=subtract(series,baseline);reconstruction={g:baseline[g]+zero[g]+first[g]-series[g] for g in G}
    return {'baseline':baseline,'zeroJetContrast':zero,'firstJet':first,'delta':delta,'full':series,'reconstructionResidual':reconstruction}


def construct_open(result,metrics,selection):
    original=result['response']['fieldOriginScattering'];ratio=result['ratio'];records={};incoming=response.lambda_series(metrics['total']['incoming'],ratio)
    for name,projector in selection['projectors'].items():
        amplitude={g:projector@a for g,a in original.items()};components=amplitude_bookkeeping(amplitude)
        parts={key:quadratic(amplitude,v['outgoing']) for key,v in metrics.items()}
        direct={key:direct_quadratic(amplitude,v['outgoing'],parts[key]) for key,v in metrics.items()}
        induced=quadratic(components['delta'],metrics['total']['outgoing']);baseline=quadratic(components['baseline'],metrics['total']['outgoing'])
        homotopy={key:response.lambda_series(value,ratio) for key,value in parts.items()}
        induced_h=response.lambda_series(induced,ratio)
        total_fraction,qres=quotient(homotopy['total'],incoming);induced_fraction,ires=quotient(induced_h,incoming)
        amplitude_h=response.lambda_series(amplitude,ratio)
        records[name]={'components':components,'currentParts':parts,'homotopyParts':homotopy,'incidentHomotopy':incoming,
            'totalFraction':total_fraction,'inducedCurrent':induced,'inducedHomotopy':induced_h,'inducedFraction':induced_fraction,
            'totalMinusBaselineVersusInduced':subtract(subtract(parts['total'],baseline),induced),
            'weakAmplitude':amplitude_h.get(1,np.zeros_like(amplitude[(0,0)])),
            'weakTotalFraction':total_fraction,'weakInducedQuadratic':induced_fraction[2],
            'residuals':{'directContractions':direct,'totalQuotient':qres,'inducedQuotient':ires,
                'currentPartJoin':subtract({g:parts['slab'][g]+parts['bulk'][g] for g in parts['total']},parts['total']),
                'amplitudeReconstruction':components['reconstructionResidual']}}
        f.require(boundary.norm(records[name]['residuals'])<1e-8,'complete channel amplitude and current contractions')
    changed={**metrics['total']['outgoing'],(0,0):2*metrics['total']['outgoing'][(0,0)]}
    mutation=subtract(quadratic(original,changed),records['allOpen']['currentParts']['total'])
    f.require(boundary.norm(mutation)>0,'actual outgoing current coefficient sensitivity')
    return records,mutation


def end_currents(result,ends):
    data=result['response'];fi=np.linalg.inv(data['fieldIncomingMap']);din=data['incomingOriginPhase'];records={}
    for offset,end in enumerate(('LEFT','RIGHT')):
        v=ends['ends'][end];incoming=np.eye(4)[2*offset:2*offset+2]
        amplitudes={g:np.vstack((a,incoming if g==(0,0) else np.zeros_like(incoming))) for g,a in data['modalBoundary'][end].items()}
        amplitudes={g:a@fi for g,a in J.multiply(amplitudes,din).items()}
        po=np.diag([float(row<5) for row in range(7)]);pi=np.eye(7)-po
        out={g:po@a for g,a in amplitudes.items()};inc={g:pi@a for g,a in amplitudes.items()}
        parts={};residual={}
        for part,current in v['currents'].items():
            signed={g:v['orientation']*a for g,a in current.items()}
            full=quadratic(amplitudes,signed);outgoing=quadratic(out,signed);incoming=quadratic(inc,signed)
            cross1=multiply(multiply(response.adjoint(out),signed),inc);cross2=multiply(multiply(response.adjoint(inc),signed),out)
            cross={g:cross1[g]+cross2[g] for g in cross1}
            parts[part]={'full':full,'outgoing':outgoing,'incoming':incoming,'interference':cross}
            residual[part]={'decomposition':subtract(full,{g:outgoing[g]+incoming[g]+cross[g] for g in full}),
                            'direct':direct_quadratic(amplitudes,signed,full)}
        f.require(boundary.norm(residual)<1e-8,'complete finite-boundary currents including closed and cross-mode terms')
        records[end]={'amplitude':amplitudes,'parts':parts,'residuals':residual}
    return records


def bulk_domains(pencils):
    records={}
    for end,p in pencils.items():
        rows=[]
        for wave in p['waves']:
            q=next(v for v in wave.free_symbols if 'Acoustic' in v.name)
            k=next(v for v in wave.free_symbols if 'Current' in v.name)
            polynomial=sp.Poly(wave,q);f.require(polynomial.degree()==2 and polynomial.nth(1)==0,'actual quadratic acoustic depth relation')
            squared=sp.cancel(-polynomial.nth(0)/polynomial.nth(2));critical=sp.solve(sp.diff(squared,k),k)
            f.require(len(critical)==1,'computed stationary momentum')
            maximum=sp.simplify(squared.subs(k,critical[0]));curvature=sp.diff(squared,k,2)
            real=sp.Dummy('realNormalMomentum',real=True)
            f.require(curvature.is_negative and (squared-maximum).subs(k,real).is_nonpositive,'computed global maximum on real momenta')
            reconstruction=sp.expand(polynomial.nth(2)*(q**2-squared)-wave)
            f.require(reconstruction==0,'acoustic polynomial reconstruction')
            rows.append({'wave':wave,'depthMomentum':q,'normalMomentum':k,'depthSquared':squared,'criticalMomentum':critical[0],
                'globalMaximum':maximum,'curvature':curvature,'reconstructionResidual':reconstruction,
                'realPropagationSet':sp.solve_univariate_inequality(squared.subs(k,real)>=0,real,relational=False),
                'domain':'real normal momentum; approved real frequency and tangential momentum'})
        records[end]=rows
    return records


def emit_result(result,r):
    mode=engine.FullPencilModes.__new__(engine.FullPencilModes);mode.r=r;mode.eta=r.symbols['eta_bg'];mode.sigma=r.symbols['sigma_W']
    eps=r.symbols['epsilon_shape'];zero=(0,0,0);current=result['currentUnit']
    def tensor(name,array,unit=zero,degree=(0,0),epsilon=0,literal=False,homotopy=None):
        a=np.asarray(array,complex)
        if a.ndim==1:a=a.reshape(-1,1)
        weight=eps**epsilon*(mode.eta**degree[0]*mode.sigma**degree[1] if homotopy is None else mode.eta**homotopy)
        body=sp.ImmutableMatrix(*a.shape,[mode.number(v)*weight for v in a.ravel()]);tag=PREFIX+'_'+name
        engine.emit(tag,body if literal else mode.compact_fingerprint(body));engine.emit('METADATA_'+tag,mode.numeric_metadata(body,lambda p:unit))
    for name,projector in result['selection']['projectors'].items():tensor('SELECTOR_'+name,projector,literal=True)
    tensor('SELECTOR_PARTITION_RESIDUAL',result['selection']['partitionResidual'],literal=True)
    for part,values in result['metrics'].items():
        for key,series in values.items():
            for g,a in series.items():tensor('FIELD_CURRENT_METRIC_'+part+'_'+key+'_'+str(g),a,current,g)
    for key,series in result['metricResiduals']['parts'].items():
        for g,a in series.items():tensor('FIELD_CURRENT_METRIC_RESIDUAL_'+key+'_'+str(g),a,current,g,literal=True)
    for g,a in result['metricResiduals']['phaseInverse'].items():tensor('PHASE_INVERSE_RESIDUAL_'+str(g),a,degree=g,literal=True)
    for name,record in result['open'].items():
        for component,series in record['components'].items():
            for g,a in series.items():tensor(name+'_AMPLITUDE_'+component+'_'+str(g),a,degree=g,epsilon=1,literal=component.endswith('Residual'))
        for part,series in record['currentParts'].items():
            for g,a in series.items():tensor(name+'_CURRENT_'+part+'_'+str(g),a,current,g,2)
        for part,series in record['homotopyParts'].items():
            for power,a in series.items():tensor(name+'_HOMOTOPY_CURRENT_'+part+'_'+str(power),a,current,epsilon=2,homotopy=power)
        for key in ('totalFraction','inducedFraction'):
            for power,a in record[key].items():tensor(name+'_'+key+'_'+str(power),a,literal=True,homotopy=power)
        for key in ('inducedCurrent','totalMinusBaselineVersusInduced'):
            for g,a in record[key].items():tensor(name+'_'+key+'_'+str(g),a,current,g,2)
        tensor(name+'_WEAK_AMPLITUDE_DERIVATIVE',record['weakAmplitude'],literal=True)
        tensor(name+'_WEAK_INDUCED_QUADRATIC_COEFFICIENT',record['weakInducedQuadratic'],literal=True)
        for key in ('totalQuotient','inducedQuotient'):
            for power,a in record['residuals'][key].items():tensor(name+'_RESIDUAL_'+key+'_'+str(power),a,current,literal=True,homotopy=power)
        for part,series in record['residuals']['directContractions'].items():
            for g,a in series.items():tensor(name+'_RESIDUAL_CONTRACTION_'+part+'_'+str(g),a,current,g,2,literal=True)
        for g,a in record['residuals']['currentPartJoin'].items():tensor(name+'_RESIDUAL_CURRENT_PART_JOIN_'+str(g),a,current,g,2,literal=True)
    for end,series in result['closedAmplitudes'].items():
        for component,values in series.items():
            for g,a in values.items():tensor(end+'_CLOSED_MATCHING_'+component+'_'+str(g),a,degree=g,epsilon=1,literal=component.endswith('Residual'))
    for end,record in result['endCurrents'].items():
        for part,values in record['parts'].items():
            for key,series in values.items():
                for g,a in series.items():tensor(end+'_NORMAL_'+part+'_'+key+'_'+str(g),a,current,g,2)
        for part,values in record['residuals'].items():
            for key,series in values.items():
                for g,a in series.items():tensor(end+'_RESIDUAL_'+part+'_'+key+'_'+str(g),a,current,g,2,literal=True)
    for g,a in result['currentMutation'].items():tensor('CURRENT_COEFFICIENT_MUTATION_'+str(g),a,current,g,2)
    for end,rows in result['bulkDomains'].items():
        for i,row in enumerate(rows):
            for key,unit in [('wave',(0,-2,0)),('depthSquared',(-2,0,0)),('criticalMomentum',(-1,0,0)),
                             ('globalMaximum',(-2,0,0)),('curvature',zero),('reconstructionResidual',(0,-2,0))]:
                tag=PREFIX+'_'+end+'_BULK_'+str(i)+'_'+key;body=row[key]
                engine.emit(tag,body);engine.emit('METADATA_'+tag,mode.numeric_metadata(body,lambda p,u=unit:u))
    boundary.structural_flags(PREFIX+'_MANIFEST',{'sourceFiles':result['sourceFiles'],'inputPackets':result['inputPackets'],
        'census':result['selection']['census'],'openRanks':result['selection']['ranks'],
        'bulkDomains':{e:[{'realPropagationSet':str(v['realPropagationSet']),'domain':v['domain']} for v in rows] for e,rows in result['bulkDomains'].items()},
        'homotopy':{'lambda':'eta_bg','sigmaOverLambda':str(result['ratio']),'weakCoefficientOrder':'lambda powers stripped after differentiation; generating coefficients emitted separately'},
        'strongEdgeHandoff':'C_strong(1) remains a downstream obligation; not evaluated here.',
        'scope':result['scope']})


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True);args=parser.parse_args()
    base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));started=time.monotonic()
    def timeout(*_):raise TimeoutError('channel current budget; preserve completed packets')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(900)
    r,source,ends,reference,pencils,pins,operands=load(base)
    selection=selectors(source,ends,reference);f.atomic_pickle(base/'channel-selectors.pickle',selection)
    metrics,metric_residuals=open_metrics(source,ends,selection);opened,mutation=construct_open(source,metrics,selection)
    f.atomic_pickle(base/'open-channel-currents.pickle',{'metrics':metrics,'records':opened,'residuals':metric_residuals,'mutation':mutation})
    closed={}
    for end in ('LEFT','RIGHT'):
        indices=[i for row in selection['census'][end] if row['kind']=='evanescent' and row['direction']=='outgoing' for i in row['columns']]
        closed[end]=amplitude_bookkeeping({g:a[indices] for g,a in source['response']['outgoingFieldAmplitudes'][end].items()})
    finite=end_currents(source,ends);f.atomic_pickle(base/'finite-end-currents.pickle',finite)
    bulk=bulk_domains(pencils);f.atomic_pickle(base/'bulk-kinematics.pickle',bulk)
    result={'selection':selection,'metrics':metrics,'metricResiduals':metric_residuals,'open':opened,'closedAmplitudes':closed,
        'endCurrents':finite,'bulkDomains':bulk,'currentMutation':mutation,'currentUnit':source['currentUnit'],'ratio':source['ratio'],
        'sourceFiles':pins,'inputPackets':operands,'settings':source['settings'],'dimensionState':dict(vars(engine.PHYSICAL_METADATA.dimensions)),
        'scope':'Finite retained-rectangle open-channel and closed matching response; full finite-end normal current is distinct from bulk-depth radiation. No physical gain/loss inferred from omitted pure second-order terms; no pole solve.'}
    f.atomic_pickle(base/'continuum-currents.pickle',result);before=f.digest(base/'continuum-currents.pickle')
    engine.EMISSION_LINES.clear();engine.PAYLOAD_ENCODER=grades.PayloadEncoder()
    with (base/'full.out').open('x') as stream,contextlib.redirect_stdout(stream):
        emit_result(result,r)
        keys={tag:'s11cdContinuumCurrents'+str(i) for i,tag in enumerate(engine.EMISSION_LINES) if not tag.startswith('PY_S11CD_METADATA_')}
        boundary.structural_flags(PREFIX+'_WRITE_KEYS',keys);index=engine.emission_index(engine.EMISSION_LINES);boundary.structural_flags(PREFIX+'_EMISSION_LINES',index)
    entries={}
    for line in grades.decoded_lines(base/'full.out'):
        tag,_,body=line.rstrip('\n').partition(': ');f.require(tag not in entries,'unique current tag');entries[tag]=grades._restore(body)
    original=engine.emit;seen=set()
    def replay(name,value):
        tag='PY_S11CD_'+name;f.require(tag not in seen and entries.get(tag)==engine.cas(value),('current emission replay',tag));seen.add(tag)
    engine.emit=replay
    try:
        emit_result(result,r);boundary.structural_flags(PREFIX+'_WRITE_KEYS',keys);boundary.structural_flags(PREFIX+'_EMISSION_LINES',index)
    finally:engine.emit=original
    f.require(seen==set(entries) and len(keys)==len(set(keys.values())) and not set(keys.values())&set(engine.IMPORT_KEYS),'current key/payload census')
    paths=0
    for tag,body in entries.items():
        if not tag.startswith('PY_S11CD_METADATA_'):continue
        structural=tag.endswith(('_MANIFEST','_WRITE_KEYS','_EMISSION_LINES'))
        for item in body:
            fields={str(k):v for k,v in (item[1] if structural else item)}
            f.require(len(fields['DIMENSION_L_T_M'])==3 and all(not v.free_symbols for v in fields['DIMENSION_L_T_M']),'resolved current dimension')
            f.require('MULTIGRADE' in fields and 'EPSILON_LAMBDA_SUPPORT' in fields,'current grade paths');paths+=1 if structural else len(fields['PATHS'])
    final='PY_S11CD_'+PREFIX+'_EMISSION_LINES';grades.restore_emission_index({str(k):v for k,v in entries[final]},list(entries)[:list(entries).index(final)])
    f.require(before==f.digest(base/'continuum-currents.pickle') and not engine.PHYSICAL_METADATA.dimensions.constraints,'current packet and dimension closure')
    f.require(all(f.digest(f.ROOT/n)==h for n,h in pins.items()) and all(f.digest(Path(n))==h for n,h in operands.items()),'current source/input post hashes')
    summary={'runDirectory':str(base),'sourceFiles':pins,'inputPackets':operands,'openRanks':selection['ranks'],
        'bulkRealPropagationSets':{e:[str(v['realPropagationSet']) for v in rows] for e,rows in bulk.items()},
        'bulkMaximumDepthSquared':{e:[str(v['globalMaximum']) for v in rows] for e,rows in bulk.items()},
        'openResidualMaximum':max(boundary.norm(v['residuals']) for v in opened.values()),
        'finiteEndResidualMaximum':max(boundary.norm(v['residuals']) for v in finite.values()),
        'currentMutationMaximum':boundary.norm(mutation),'tagCount':len(entries),'writeKeys':len(keys),'metadataPaths':paths,
        'packetSha256BeforeEmission':before,'packetSha256AfterEmission':f.digest(base/'continuum-currents.pickle'),
        'artifacts':{p.name:{'bytes':p.stat().st_size,'sha256':f.digest(p)} for p in base.iterdir() if p.suffix in ('.pickle','.out')},
        'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':result['scope']}
    f.save(base/'checks.json',summary);signal.alarm(0);print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
