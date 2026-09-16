#!/usr/bin/env python3
"""Source-derived momentum boxes and their literal Fourier quadrature operands."""
import argparse,copy,json,resource,shutil,time,types
from pathlib import Path
import numpy as np
import sympy as sp
import S11c_d_position_domain as position

parallel,native,momentum,engine=position.parallel,position.native,position.momentum,position.engine
ROOT,STORE=position.ROOT,position.STORE
require,digest,save,atomic_pickle,unpickle=position.require,position.digest,position.save,position.atomic_pickle,position.unpickle
M=ROOT/'_measurements';CHECKPOINT=M/'S11c_d_position_domain_checkpoint.json'
PLAN=M/'S11c_d_momentum_domain_plan.md';PREFIX='MOMENTUM_DOMAIN_PREPARE_LAB_HELD_RHO4_CONSTANT'
CUTOFFS=(2.,3.,4.)


def artifact(base,path):return position.artifact(base,path)


def assignments(frequency,momenta,cutoff,extra=()):
    box={k:(-sp.Rational(str(cutoff)),sp.Rational(str(cutoff))) for k in momenta}
    info=engine.BoundedSourceFourierQuadrature.affine_range(frequency,box)
    lo,hi=map(float,info['bounds']);values=sorted(set([*np.linspace(lo,hi,65),*(float(v) for v in (0.,*extra) if lo<=float(v)<=hi)]))
    points=[]
    for value in values:
        t=0. if hi==lo else (value-lo)/(hi-lo)
        points.append([(-cutoff+2*cutoff*t if info['coefficients'][k]>0 else cutoff-2*cutoff*t if info['coefficients'][k]<0 else 0.) for k in momenta])
    points=np.asarray(points);fn=sp.lambdify(momenta,frequency,'numpy')
    residual=np.asarray([float(fn(*v))-nu for v,nu in zip(points,values)])
    require(np.max(abs(residual))<1e-12,'affine frequency assignment')
    return {'box':box,'range':info,'frequencies':np.asarray(values),'assignments':points,'assignmentResidual':residual}


def rebind(r,bound,cutoff):
    momenta=tuple(r.normal_map[g[2]] for g in r.momentum_groups);cut=sp.Rational(str(cutoff))
    symbols={limit[2] for row in bound['rows'] for limit in row['symbolicLimits']};mapping=dict(bound['cutoffBindings'])
    require(len(symbols)==3 and all(mapping[s]==2 for s in symbols),'native momentum cutoff symbols')
    for s in symbols:mapping[s]=cut
    rows=[];row_joins=[]
    for old in bound['rows']:
        limits=tuple(engine.memo_xreplace(v,mapping) for v in old['symbolicLimits'])
        restored=tuple(engine.memo_xreplace(v,bound['cutoffBindings']) for v in old['symbolicLimits'])
        require(restored==old['limits'] and all(tuple(v[1:])==(-cut,cut) for v in limits),'native limit/reverse join')
        require(tuple(v[0] for v in limits)==tuple(v[0] for v in old['limits']),'ordered native momentum variables')
        rows.append(dict(old,limits=limits));row_joins.append({'row':old['index'],'old':old['limits'],'new':limits,'restored':restored,'sourceLimit':old['sourceLimit']})
    sources={};source_joins=[]
    for key,old in bound['sources'].items():
        info=assignments(old['frequency'],momenta,cutoff,[bound['tests'][old['test']][1]])
        fresh=dict(old,**{k:v for k,v in info.items() if k!='box'})
        sources[key]=fresh;source_joins.append({'test':key[0],'sourceIndex':key[1],'uses':old['uses'],'frequency':old['frequency'],**info,'unit':bound['momentumUnit']})
        require(fresh['boundSource']==old['boundSource'] and fresh['boundAmplitude']==old['boundAmplitude'] and fresh['originalSourceIntegral']==old['originalSourceIntegral'],'unchanged source integrand')
    profiles=[]
    for j,p in enumerate(bound['profiles']):
        integral=next(v for v in bound['profileUnits'] if v.function==p['bound'].function)
        require(integral.limits==((r.xi,-14,14),),'accepted finite profile domain')
        info=assignments(p['transfer'],momenta,cutoff)
        amplitude=engine.memo_xreplace(p['normalizedIntegrand'],{r.transfer:sp.S.Zero})
        residual=sp.simplify(engine.memo_xreplace(p['normalizedIntegrand'],{r.transfer:p['transfer']})-integral.function)
        require(residual==0,'native profile transfer reconstruction')
        profiles.append({'index':j,'integral':integral,'source':p['source'],'transfer':p['transfer'],'normalized':p['normalizedIntegrand'],'amplitude':amplitude,'reconstructionResidual':residual,'unit':p['unit'],**info})
    actual=set().union(*(f['coefficient'].atoms(sp.Integral) for row in rows for f in row['factors']))
    require(actual==set(bound['profileUnits']),'all six profiles retained')
    return dict(bound,rows=rows,sources=sources,cutoffBindings=mapping),{'cutoff':cutoff,'momenta':momenta,'rows':row_joins,'sources':source_joins,'profiles':profiles,'cutoffBindings':mapping}


def load(base):
    accepted,previous=momentum.source.accepted(CHECKPOINT)
    for n,h in accepted['sourceFiles'].items():require(digest(ROOT/n)==h==digest(previous/'source'/n),'accepted source')
    for key in ('recordArtifacts','workerArtifacts'):
        for a in accepted[key]:require(digest(previous/a['path'])==a['sha256'],'accepted layout/partial/worker')
    bp=unpickle(previous/'bound-position-domains.pickle');packet=unpickle(previous/'position-domain.pickle')
    require(bp['provenance']==packet['provenance']==accepted['provenance'],'accepted position packet joins')
    _,sb=momentum.source.accepted(momentum.source.SOURCE);r,dimensions=momentum.source.native.source.restore_context(unpickle(sb/'reduced-action.pickle'));dimensions.__dict__.update(packet['dimensionState'])
    bound=bp['domains'][2];domains={};joins={}
    for i,cut in enumerate(CUTOFFS):domains[i],joins[i]=rebind(r,bound,cut)
    require(domains[0]['rows']==bound['rows'] and domains[0]['cutoffBindings']==bound['cutoffBindings'],'baseline native binding')
    for key,old in bound['sources'].items():
        v=domains[0]['sources'][key];require(v['range']==old['range'] and np.array_equal(v['frequencies'],old['frequencies']) and np.array_equal(v['assignments'],old['assignments']),'baseline source-range identity')
    pins=dict(accepted['sourceFiles'])
    for p in (CHECKPOINT,PLAN,Path(__file__).resolve()):pins[str(p.relative_to(ROOT))]=digest(p)
    for n in pins:
        p=base/'source'/n;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/n,p)
    for old,new in (('bound-position-domains.pickle','accepted-bound-position-domains.pickle'),('position-domain.pickle','accepted-position-domain.pickle')):shutil.copyfile(previous/old,base/new)
    provenance={'POSITION_DOMAIN_CHECKPOINT_SHA256':digest(CHECKPOINT),'POSITION_DOMAIN_BOUND_PACKET_SHA256':digest(previous/'bound-position-domains.pickle'),'POSITION_DOMAIN_PACKET_SHA256':digest(previous/'position-domain.pickle'),'APPROVED_INPUT_SHA256':accepted['provenance']['APPROVED_INPUT_SHA256'],'INSTRUMENT_SHA256':digest(Path(__file__)),'PLAN_SHA256':digest(PLAN)}
    data={'r':r,'bound':bound,'domains':domains,'joins':joins,'provenance':provenance,'sourceFiles':pins}
    atomic_pickle(base/'momentum-domains.pickle',{'domains':domains,'joins':joins,'provenance':provenance,'sourceFiles':pins,'dimensionState':dict(vars(dimensions))})
    save(base/'preflight.json',{'sourceFiles':pins,'provenance':provenance,'momentumCutoffs':CUTOFFS,'sourceCutoff':48.,'profileCutoff':14.,'regulator':.2,'rows':80,'sources':70,'profiles':6,'baselineRangeAndLimitJoins':True})
    return data


def record_guard(record):
    scale=1+abs(record['adaptiveValues'])
    require(np.all(np.isfinite(record['gaussValues'])) and np.all(np.isfinite(record['adaptiveValues'])),'finite transforms')
    require(np.max(abs(record['literalResidual'])/scale)<1e-9,'literal original transform residual')
    require(np.max(abs(record['adaptiveResidual'])/scale)<1e-9,'finest transform/adaptive residual')
    require(np.max(abs(record['measureMutationResidual']))>0,'actual transform measure response')


def evaluate_task(directory,data,task,*,max_batches=None,resume_state=None):
    require(max_batches is None and resume_state is None,'preparation task')
    ti,index=task;require(ti in (0,1) and index in (1,2),'preparation assignment')
    r=data['r'];bound=data['domains'][index];join=data['joins'][index];momenta=join['momenta'];records=[];inventory=[];started=time.monotonic()
    for key,source in sorted(bound['sources'].items()):
        if key[0]!=ti:continue
        width=float(source['testWidth']);ell=source['profileWidth'];points=sorted({-48.,48.,0.,-width,width,-ell,ell});orders=(128,256,384)
        q=engine.BoundedSourceFourierQuadrature(r.zp,source['boundAmplitude']);freq=source['frequencies']
        values=np.asarray([q.gauss(freq,points,n) for n in orders]);adaptive_values,errors=q.adaptive(freq,points)
        x,w=q.rule(points,orders[-1]);fn=sp.lambdify((r.zp,*momenta),source['boundSource'],'numpy',cse=True)
        direct=np.asarray([np.dot(np.broadcast_to(np.asarray(fn(x,*p),complex),x.shape),w) for p in source['assignments']])
        mutated=q.gauss(freq,points,orders[-1],weight_scale=1.001)
        v={'kind':'source','test':ti,'index':index,'sourceIndex':key[1],'orders':orders,'points':points,'frequencies':freq,'assignments':source['assignments'],'original':source['boundSource'],'originalUnit':source['amplitudeUnit'],'unit':source['integralUnit'],
            'gaussValues':values,'adaptiveValues':adaptive_values,'adaptiveErrorEstimates':errors,'literalValues':direct,'literalResidual':values[-1]-direct,'adaptiveResidual':values[-1]-adaptive_values,'orderChanges':np.diff(values,axis=0),'measureMutationValues':mutated,'measureMutationResidual':mutated-values[-1],'peakWorkspaceEstimateBytes':q.peak_workspace_bytes}
        path=directory/'records'/f'source-{key[1]:02}.pickle';path.parent.mkdir(exist_ok=True);atomic_pickle(path,v);inventory.append(artifact(directory,path));save(directory/'record-inventory.json',inventory);record_guard(v);records.append(v)
    if ti==0:
        for p in join['profiles']:
            orders=(256,384,512);points=(-14.,0.,14.);q=engine.BoundedSourceFourierQuadrature(r.xi,p['amplitude']);freq=p['frequencies']
            values=np.asarray([q.gauss(freq,points,n) for n in orders]);adaptive_values,errors=q.adaptive(freq,points)
            x,w=q.rule(points,orders[-1]);fn=sp.lambdify((r.xi,*momenta),p['integral'].function,'numpy',cse=True)
            direct=np.asarray([np.dot(np.broadcast_to(np.asarray(fn(x,*a),complex),x.shape),w) for a in p['assignments']]);mutated=q.gauss(freq,points,orders[-1],weight_scale=1.001)
            v={'kind':'profile','test':ti,'index':index,'profileIndex':p['index'],'orders':orders,'points':points,'frequencies':freq,'assignments':p['assignments'],'original':p['integral'],'originalUnit':p['unit'],'unit':p['unit'],
                'gaussValues':values,'adaptiveValues':adaptive_values,'adaptiveErrorEstimates':errors,'literalValues':direct,'literalResidual':values[-1]-direct,'adaptiveResidual':values[-1]-adaptive_values,'orderChanges':np.diff(values,axis=0),'measureMutationValues':mutated,'measureMutationResidual':mutated-values[-1],'peakWorkspaceEstimateBytes':q.peak_workspace_bytes}
            path=directory/'records'/f'profile-{p["index"]}.pickle';atomic_pickle(path,v);inventory.append(artifact(directory,path));save(directory/'record-inventory.json',inventory);record_guard(v);records.append(v)
    # Explicit corners, center, transfer diagonals and width offsets, not a
    # random witness for the enlarged box. This is only sampled execution.
    cut=CUTOFFS[index];width=float(bound['abel']['width'].subs(r.regulator,.2));coordinates={tuple(v) for v in __import__('itertools').product((-cut,0.,cut),repeat=3)}
    for a,b in bound['pairs']:
        for center in (-cut/2,0.,cut/2):
            for offset in (-width,0.,width):
                point=[0.]*3;point[momenta.index(a)]=center;point[momenta.index(b)]=center+offset;coordinates.add(tuple(point))
    coordinates=np.asarray(sorted(coordinates));env={k:coordinates[:,j] for j,k in enumerate(momenta)};env[r.regulator]=np.full(len(coordinates),.2)
    worker=engine.BoundedSourceFourierQuadrature.ThreeMomentum(bound['rows'],bound['sources'],r);setting={'kind':'concentrated','sourceBound':48.,'sourceNodes':384,'profileBound':14.,'profileNodes':512,'regulator':.2,'momentumBound':cut};cache={};coefficients=[]
    for row in bound['rows']:
        for j,f in enumerate(row['factors']):coefficients.append({'row':row['index'],'factor':j,'sourceIndex':f['sourceIndex'],'unit':f['unit'],'values':worker.coefficient_value(f['coefficient'],env,np.asarray(bound['positions']),setting,cache)})
    probe={'momenta':momenta,'coordinates':coordinates,'pairs':bound['pairs'],'width':width,'setting':setting,'coefficients':coefficients,'positions':bound['positions'],'evaluatedProfiles':worker.profile_integrals}
    atomic_pickle(directory/'coefficient-probes.pickle',probe);inventory.append(artifact(directory,directory/'coefficient-probes.pickle'))
    require(worker.profile_integrals==set(bound['profileUnits']) and len(coefficients)==80,'all native coefficient/profile probes')
    packet={'task':task,'records':records,'coefficientProbes':probe,'recordArtifacts':inventory,'sourceFiles':data['sourceFiles'],'provenance':data['provenance'],'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    atomic_pickle(directory/'result.pickle',packet);require(data['sourceFiles']=={n:digest(ROOT/n) for n in data['sourceFiles']},'worker source identity')
    save(directory/'checks.json',{'task':task,'resultSha256':digest(directory/'result.pickle'),'wallSeconds':packet['wallSeconds'],'peakRssKiB':packet['peakRssKiB']});return packet


def dispatch(base,data):
    old=parallel.evaluate_task;parallel.evaluate_task=evaluate_task
    try:return parallel.dispatch(base,data,[(ti,index) for ti in range(2) for index in (1,2)])
    finally:parallel.evaluate_task=old


def emit(result,r,bound,rows,finest,provenance):
    dimensions=engine.PHYSICAL_METADATA.dimensions;zero=dimensions.zero;kun=bound['momentumUnit'];metadata=engine.FullPencilModes.__new__(engine.FullPencilModes);metadata.r=r
    def numeric(name,value,unit,literal=False):
        body=momentum.number(value);engine.emit(PREFIX+'_'+name,body if literal else metadata.compact_fingerprint(body));engine.emit('METADATA_'+PREFIX+'_'+name,metadata.numeric_metadata(body,unit))
    engine.physical(PREFIX+'_PROVENANCE',provenance);engine.physical(PREFIX+'_ABEL_PAIRS',bound['pairs']);engine.physical(PREFIX+'_ABEL_WIDTH',bound['abel']['width'])
    numeric('SOURCE_CUTOFF',48.,lambda p:dimensions.measure(r.zp),True);numeric('PROFILE_CUTOFF',14.,lambda p:dimensions.measure(r.xi),True)
    for index,join in result['joins'].items():
        numeric('MOMENTUM_CUTOFF_'+str(index),join['cutoff'],lambda p:kun,True)
        for row in join['rows']:
            for name in ('old','new','restored'):numeric('LIMIT_'+name.upper()+'_'+str(index)+'_'+str(row['row']),row[name],lambda p:kun,True)
        for kind in ('sources','profiles'):
            for j,v in enumerate(join[kind]):
                suffix=str(index)+'_'+kind.upper()+'_'+str(j);freq_unit=kun if kind=='sources' else zero;coeff_unit=tuple(a-b for a,b in zip(freq_unit,kun))
                if kind=='sources':
                    numeric('SOURCE_IDENTITY_'+suffix,(v['test'],v['sourceIndex']),lambda p:zero,True)
                    numeric('SOURCE_ROW_USES_'+suffix,v['uses'],lambda p:zero,True)
                else:
                    numeric('PROFILE_IDENTITY_'+suffix,v['index'],lambda p:zero,True)
                    engine.emit(PREFIX+'_PROFILE_SOURCE_'+suffix,engine.carrier_fingerprint(v['source']))
                    engine.emit('METADATA_'+PREFIX+'_PROFILE_SOURCE_'+suffix,metadata.numeric_metadata(v['source'],lambda p:v['unit']))
                    numeric('PROFILE_RECONSTRUCTION_RESIDUAL_'+suffix,v['reconstructionResidual'],lambda p:v['unit'],True)
                frequency=v['frequency'] if kind=='sources' else v['transfer']
                engine.emit(PREFIX+'_FREQUENCY_EXPRESSION_'+suffix,engine.carrier_fingerprint(frequency))
                engine.emit('METADATA_'+PREFIX+'_FREQUENCY_EXPRESSION_'+suffix,metadata.numeric_metadata(frequency,lambda p:freq_unit))
                numeric('RANGE_ORIGIN_'+suffix,v['range']['origin'],lambda p:freq_unit,True);numeric('RANGE_COEFFICIENTS_'+suffix,list(v['range']['coefficients'].values()),lambda p:coeff_unit,True)
                numeric('RANGE_BOUNDARIES_'+suffix,v['range']['bounds'],lambda p:freq_unit,True);numeric('RANGE_RESIDUAL_'+suffix,v['range']['residual'],lambda p:freq_unit,True)
                numeric('ASSIGNMENTS_'+suffix,v['assignments'],lambda p:kun);numeric('FREQUENCIES_'+suffix,v['frequencies'],lambda p:freq_unit,True);numeric('ASSIGNMENT_RESIDUAL_'+suffix,v['assignmentResidual'],lambda p:freq_unit,True)
    for item in result['workers']:
        ti,index=item['task']
        for v in item['records']:
            suffix=f'{ti}_{index}_{v["kind"].upper()}_{v.get("sourceIndex",v.get("profileIndex"))}';coordinate_unit=dimensions.measure(r.zp if v['kind']=='source' else r.xi);freq_unit=tuple(-x for x in coordinate_unit)
            engine.emit(PREFIX+'_ORIGINAL_'+suffix,engine.carrier_fingerprint(v['original']));engine.emit('METADATA_'+PREFIX+'_ORIGINAL_'+suffix,metadata.numeric_metadata(v['original'],lambda p:v['originalUnit']))
            numeric('ORDERS_'+suffix,v['orders'],lambda p:zero,True);numeric('PANELS_'+suffix,v['points'],lambda p:coordinate_unit,True);numeric('SAMPLED_FREQUENCIES_'+suffix,v['frequencies'],lambda p:freq_unit,True)
            for name in ('gaussValues','adaptiveValues','adaptiveErrorEstimates','literalValues','literalResidual','adaptiveResidual','orderChanges','measureMutationValues','measureMutationResidual'):numeric(name.upper()+'_'+suffix,v[name],lambda p:v['unit'],name.endswith('Residual') or name in ('orderChanges','adaptiveErrorEstimates'))
        probe=item['coefficientProbes'];suffix=f'{ti}_{index}'
        numeric('PROBE_MOMENTA_'+suffix,probe['coordinates'],lambda p:kun);numeric('PROBE_POSITIONS_'+suffix,probe['positions'],lambda p:dimensions.measure(r.z),True)
        for v in probe['coefficients']:numeric('COEFFICIENT_VALUES_'+suffix+'_'+str(v['row'])+'_'+str(v['factor']),v['values'],lambda p:v['unit'])


def main():
    p=argparse.ArgumentParser();p.add_argument('--run-directory',type=Path,required=True);a=p.parse_args();base=a.run_directory.resolve();base.relative_to(STORE);base.mkdir(parents=True,exist_ok=False);started=time.monotonic();data=load(base)
    workers=dispatch(base,data);artifacts=[]
    for task,w in sorted(workers.items()):
        require(w['sourceFiles']==data['sourceFiles'] and w['provenance']==data['provenance'],'worker source/provenance')
        directory=base/f'worker-{task[0]}-{task[1]}';artifacts.append(artifact(base,directory/'result.pickle'))
        for v in w['recordArtifacts']:
            require(digest(directory/v['path'])==v['sha256'],'worker record hash');artifacts.append(dict(v,path=str((directory/v['path']).relative_to(base))))
    result={'joins':data['joins'],'workers':[w for _,w in sorted(workers.items())]};atomic_pickle(base/'momentum-domain-preparation.pickle',{'result':result,'provenance':data['provenance'],'sourceFiles':data['sourceFiles'],'dimensionState':dict(vars(engine.PHYSICAL_METADATA.dimensions))});before={p.name:digest(p) for p in base.glob('*.pickle')}
    # No cache tags use the position-domain address schema. Reuse the unchanged
    # native manifest, payload and metadata validator with this emitter/prefix.
    replay=types.FunctionType(native.emit_and_replay.__code__,dict(native.emit_and_replay.__globals__,emit=emit,PREFIX=PREFIX),native.emit_and_replay.__name__,native.emit_and_replay.__defaults__,native.emit_and_replay.__closure__)
    engine.EMISSION_LINES.clear();engine.PAYLOAD_ENCODER=engine.PayloadEncoder();entries,keys,paths=replay(base,result,data['r'],data['bound'],data['bound']['rows'],{},data['provenance'])
    norms=[{'task':w['task'],'kind':v['kind'],'sourceIndex':v.get('sourceIndex'),'profileIndex':v.get('profileIndex'),'literalResidual':float(np.max(abs(v['literalResidual']))),'adaptiveResidual':float(np.max(abs(v['adaptiveResidual']))),'orderChanges':[float(np.max(abs(row))) for row in v['orderChanges']],'adaptiveErrorEstimate':float(np.max(v['adaptiveErrorEstimates']))} for w in result['workers'] for v in w['records']]
    summary={'runDirectory':str(base),'sourceFiles':data['sourceFiles'],'provenance':data['provenance'],'rows':80,'sources':70,'profiles':6,'domainChoices':CUTOFFS,'transformRecords':len(norms),'norms':norms,'workerArtifacts':artifacts,'workerManifest':json.loads((base/'workers.json').read_text()),'tagCount':len(entries),'writeKeyCount':len(keys),'metadataPaths':paths,'packetHashesBeforeEmission':before,'packetHashesAfterEmission':{n:digest(base/n) for n in before},'artifacts':{p.name:{'bytes':p.stat().st_size,'sha256':digest(p)} for p in base.iterdir() if p.suffix in ('.pickle','.out')},'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':'Momentum cutoff binding and sampled source/profile quadrature and coefficient execution only. No wider-box complete action, uniform or exceptional-locus coverage, infinite-tail bound, Abel limit, scattering or pole.'}
    require(before==summary['packetHashesAfterEmission'] and not engine.PHYSICAL_METADATA.dimensions.constraints and data['sourceFiles']=={n:digest(ROOT/n) for n in data['sourceFiles']},'final packet/source/dimension guard')
    for v in artifacts:require(digest(base/v['path'])==v['sha256'],'final record hash')
    save(base/'checks.json',summary);print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
