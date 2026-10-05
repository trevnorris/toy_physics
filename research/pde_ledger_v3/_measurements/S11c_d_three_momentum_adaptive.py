#!/usr/bin/env python3
"""Source-conditioned native inner sums with independent adaptive outer integration."""
import argparse,ast,copy,hashlib,inspect,json,resource,shutil,textwrap,time
from pathlib import Path
import numpy as np
from scipy.integrate import quad_vec
import S11c_d_parallel_momentum as parallel

native=parallel.native
engine=parallel.engine
ROOT,STORE=parallel.ROOT,parallel.STORE
require,digest,save,atomic_pickle,unpickle=parallel.require,parallel.digest,parallel.save,parallel.atomic_pickle,parallel.unpickle
M=ROOT/'_measurements'
CHECKPOINT=M/'S11c_d_three_momentum_source_checkpoint.json'
PLAN=M/'S11c_d_three_momentum_adaptive_plan.md'
PREFIX='THREE_MOMENTUM_ADAPTIVE_LAB_HELD_RHO4_CONSTANT'


def derived_methods():
    old=ast.parse(textwrap.dedent(inspect.getsource(engine.BoundedSourceFourierQuadrature.ThreeMomentum.batches)))
    tree=copy.deepcopy(old)
    loop=next(n for n in tree.body[0].body if isinstance(n,ast.For))
    require(isinstance(loop.iter,ast.Call) and getattr(loop.iter.func,'id',None)=='descend','native recursion entry')
    original_args=copy.deepcopy(loop.iter.args)
    loop.iter.args=ast.parse('descend(len(variables)-2,{variables[-1]:self.fixed_outer},1.)',mode='eval').body.args
    proof=copy.deepcopy(tree);next(n for n in proof.body[0].body if isinstance(n,ast.For)).iter.args=original_args
    require(ast.dump(proof)==ast.dump(old),'conditional recursion AST join')
    ns=dict(vars(engine));exec(compile(ast.fix_missing_locations(tree),'<native-conditioned-recursion>','exec'),ns)
    batches=ns['batches'];batch_hash=hashlib.sha256(ast.dump(old).encode()).hexdigest()
    old=ast.parse(textwrap.dedent(inspect.getsource(engine.BoundedSourceFourierQuadrature.FiniteMomentum.group)))
    tree=copy.deepcopy(old)
    assign=next(n for n in tree.body[0].body if isinstance(n,ast.Assign) and getattr(n.targets[0],'id',None)=='expected_volume')
    prior=copy.deepcopy(assign.value);require(isinstance(prior,ast.BinOp) and isinstance(prior.op,ast.Pow),'native volume expression')
    assign.value.right=ast.BinOp(left=assign.value.right,op=ast.Sub(),right=ast.Constant(1))
    proof=copy.deepcopy(tree);next(n for n in proof.body[0].body if isinstance(n,ast.Assign) and getattr(n.targets[0],'id',None)=='expected_volume').value=prior
    require(ast.dump(proof)==ast.dump(old),'conditional group AST join')
    exec(compile(ast.fix_missing_locations(tree),'<native-conditioned-group>','exec'),ns)
    return batches,ns['group'],{'nativeBatchesAstSha256':batch_hash,'nativeGroupAstSha256':hashlib.sha256(ast.dump(old).encode()).hexdigest()}


BATCHES,GROUP,METHOD_JOINS=derived_methods()


class ConditionalMomentum(engine.BoundedSourceFourierQuadrature.ThreeMomentum):
    batches=BATCHES
    group=GROUP
    def __init__(self,*args,outer_value,**kwargs):
        super().__init__(*args,**kwargs);self.fixed_outer=float(outer_value)


def load(base):
    accepted,previous=native.momentum.source.accepted(CHECKPOINT)
    for n,h in accepted['sourceFiles'].items():
        require(digest(ROOT/n)==h==digest(previous/'source'/n),('accepted source',n))
    for key in ('recordArtifacts','partialArtifacts','workerArtifacts'):
        for a in accepted[key]:require(digest(previous/a['path'])==a['sha256'],'accepted saved operand')
    # This restores the unchanged source context and original accepted layouts.
    r,bound,rows,variables,finest,provenance,pins=native.load(base)
    packet=unpickle(previous/'three-momentum-source.pickle')
    require(packet['boundPacketSha256']==provenance['BOUND_PACKET_SHA256'],'accepted refinement bound join')
    for ti in range(2):
        item=next(v for v in packet['result']['records'] if v['test']==ti and v['index']==2)
        setting=item['setting']
        require((setting['outerOrder'],tuple(setting['innerOrders']),setting['sourceNodes'],setting['profileNodes'])==(144,(24,24),256,256),'accepted final setting')
        original=finest[ti]
        finest[ti]=dict(original,setting=setting,action=item['action'],contributions=[{'index':v['index'],'terms':v['values']} for v in item['terms']],
            groups=[item['group'] if len(g['variables'])==3 else g for g in original['groups']],layoutSettings=original['layoutSettings'].copy())
        finest[ti]['layoutSettings'][3]=setting
        require(all(tuple(l[0] for l in row['limits'])==variables for row in rows),'complete native limit join')
        require(all(tuple(map(float,l[1:]))==(-setting['momentumBound'],setting['momentumBound']) for row in rows for l in row['limits']),'native finite domain join')
    paths=[*(ROOT/n for n in accepted['sourceFiles']),CHECKPOINT,PLAN,Path(__file__).resolve(),M/'S11c_d_three_momentum_adaptive_preflight.py']
    for path in paths:
        n=str(path.relative_to(ROOT));pins[n]=digest(path);target=base/'source'/n;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,target)
    shutil.copyfile(previous/'three-momentum-source.pickle',base/'accepted-three-momentum-source.pickle')
    provenance=dict(provenance,SOURCE_REFINEMENT_CHECKPOINT_SHA256=digest(CHECKPOINT),SOURCE_REFINEMENT_PACKET_SHA256=digest(base/'accepted-three-momentum-source.pickle'),ADAPTIVE_INSTRUMENT_SHA256=digest(Path(__file__)),ADAPTIVE_PLAN_SHA256=digest(PLAN))
    data=(r,bound,rows,variables,finest,provenance,pins)
    save(base/'preflight.json',{'sourceFiles':pins,'provenance':provenance,'methodJoins':METHOD_JOINS,'acceptedRunDirectory':str(previous)})
    return data


def evaluate_task(directory,data,task,*,max_batches=None,resume_state=None):
    require(max_batches is None and resume_state is None,'adaptive worker request')
    r,bound,rows,variables,finest,provenance,pins=data;ti,half=task
    require(half in (0,1) and ti in (0,1),'adaptive task')
    setting=dict(finest[ti]['setting'],**finest[ti].get('adaptiveOverrides',{}));cut=setting['momentumBound'];interval=(-cut,0.) if half==0 else (0.,cut)
    tolerance=setting.get('adaptiveTolerance',5e-11)
    width=float(bound['abel']['width'].subs(r.regulator,setting['regulator']))
    worker=ConditionalMomentum(bound['rows'],bound['sources'],r,outer_value=0.)
    artifacts=[];point_artifacts=[];point_cache={};calls=0;cache_hits=0;started=time.monotonic();n=len(rows)*len(bound['positions'])
    def function(k):
        nonlocal calls,cache_hits
        calls+=1;key=float(k).hex()
        if key in point_cache:cache_hits+=1;return point_cache[key].copy()
        index=len(point_cache);worker.fixed_outer=float(k)
        def partial(state):
            if state['batchCount']%64==0:
                path=directory/'partials'/f'{index:05}-{state["batchCount"]:06}.pickle';path.parent.mkdir(exist_ok=True)
                atomic_pickle(path,{'outerValue':float(k),'state':state})
                artifacts.append({'path':str(path.relative_to(directory)),'bytes':path.stat().st_size,'sha256':digest(path)})
                save(directory/'partial-inventory.json',artifacts)
        group=worker.group(ti,variables,setting,bound['pairs'],width,bound['positions'],partial)
        packet={'task':task,'outerValue':float(k),'outerUnit':bound['momentumUnit'],
            'conditionalUnits':[tuple(x-y for x,y in zip(row['unit'],bound['momentumUnit'])) for row in rows],
            'conditionalMassUnit':tuple((len(variables)-1)*v for v in bound['momentumUnit']),
            'setting':setting,'group':group,'provenance':provenance,'sourceFiles':pins,'methodJoins':METHOD_JOINS}
        path=directory/'points'/f'{index:05}.pickle';path.parent.mkdir(exist_ok=True);atomic_pickle(path,packet)
        point_artifacts.append({'path':str(path.relative_to(directory)),'bytes':path.stat().st_size,'sha256':digest(path)})
        save(directory/'point-inventory.json',point_artifacts)
        require(abs(group['volumeResidual'])<1e-10*(1+abs(group['boxVolume'])),'conditional finite area')
        require(np.all(np.isfinite(group['values'])) and np.all(np.isfinite(group['measureMutationValues'])),'conditional finite values')
        value=np.concatenate((group['values'].ravel(),group['measureMutationValues'].ravel(),[group['quadratureMass']]))
        point_cache[key]=value
        return value.copy()
    value,error,info=quad_vec(function,*interval,epsabs=tolerance,epsrel=0.,norm='max',quadrature='gk21',workers=1,cache_size=8*1024*1024,limit=1000,full_output=True)
    shape=(len(rows),len(bound['positions']))
    group={'test':ti,'rowIndices':[row['index'] for row in rows],'variables':variables,'values':value[:n].reshape(shape),
        'measureMutationValues':value[n:2*n].reshape(shape),'measureMutationResidual':(value[n:2*n]-value[:n]).reshape(shape),
        'quadratureMass':value[-1],'boxVolume':(interval[1]-interval[0])*(2*cut)**(len(variables)-1),'interval':interval,
        'unitFrameErrorEstimate':float(error),'unitFrameTolerance':tolerance,'relativeTolerance':0.,'evaluations':int(info.neval),
        'success':bool(info.success),'status':int(info.status),'message':info.message,'intervals':info.intervals,
        'intervalPackedValues':info.integrals,'unitFrameIntervalErrors':info.errors,'calls':calls,'cacheHits':cache_hits,'uniquePoints':len(point_cache)}
    group['volumeResidual']=group['quadratureMass']-group['boxVolume']
    packet={'task':task,'kind':'adaptive-complete','setting':setting,'group':group,'telemetry':parallel.cache_state(worker),
        'pointArtifacts':point_artifacts,'partialArtifacts':artifacts,'provenance':provenance,'sourceFiles':pins,'methodJoins':METHOD_JOINS,
        'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    atomic_pickle(directory/'result.pickle',packet)
    require(group['success'] and np.isfinite(error) and error<=tolerance and np.all(np.isfinite(value)),'adaptive outer unfinished')
    require(abs(group['volumeResidual'])<1e-10*(1+abs(group['boxVolume'])),'adaptive finite mass')
    require(pins=={n:digest(ROOT/n) for n in pins},'adaptive worker source changed')
    save(directory/'checks.json',{'task':task,'kind':packet['kind'],'resultSha256':digest(directory/'result.pickle'),'wallSeconds':packet['wallSeconds'],'peakRssKiB':packet['peakRssKiB']})
    return packet


def dispatch(base,data):
    original=parallel.evaluate_task;parallel.evaluate_task=evaluate_task
    try:return parallel.dispatch(base,data,[(ti,half) for ti in range(2) for half in range(2)])
    finally:parallel.evaluate_task=original


def combine(base,data,packets):
    r,bound,rows,variables,finest,provenance,pins=data
    records=[];inventory=[];census={};profiles=set();cache_checks={};phase=batch=0;point_artifacts=[];partial_artifacts=[]
    for task,packet in sorted(packets.items()):
        require(packet['provenance']==provenance and packet['sourceFiles']==pins and packet['methodJoins']==METHOD_JOINS,'adaptive worker provenance')
        require(packet['setting']==dict(finest[task[0]]['setting'],**finest[task[0]].get('adaptiveOverrides',{})), 'adaptive worker setting join')
        require(tuple(packet['group']['variables'])==variables and packet['group']['rowIndices']==[row['index'] for row in rows], 'adaptive worker row/limit join')
        g=packet['group'];require(g['success'] and g['unitFrameErrorEstimate']<=g['unitFrameTolerance'],'adaptive status')
        for key,target in (('pointArtifacts',point_artifacts),('partialArtifacts',partial_artifacts)):
            for a in packet[key]:
                path=base/f'worker-{task[0]}-{task[1]}'/a['path'];require(digest(path)==a['sha256'],'adaptive saved artifact')
                target.append(dict(a,path=str(path.relative_to(base))))
        telemetry=packet['telemetry']
        profiles.update(telemetry['evaluatedProfileIntegrals'])
        for k,(lo,hi,count) in telemetry['sourceFrequencyCensus'].items():
            old=census.get(k,(float('inf'),float('-inf'),0));census[k]=(min(old[0],lo),max(old[1],hi),old[2]+count)
        for cache in telemetry['cacheChecks']:cache_checks[(tuple(cache['points']),cache['order'],cache['profileRule'])]=cache
        phase=max(phase,telemetry['phaseWorkspaceEstimateBytes']);batch=max(batch,telemetry['batchCacheEstimateBytes'])
    for ti in range(2):
        ref=next(g for g in finest[ti]['groups'] if g['variables']==variables)
        baseline={'test':ti,'index':0,'method':'retained','setting':finest[ti]['setting'],'group':ref,**native.full_action(r,bound,finest[ti],ref)}
        records.append(baseline)
        halves=[packets[(ti,h)]['group'] for h in range(2)]
        cut=finest[ti]['setting']['momentumBound']
        require([tuple(g['interval']) for g in halves]==[(-cut,0.),(0.,cut)],'outer half-domain partition')
        values=sum(g['values'] for g in halves);mutated=sum(g['measureMutationValues'] for g in halves);mass=sum(g['quadratureMass'] for g in halves)
        group={'test':ti,'variables':variables,'rowIndices':[v['index'] for v in rows],'values':values,'measureMutationValues':mutated,
            'measureMutationResidual':mutated-values,'quadratureMass':mass,'boxVolume':(2*cut)**len(variables),'volumeResidual':mass-(2*cut)**len(variables),
            'unitFrameErrorEstimate':sum(g['unitFrameErrorEstimate'] for g in halves),'unitFrameTolerance':sum(g['unitFrameTolerance'] for g in halves),
            'evaluations':sum(g['evaluations'] for g in halves),'success':all(g['success'] for g in halves),'halves':halves}
        item={'test':ti,'index':1,'method':'adaptive','setting':packets[(ti,0)]['setting'],'group':group,**native.full_action(r,bound,finest[ti],group),'referenceIntegralDifference':values-ref['values']}
        records.append(item)
        require(abs(group['volumeResidual'])<1e-10*(1+group['boxVolume']),'combined mass')
        require(not np.any((abs(values)>1e-9)&(abs(mutated-values)<1e-12)),'combined measure sensitivity')
    for item in records:
        for term in item['terms']:require(all(v==0 for changed,v in zip(term['changed'],term['differences']) if not changed),'adaptive held term')
        path=base/'records'/f'{item["test"]}-{item["index"]}.pickle';path.parent.mkdir(exist_ok=True);atomic_pickle(path,item)
        inventory.append({'path':str(path.relative_to(base)),'bytes':path.stat().st_size,'sha256':digest(path)})
    require(set(census)=={(ti,f['sourceIndex']) for ti in range(2) for row in rows for f in row['factors']},'adaptive source census')
    require(profiles==set().union(*(f['coefficient'].atoms(native.sp.Integral) for row in rows for f in row['factors'])),'adaptive profile census')
    result={'records':records,'recordArtifacts':inventory,'pointArtifacts':point_artifacts,'partialArtifacts':partial_artifacts,'cacheChecks':list(cache_checks.values()),
        'sourceFrequencyCensus':census,'evaluatedProfileIntegrals':tuple(sorted(profiles,key=native.sp.default_sort_key)),
        'sourceRuleCacheBytes':sum(v['bytes'] for v in cache_checks.values()),'phaseWorkspaceEstimateBytes':phase,'batchCacheEstimateBytes':batch,'workspaceBudgetBytes':32*1024*1024}
    for name,key in (('record-inventory.json','recordArtifacts'),('point-inventory.json','pointArtifacts'),('partial-inventory.json','partialArtifacts')):save(base/name,result[key])
    atomic_pickle(base/'three-momentum-adaptive.pickle',{'result':result,'provenance':provenance,'boundPacketSha256':digest(base/'accepted-bound-momentum.pickle'),'dimensionState':dict(vars(engine.PHYSICAL_METADATA.dimensions))})
    return result


def emit_and_replay(base,result,data):
    r,bound,rows,variables,finest,provenance,pins=data
    # Reuse the native emitter, replacing only its misleading Gauss outer-order label.
    tree=ast.parse(textwrap.dedent(inspect.getsource(native.emit)));old_tree=copy.deepcopy(tree)
    calls=[n for n in ast.walk(tree) if isinstance(n,ast.Call) and getattr(n.func,'id',None)=='numeric' and isinstance(n.args[0],ast.BinOp) and getattr(n.args[0].left,'value',None)=='CHANGED_RULE_ORDERS_']
    require(len(calls)==1,'native quadrature-order emission site')
    call=calls[0];old_args=copy.deepcopy(call.args)
    call.args[0].left.value='HELD_INNER_SOURCE_PROFILE_ORDERS_';call.args[1]=ast.parse("(*setting['innerOrders'],setting['sourceNodes'],setting['profileNodes'])",mode='eval').body
    proof=copy.deepcopy(tree);pc=next(n for n in ast.walk(proof) if isinstance(n,ast.Call) and getattr(n.func,'id',None)=='numeric' and isinstance(n.args[0],ast.BinOp) and getattr(n.args[0].left,'value',None)=='HELD_INNER_SOURCE_PROFILE_ORDERS_');pc.args=old_args
    require(ast.dump(proof)==ast.dump(old_tree),'native emitter order-label-only join')
    ns=dict(vars(native),PREFIX=PREFIX);exec(compile(ast.fix_missing_locations(tree),str(Path(__file__))+'<native-emitter>','exec'),ns);common=ns['emit']
    def emit(result,r,bound,rows,finest,provenance):
        common(result,r,bound,rows,finest,provenance)
        metadata=engine.FullPencilModes.__new__(engine.FullPencilModes);metadata.r=r
        def numeric(name,value,units):
            body=native.number(value);engine.emit(PREFIX+'_'+name,body);engine.emit('METADATA_'+PREFIX+'_'+name,metadata.numeric_metadata(body,units))
        zero=(0,0,0);mom=bound['momentumUnit']
        for item in result['records']:
            if item['method']!='adaptive':continue
            g=item['group'];suffix=str(item['test']);unit=lambda p:bound['rows'][g['rowIndices'][p[0]]]['unit']
            engine.physical(PREFIX+'_OUTER_RULE_'+suffix,'adaptive GK21; disjoint halves; epsrel=0; inner/source/profile rules held')
            numeric('ADAPTIVE_ERROR_TOLERANCE_'+suffix,(g['unitFrameErrorEstimate'],g['unitFrameTolerance']),lambda p:zero)
            numeric('ADAPTIVE_EVALUATIONS_SUCCESS_'+suffix,(g['evaluations'],int(g['success'])),lambda p:zero)
            numeric('ADAPTIVE_MASS_'+suffix,(g['quadratureMass'],g['boxVolume'],g['volumeResidual']),lambda p:tuple(len(variables)*v for v in mom))
            numeric('ADAPTIVE_MEASURE_MUTATED_'+suffix,g['measureMutationValues'],unit)
            numeric('ADAPTIVE_MEASURE_RESIDUAL_'+suffix,g['measureMutationResidual'],unit)
            for hi,h in enumerate(g['halves']):
                key=suffix+'_'+str(hi);n=len(rows)*len(bound['positions'])
                numeric('HALF_INTERVALS_'+key,h['intervals'],lambda p:mom)
                numeric('HALF_INTERVAL_ERRORS_'+key,h['unitFrameIntervalErrors'],lambda p:zero)
                numeric('HALF_ERROR_STATUS_'+key,(h['unitFrameErrorEstimate'],h['unitFrameTolerance'],h['status'],int(h['success']),h['calls'],h['cacheHits'],h['uniquePoints']),lambda p:zero)
                numeric('HALF_INTERVAL_VALUES_'+key,h['intervalPackedValues'],lambda p:tuple(len(variables)*v for v in mom) if p[1]==2*n else bound['rows'][g['rowIndices'][(p[1]%n)//len(bound['positions'])]]['unit'])
    previous_emit,previous_prefix=native.emit,native.PREFIX;native.emit,native.PREFIX=emit,PREFIX
    engine.PAYLOAD_ENCODER=engine.PayloadEncoder();engine.EMISSION_LINES.clear()
    try:return native.emit_and_replay(base,result,r,bound,rows,finest,provenance)
    finally:native.emit,native.PREFIX=previous_emit,previous_prefix


def finish(base,data,packets):
    result=combine(base,data,packets);before={p.name:digest(p) for p in base.glob('*.pickle')}
    entries,keys,paths=emit_and_replay(base,result,data)
    r,bound,rows,variables,finest,provenance,pins=data
    worker_manifest=json.loads((base/'workers.json').read_text())
    worker_artifacts=[]
    for task in sorted(packets):
        path=base/f'worker-{task[0]}-{task[1]}'/'result.pickle'
        require(digest(path)==worker_manifest['resultFiles'][f'{task[0]}-{task[1]}'],'final worker packet hash')
        worker_artifacts.append({'path':str(path.relative_to(base)),'bytes':path.stat().st_size,'sha256':digest(path)})
    summary={'runDirectory':str(base),'sourceFiles':pins,'provenance':provenance,'methodJoins':METHOD_JOINS,'rowCount':len(rows),'records':len(result['records']),
        'norms':[{'test':v['test'],'integralDifference':float(np.max(abs(v['referenceIntegralDifference']))),'actionDifference':float(np.max(abs(v['differenceFromReference']))),'unitFrameErrorEstimate':v['group']['unitFrameErrorEstimate'],'evaluations':v['group']['evaluations']} for v in result['records'] if v['method']=='adaptive'],
        'tagCount':len(entries),'writeKeyCount':len(keys),'metadataPaths':paths,'recordArtifacts':result['recordArtifacts'],'pointArtifacts':result['pointArtifacts'],'partialArtifacts':result['partialArtifacts'],
        'workerManifest':worker_manifest,'workerArtifacts':worker_artifacts,'workerPeakRssKiB':{str(k):v['peakRssKiB'] for k,v in packets.items()},
        'packetHashesBeforeEmission':before,'packetHashesAfterEmission':{n:digest(base/n) for n in before},'artifacts':{p.name:{'bytes':p.stat().st_size,'sha256':digest(p)} for p in base.iterdir() if p.suffix in ('.pickle','.out')},
        'scope':'Independent adaptive outer quadrature only; accepted inner/source/profile rules held. No uniform or independent-grade convergence, physical tails/interchange, Abel limits, scattering or poles.'}
    for key in ('recordArtifacts','pointArtifacts','partialArtifacts','workerArtifacts'):
        for a in summary[key]:require(digest(base/a['path'])==a['sha256'],'final saved hash')
    require(before==summary['packetHashesAfterEmission'] and pins=={n:digest(ROOT/n) for n in pins} and not engine.PHYSICAL_METADATA.dimensions.constraints,'final source/packet/dimension guard')
    save(base/'checks.json',summary)
    return summary


def main():
    p=argparse.ArgumentParser();p.add_argument('--run-directory',type=Path,required=True);a=p.parse_args()
    base=a.run_directory.resolve();base.relative_to(STORE);base.mkdir(parents=True,exist_ok=False)
    started=time.monotonic();data=load(base);packets=dispatch(base,data);summary=finish(base,data,packets)
    summary['wallSeconds']=time.monotonic()-started;summary['coordinatorPeakRssKiB']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    save(base/'checks.json',summary);print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
