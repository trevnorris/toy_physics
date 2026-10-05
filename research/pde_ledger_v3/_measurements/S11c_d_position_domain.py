#!/usr/bin/env python3
"""Finite position-domain action changes from accepted native operators."""
import argparse,ast,contextlib,copy,hashlib,inspect,json,resource,shutil,textwrap,time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import sympy as sp
import S11c_d_three_momentum_adaptive as adaptive

parallel=adaptive.parallel;native=adaptive.native;momentum=native.momentum
engine=adaptive.engine
ROOT,STORE=adaptive.ROOT,adaptive.STORE
require,digest,save,atomic_pickle,unpickle=adaptive.require,adaptive.digest,adaptive.save,adaptive.atomic_pickle,adaptive.unpickle
M=ROOT/'_measurements';CHECKPOINT=M/'S11c_d_three_momentum_adaptive_checkpoint.json'
PLAN=M/'S11c_d_position_domain_plan.md';PREFIX='POSITION_DOMAIN_LAB_HELD_RHO4_CONSTANT'
CHOICES=((32.,10.),(48.,10.),(48.,14.))


def artifact(base,path):return {'path':str(path.relative_to(base)),'bytes':path.stat().st_size,'sha256':digest(path)}


def rebind(r,bound,choice):
    source_cut,profile_cut=map(lambda x:sp.Rational(str(x)),choice)
    replacements={};units={};proofs=[]
    for old,unit in bound['profileUnits'].items():
        require(old.limits==((r.xi,-10,10),),'accepted profile limits')
        fresh=sp.Integral(old.function,(r.xi,-profile_cut,profile_cut))
        replacements[old]=fresh;units[fresh]=unit
        proofs.append({'old':old,'new':fresh,'integrandResidual':fresh.function-old.function,
            'restored':sp.Integral(fresh.function,*old.limits),'unit':unit})
    require(len(replacements)==len(units)==6,'complete profile map')
    reverse={v:k for k,v in replacements.items()};rows=[];row_proofs=[]
    for old in bound['rows']:
        require(tuple(old['sourceLimit'])==(r.zp,-32,32),'accepted source limit')
        fresh=dict(old,sourceLimit=sp.Tuple(r.zp,-source_cut,source_cut),factors=[])
        for f in old['factors']:
            coefficient=engine.memo_xreplace(f['coefficient'],replacements)
            restored=engine.memo_xreplace(coefficient,reverse)
            require(restored==f['coefficient'],'coefficient reverse join')
            fresh['factors'].append(dict(f,coefficient=coefficient))
        require(fresh['limits']==old['limits'] and fresh['original']==old['original'],'native remaining limits/original')
        rows.append(fresh);row_proofs.append({'row':old['index'],'oldSourceLimit':old['sourceLimit'],
            'newSourceLimit':fresh['sourceLimit'],'remainingLimits':fresh['limits'],
            'coefficientReverseResiduals':[engine.memo_xreplace(f['coefficient'],reverse)-g['coefficient'] for f,g in zip(fresh['factors'],old['factors'])]})
    source_symbols={v['originalSourceIntegral'].limits[0][2] for v in bound['sources'].values()}
    profile_symbols={i.limits[0][2] for row in bound['rows'] for f in row['factors'] for i in f['symbolicCoefficient'].atoms(sp.Integral) if i.limits[0][0]==r.xi}
    require(len(source_symbols)==len(profile_symbols)==1 and source_symbols.isdisjoint(profile_symbols),'symbolic cutoff census')
    cutoffs=dict(bound['cutoffBindings'])
    for symbols,value,old_value in ((source_symbols,source_cut,32),(profile_symbols,profile_cut,10)):
        for symbol in symbols:require(cutoffs[symbol]==old_value,'accepted symbolic cutoff binding');cutoffs[symbol]=value
    actual=set().union(*(f['coefficient'].atoms(sp.Integral) for row in rows for f in row['factors']))
    require(actual==set(units),'rebound profile census')
    result=dict(bound,rows=rows,profileUnits=units,cutoffBindings=cutoffs)
    return result,{'choice':choice,'profiles':proofs,'rows':row_proofs,'cutoffBindings':cutoffs}


def assemble(r,bound,reference,groups):
    row_indices=[i for g in groups for i in g['rowIndices']]
    require(len(row_indices)==len(set(row_indices))==80 and set(row_indices)==set(range(80)),'all native rows')
    merged={'rowIndices':row_indices,'values':np.concatenate([g['values'] for g in groups],axis=0)}
    # Re-enter the already validated native full-action contraction, replacing all
    # native rows at once. No unchanged nonlocal contribution can be substituted.
    result=native.full_action(r,bound,reference,merged)
    require(all(all(v['changed']) for v in result['terms']),'no held nonlocal term')
    return result


def settings(reference,choice,smoke=False):
    result={}
    for n,value in reference['layoutSettings'].items():
        s=dict(value,sourceBound=choice[0],profileBound=choice[1],sourceNodes=256,profileNodes=256)
        s['innerOrders']=tuple(value.get('innerOrders',(value['panelOrder'],)*(n-1)))
        if smoke:s.update(outerOrder=2,panelOrder=1,innerOrders=(1,)*(n-1),sourceNodes=16,profileNodes=16)
        result[n]=s
    return result


def load(base):
    accepted,previous=momentum.source.accepted(CHECKPOINT)
    for n,h in accepted['sourceFiles'].items():require(digest(ROOT/n)==h==digest(previous/'source'/n),('accepted source',n))
    for key in ('recordArtifacts','pointArtifacts','partialArtifacts','workerArtifacts'):
        for a in accepted[key]:require(digest(previous/a['path'])==a['sha256'],'accepted adaptive artifact')
    r,bound,_,_,finest,provenance,pins=adaptive.load(base)
    shutil.copyfile(previous/'three-momentum-adaptive.pickle',base/'accepted-three-momentum-adaptive.pickle')
    packet=unpickle(base/'accepted-three-momentum-adaptive.pickle')
    for ti,ref in finest.items():
        held=next(v for v in packet['result']['records'] if v['test']==ti and v['method']=='retained')
        require(np.array_equal(held['action'],ref['action']),'accepted adaptive/Gauss full-action join')
        rebuilt=assemble(r,bound,ref,ref['groups']);require(np.array_equal(rebuilt['action'],ref['action']),'all-layout action reconstruction')
        for v in rebuilt['terms']:require(not np.any(v['differences']),'baseline native contribution reconstruction')
    domains={};joins={}
    for index,choice in enumerate(CHOICES):domains[index],joins[index]=rebind(r,bound,choice)
    require(domains[0]['rows']==bound['rows'] and domains[0]['profileUnits']==bound['profileUnits'],'exact baseline rebound')
    paths=[*(ROOT/n for n in accepted['sourceFiles']),CHECKPOINT,PLAN,Path(__file__).resolve(),M/'S11c_d_position_domain_preflight.py']
    for path in paths:
        n=str(path.relative_to(ROOT));pins[n]=digest(path);target=base/'source'/n;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,target)
    provenance=dict(provenance,ADAPTIVE_CHECKPOINT_SHA256=digest(CHECKPOINT),POSITION_DOMAIN_INSTRUMENT_SHA256=digest(Path(__file__)),POSITION_DOMAIN_PLAN_SHA256=digest(PLAN))
    atomic_pickle(base/'bound-position-domains.pickle',{'domains':domains,'joins':joins,'references':finest,'provenance':provenance,'sourceFiles':pins,'dimensionState':dict(vars(engine.PHYSICAL_METADATA.dimensions))})
    data={'r':r,'bound':bound,'domains':domains,'joins':joins,'finest':finest,'provenance':provenance,'sourceFiles':pins,'smoke':False}
    save(base/'preflight.json',{'sourceFiles':pins,'provenance':provenance,'choices':CHOICES,'rows':len(bound['rows']),'sources':len(bound['sources']),'profiles':len(bound['profileUnits']),'baselineActionReconstruction':True,'acceptedRunDirectory':str(previous)})
    return data


def evaluate_task(directory,data,task,*,max_batches=None,resume_state=None):
    require(max_batches is None and resume_state is None,'position-domain task request')
    ti,index=task;require(ti in (0,1) and index in (1,2),'position-domain task')
    r=data['r'];bound=data['domains'][index];reference=data['finest'][ti]
    rules=settings(reference,CHOICES[index],data['smoke'])
    worker=engine.BoundedSourceFourierQuadrature.ThreeMomentum(bound['rows'],bound['sources'],r)
    groups=[];artifacts=[];started=time.monotonic()
    for original in reference['groups']:
        variables=tuple(original['variables']);n=len(variables);setting=rules[n]
        width=float(bound['abel']['width'].subs(r.regulator,setting['regulator']))
        def partial(state):
            if state['batchCount']%64==0:
                path=directory/'partials'/f'{n}-{state["batchCount"]:07}.pickle';path.parent.mkdir(exist_ok=True)
                atomic_pickle(path,state);artifacts.append(artifact(directory,path));save(directory/'partial-inventory.json',artifacts)
        group=worker.group(ti,variables,setting,bound['pairs'],width,bound['positions'],partial)
        path=directory/f'layout-{n}.pickle';atomic_pickle(path,{'setting':setting,'group':group});artifacts.append(artifact(directory,path));save(directory/'partial-inventory.json',artifacts)
        require(group['rowIndices']==original['rowIndices'],'all native rows/limit order')
        require(abs(group['volumeResidual'])<1e-10*(1+abs(group['boxVolume'])),'finite mass')
        require(np.all(np.isfinite(group['values'])) and np.all(np.isfinite(group['measureMutationValues'])),'finite integrals')
        require(not np.any((abs(group['values'])>1e-9)&(abs(group['measureMutationResidual'])<1e-12)),'measure sensitivity')
        groups.append(group)
    telemetry=parallel.cache_state(worker)
    require(set(telemetry['sourceFrequencyCensus'])=={(ti,si) for t,si in bound['sources'] if t==ti},'all source census')
    require(telemetry['evaluatedProfileIntegrals']==set(bound['profileUnits']),'all profile census')
    result={'task':task,'groups':groups,'settings':rules,**assemble(r,bound,reference,groups),'telemetry':telemetry,
        'partialArtifacts':artifacts,'sourceFiles':data['sourceFiles'],'provenance':data['provenance'],
        'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    atomic_pickle(directory/'result.pickle',result)
    require(data['sourceFiles']=={n:digest(ROOT/n) for n in data['sourceFiles']},'worker source change')
    save(directory/'checks.json',{'task':task,'resultSha256':digest(directory/'result.pickle'),'wallSeconds':result['wallSeconds'],'peakRssKiB':result['peakRssKiB']})
    return result


def dispatch(base,data):
    old=parallel.evaluate_task;parallel.evaluate_task=evaluate_task
    try:return parallel.dispatch(base,data,[(ti,i) for ti in range(2) for i in (1,2)])
    finally:parallel.evaluate_task=old


def combine(base,data,workers):
    records=[];inventory=[];artifacts=[]
    for ti in range(2):
        reference=data['finest'][ti]
        baseline={'test':ti,'index':0,'method':'retained','groups':reference['groups'],
            'settings':{n:dict(s,innerOrders=tuple(s.get('innerOrders',(s['panelOrder'],)*(n-1)))) for n,s in reference['layoutSettings'].items()},
            **assemble(data['r'],data['bound'],reference,reference['groups'])}
        previous=baseline
        for index in range(3):
            if index==0:item=baseline
            else:
                w=workers[(ti,index)];directory=base/f'worker-{ti}-{index}'
                require(w['task']==(ti,index) and w['settings']==settings(reference,CHOICES[index],data['smoke']),'worker assignment/settings')
                require(w['sourceFiles']==data['sourceFiles'] and w['provenance']==data['provenance'],'worker source joins')
                artifacts.append(artifact(base,directory/'result.pickle'))
                for a in w['partialArtifacts']:
                    require(digest(directory/a['path'])==a['sha256'],'worker partial hash')
                    artifacts.append(dict(a,path=str((directory/a['path']).relative_to(base))))
                item=dict(w,test=ti,index=index,method='gauss',precedingActionDifference=w['action']-previous['action'])
                item['integralChanges']=[b['values']-a['values'] for a,b in zip(previous['groups'],w['groups'])]
                item['termChanges']=[np.asarray(b['values'])-np.asarray(a['values']) for a,b in zip(previous['terms'],w['terms'])]
                for a,b in zip(previous['terms'],w['terms']):require(a['index']==b['index'] and a['rows']==b['rows'],'native term joins')
            path=base/'records'/f'{ti}-{index}.pickle';path.parent.mkdir(exist_ok=True);atomic_pickle(path,item);inventory.append(artifact(base,path));save(base/'record-inventory.json',inventory)
            records.append(item);previous=item
    result={'records':records,'recordArtifacts':inventory,'workerArtifacts':artifacts,'joins':data['joins']}
    atomic_pickle(base/'position-domain.pickle',{'result':result,'provenance':data['provenance'],'sourceFiles':data['sourceFiles'],'dimensionState':dict(vars(engine.PHYSICAL_METADATA.dimensions))})
    return result


def emit(result,r,bound,rows,finest,provenance):
    dimensions=engine.PHYSICAL_METADATA.dimensions;zero=dimensions.zero
    metadata=engine.FullPencilModes.__new__(engine.FullPencilModes);metadata.r=r
    def numeric(name,value,unit,literal=False):
        body=momentum.number(value);engine.emit(PREFIX+'_'+name,body if literal else metadata.compact_fingerprint(body))
        engine.emit('METADATA_'+PREFIX+'_'+name,metadata.numeric_metadata(body,unit))
    engine.physical(PREFIX+'_PROVENANCE',provenance)
    engine.physical(PREFIX+'_ABEL_PAIRS',bound['pairs']);engine.physical(PREFIX+'_ABEL_WIDTH',bound['abel']['width'])
    for row in bound['rows']:
        engine.fingerprinted(PREFIX+'_ORIGINAL_'+str(row['index']),row['original'])
        numeric('REMAINING_LIMITS_'+str(row['index']),row['limits'],lambda p:bound['momentumUnit'],True)
    for index,join in result['joins'].items():
        numeric('SOURCE_BOUND_'+str(index),join['choice'][0],lambda p:dimensions.measure(r.zp),True)
        numeric('PROFILE_BOUND_'+str(index),join['choice'][1],lambda p:dimensions.measure(r.xi),True)
        for j,profile in enumerate(join['profiles']):
            for key in ('old','new','restored','integrandResidual'):
                value=profile[key];name=f'{index}_{j}_{key.upper()}'
                engine.emit(PREFIX+'_PROFILE_'+name,engine.carrier_fingerprint(value))
                engine.emit('METADATA_'+PREFIX+'_PROFILE_'+name,metadata.numeric_metadata(value,lambda p:profile['unit']))
        for row in join['rows']:
            numeric('REBOUND_SOURCE_LIMIT_'+str(index)+'_'+str(row['row']),row['newSourceLimit'],lambda p:dimensions.measure(r.zp),True)
            numeric('COEFFICIENT_REVERSE_RESIDUAL_'+str(index)+'_'+str(row['row']),row['coefficientReverseResiduals'],lambda p:bound['rows'][row['row']]['factors'][p[0]]['unit'],True)
    for item in result['records']:
        key=f'{item["test"]}_{item["index"]}'
        engine.physical(PREFIX+'_METHOD_'+key,item['method'])
        for name in ('action','differenceFromReference','precedingActionDifference'):
            if name in item:numeric(name.upper()+'_'+key,item[name],lambda p:bound['equationUnits'][p[2]],name!='action')
        for j,group in enumerate(item['groups']):
            n=len(group['variables']);suffix=key+'_'+str(j);setting=item['settings'][n]
            numeric('ROWS_'+suffix,group['rowIndices'],lambda p:zero,True)
            engine.physical(PREFIX+'_VARIABLES_'+suffix,group['variables'])
            numeric('LAYOUT_ORDERS_'+suffix,(setting['outerOrder'],*setting['innerOrders'],setting['sourceNodes'],setting['profileNodes']),lambda p:zero,True)
            for name,unit in (('sourceBound',dimensions.measure(r.zp)),('profileBound',dimensions.measure(r.xi)),('momentumBound',bound['momentumUnit']),('regulator',dimensions.measure(r.regulator))):numeric(name.upper()+'_'+suffix,setting[name],lambda p:unit,True)
            unit=lambda p:bound['rows'][group['rowIndices'][p[0]]]['unit']
            for name in ('values','measureMutationValues','measureMutationResidual'):numeric(name.upper()+'_'+suffix,group[name],unit,name.endswith('Residual'))
            if 'integralChanges' in item:numeric('INTEGRAL_CHANGE_'+suffix,item['integralChanges'][j],unit,True)
            numeric('MASS_'+suffix,(group['quadratureMass'],group['boxVolume'],group['volumeResidual']),lambda p:tuple(n*v for v in bound['momentumUnit']),True)
        for j,term in enumerate(item['terms']):
            suffix=key+'_'+str(j);unit=lambda p:bound['equationUnits'][term['index'][2]]
            numeric('CELL_'+suffix,term['index'],lambda p:zero,True);numeric('TERM_ROWS_'+suffix,term['rows'],lambda p:zero,True)
            for name in ('reference','values','differences'):numeric('TERM_'+name.upper()+'_'+suffix,term[name],unit,name=='differences')
            if 'termChanges' in item:numeric('TERM_CHANGE_'+suffix,item['termChanges'][j],unit,True)
        if 'telemetry' in item:
            t=item['telemetry']
            for j,c in enumerate(t['cacheChecks']):
                suffix=key+'_'+str(j);u=dimensions.measure(r.xi if c['profileRule'] else r.zp)
                numeric('CACHE_KIND_'+suffix,int(c['profileRule']),lambda p:zero,True)
                numeric('CACHE_PANELS_'+suffix,c['points'],lambda p:u,True)
                numeric('CACHE_ORDER_'+suffix,c['order'],lambda p:zero,True)
                for name in ('nodesResidual','weightsResidual'):numeric('CACHE_'+name.upper()+'_'+suffix,c[name],lambda p:u,True)
                numeric('CACHE_READONLY_'+suffix,int(c['readOnly']),lambda p:zero,True)
            for (ti,si),v in sorted(t['sourceFrequencyCensus'].items()):numeric('FREQUENCIES_'+key+'_'+str(si),v,lambda p:zero if p[0]==2 else bound['momentumUnit'],True)
            numeric('WORKSPACE_'+key,[t[k] for k in ('phaseWorkspaceEstimateBytes','batchCacheEstimateBytes','workspaceBudgetBytes')],lambda p:zero,True)


def cache_for_tag(result,tag):
    prefixes=tuple('PY_S11CD_METADATA_'+PREFIX+'_CACHE_'+name+'_' for name in ('NODESRESIDUAL','WEIGHTSRESIDUAL'))
    prefix=next((p for p in prefixes if tag.startswith(p)),None)
    require(prefix is not None,'position-domain cache metadata tag')
    suffix=tag.removeprefix(prefix).split('_')
    require(len(suffix)==3 and all(v.isdigit() for v in suffix),'position-domain cache metadata coordinates')
    test,index,cache=map(int,suffix)
    records=[v for v in result['records'] if v['test']==test and v['index']==index and 'telemetry' in v]
    require(len(records)==1,'position-domain cache metadata record')
    caches=records[0]['telemetry']['cacheChecks']
    require(0<=cache<len(caches),'position-domain cache metadata index')
    return caches[cache]


def replay_adapter():
    original=ast.parse(textwrap.dedent(inspect.getsource(native.emit_and_replay)))
    tree=copy.deepcopy(original)
    assignments=[n for n in ast.walk(tree) if isinstance(n,ast.Assign) and len(n.targets)==1 and isinstance(n.targets[0],ast.Name) and n.targets[0].id=='cache']
    require(len(assignments)==1 and ast.unparse(assignments[0].value)=="result['cacheChecks'][int(tag.rsplit('_', 1)[1])]",'native cache lookup census')
    previous=copy.deepcopy(assignments[0].value)
    assignments[0].value=ast.parse('cache_for_tag(result,tag)',mode='eval').body
    restored=copy.deepcopy(tree)
    next(n for n in ast.walk(restored) if isinstance(n,ast.Assign) and len(n.targets)==1 and isinstance(n.targets[0],ast.Name) and n.targets[0].id=='cache').value=previous
    require(ast.dump(restored)==ast.dump(original),'full native replayer AST join')
    namespace=dict(vars(native),emit=emit,PREFIX=PREFIX,cache_for_tag=cache_for_tag)
    exec(compile(ast.fix_missing_locations(tree),'<position-domain-cache-replay>','exec'),namespace)
    return namespace['emit_and_replay'],{'nativeReplayAstSha256':hashlib.sha256(ast.dump(original).encode()).hexdigest(),'restoredReplayAstSha256':hashlib.sha256(ast.dump(restored).encode()).hexdigest()}


def finish(base,data,workers,started):
    result=combine(base,data,workers);before={p.name:digest(p) for p in base.glob('*.pickle')}
    replay,method_join=replay_adapter()
    engine.EMISSION_LINES.clear();engine.PAYLOAD_ENCODER=engine.PayloadEncoder()
    entries,keys,paths=replay(base,result,data['r'],data['bound'],data['bound']['rows'],data['finest'],data['provenance'])
    norms=[{'test':v['test'],'index':v['index'],'rawIntegralChanges':[float(np.max(abs(a))) for a in v['integralChanges']],
        'actionChange':float(np.max(abs(v['precedingActionDifference']))),'rawTermChange':max(float(np.max(abs(t))) if len(t) else 0. for t in v['termChanges'])} for v in result['records'] if 'integralChanges' in v]
    summary={'runDirectory':str(base),'sourceFiles':data['sourceFiles'],'provenance':data['provenance'],'smoke':data['smoke'],
        'rows':80,'sources':70,'profiles':6,'records':len(result['records']),'norms':norms,'tagCount':len(entries),'writeKeyCount':len(keys),'metadataPaths':paths,'replayMethodJoin':method_join,
        'recordArtifacts':result['recordArtifacts'],'workerArtifacts':result['workerArtifacts'],'workerManifest':json.loads((base/'workers.json').read_text()),
        'packetHashesBeforeEmission':before,'packetHashesAfterEmission':{n:digest(base/n) for n in before},
        'artifacts':{p.name:{'bytes':p.stat().st_size,'sha256':digest(p)} for p in base.iterdir() if p.suffix in ('.pickle','.out')},
        'wallSeconds':time.monotonic()-started,'coordinatorPeakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'scope':'Finite source/profile domain changes only at held finite momentum domain and positive regulator. No infinite tail bound, uniform/independent-grade coverage, Abel limit, scattering or pole result.'}
    for a in (*result['recordArtifacts'],*result['workerArtifacts']):require(digest(base/a['path'])==a['sha256'],'final artifact hash')
    require(before==summary['packetHashesAfterEmission'] and data['sourceFiles']=={n:digest(ROOT/n) for n in data['sourceFiles']} and not engine.PHYSICAL_METADATA.dimensions.constraints,'final packet/source/dimension guard')
    save(base/'checks.json',summary);return summary


def main():
    p=argparse.ArgumentParser();p.add_argument('--run-directory',type=Path,required=True);a=p.parse_args()
    base=a.run_directory.resolve();base.relative_to(STORE);base.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    data=load(base);workers=dispatch(base,data);summary=finish(base,data,workers,started);print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
