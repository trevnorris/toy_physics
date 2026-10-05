#!/usr/bin/env python3
"""Refine the retained native triple layout on both accepted wider boxes."""
import argparse,ast,copy,hashlib,inspect,json,shutil,textwrap,time,types
from pathlib import Path
import numpy as np
import S11c_d_wide_momentum_refinement as wide

ROOT,STORE,M=wide.ROOT,wide.STORE,wide.M
CHECKPOINT=M/'S11c_d_wide_momentum_adaptive_checkpoint.json'
PLAN=M/'S11c_d_wide_three_momentum_plan.md'
PREFIX='WIDE_THREE_MOMENTUM_LAB_HELD_RHO4_CONSTANT'
engine,position,parallel,domain=wide.engine,wide.position,wide.parallel,wide.domain
require,digest,save,atomic_pickle,unpickle,artifact=wide.require,wide.digest,wide.save,wide.atomic_pickle,wide.unpickle,wide.artifact
TASKS=wide.TASKS


def settings(reference,stage,smoke=False):
    require(stage in (0,1,2,3),'triple refinement stage')
    result=copy.deepcopy(reference['layoutSettings'])
    if stage:
        result[3]['outerOrder']=2 if smoke else 216
        result[3]['innerOrders']=((1 if stage<2 else 2),(1 if stage<3 else 2)) if smoke else ((24 if stage<2 else 32),(24 if stage<3 else 32))
        if smoke:result[3].update(sourceNodes=16,profileNodes=16)
    return result


def derived(function,replacements,namespace):
    original=ast.parse(textwrap.dedent(inspect.getsource(function)));tree=copy.deepcopy(original);pairs=[]
    for old,new,mode in replacements:
        a=ast.parse(old,mode=mode);b=ast.parse(new,mode=mode)
        pairs.append((a.body if mode=='eval' else a.body[0],b.body if mode=='eval' else b.body[0]))
    class Replace(ast.NodeTransformer):
        def __init__(self,pairs):self.pairs=pairs;self.counts=[0]*len(pairs)
        def generic_visit(self,node):
            for i,(old,new) in enumerate(self.pairs):
                if ast.dump(node)==ast.dump(old):self.counts[i]+=1;return copy.deepcopy(new)
            return super().generic_visit(node)
    change=Replace(pairs);tree=change.visit(tree);require(change.counts==[1]*len(pairs),'exact worker adaptation census')
    reverse=Replace([(b,a) for a,b in pairs]);restored=reverse.visit(copy.deepcopy(tree))
    require(reverse.counts==[1]*len(pairs) and ast.dump(restored)==ast.dump(original),'whole function reverse AST join')
    ns=dict(vars(wide),**namespace);exec(compile(ast.fix_missing_locations(tree),'<wide-triple-native-adapter>','exec'),ns)
    return ns[function.__name__],{'nativeAstSha256':hashlib.sha256(ast.dump(original).encode()).hexdigest(),'restoredAstSha256':hashlib.sha256(ast.dump(restored).encode()).hexdigest(),'adaptationCount':len(pairs)}


make_record,MAKE_JOIN=derived(wide.make_record,[('n = 0 if stage == 0 else 1 if stage <= 2 else 2','n = 0 if stage == 0 else 3','exec')],{'settings':settings})
guard_record,GUARD_JOIN=derived(wide.guard_record,[("rows = [v for v in bound['rows'] if len(v['limits']) <= max(n,1)]","rows = [v for v in bound['rows'] if len(v['limits']) == n]",'exec')],{})
evaluate_task,WORKER_JOIN=derived(wide.evaluate_task,[('range(7)','range(4)','eval'),('n = 1 if stage <= 2 else 2','n = 3','exec')],{'settings':settings,'make_record':make_record,'guard_record':guard_record})


def load(base):
    accepted,previous=domain.momentum.source.accepted(CHECKPOINT)
    require(accepted['status']=='PUBLISHED_ANNEX_VERIFIED','accepted independent outer publication')
    for n,h in accepted['sourceFiles'].items():require(digest(ROOT/n)==h==digest(previous/'source'/n),'accepted current/frozen source')
    for category in ('recordArtifacts','workerArtifacts'):
        for a in accepted[category]:require(digest(previous/a['path'])==a['sha256'],'accepted saved artifact')
    bp=unpickle(previous/'bound-wide-adaptive.pickle');packet=unpickle(previous/'wide-adaptive.pickle')
    require(bp['sourceFiles']==packet['sourceFiles']==accepted['sourceFiles'] and bp['provenance']==packet['provenance']==accepted['provenance'],'accepted packet provenance')
    _,sb=domain.momentum.source.accepted(domain.momentum.source.SOURCE)
    r,dimensions=domain.momentum.source.native.source.restore_context(unpickle(sb/'reduced-action.pickle'));dimensions.__dict__.update(packet['dimensionState'])
    refs=bp['refs'];domain_checkpoint=M/'S11c_d_momentum_domain_action_checkpoint.json'
    require(digest(domain_checkpoint)==accepted['sourceFiles'][str(domain_checkpoint.relative_to(ROOT))],'pinned complete-domain checkpoint')
    original=Path(json.loads(domain_checkpoint.read_text())['runDirectory']);old=unpickle(original/'momentum-domain-action.pickle')
    original_cp=json.loads(domain_checkpoint.read_text())
    require(digest(original/'momentum-domain-action.pickle')==original_cp['artifacts']['momentum-domain-action.pickle']['sha256'],'original domain result hash')
    for task in TASKS:
        baseline=next(v for v in packet['result']['records'] if (v['test'],v['domain'],v['stage'])==(*task,0));ref=refs[task]
        rebuilt=position.assemble(r,bp['domains'][task[1]],ref,ref['groups'])
        require(np.array_equal(ref['action'],baseline['action']) and np.array_equal(rebuilt['action'],baseline['action']) and all(not np.any(v['differences']) for v in rebuilt['terms']),'accepted Gauss full-action reference')
        require(ref['layoutSettings']==baseline['settings'],'accepted Gauss rule identity')
        now=next(g for g in baseline['groups'] if len(g['variables'])==3);oldrec=next(v for v in old['result']['records'] if (v['test'],v['index'])==task);prior=next(g for g in oldrec['groups'] if len(g['variables'])==3)
        require(all(np.array_equal(now[k],prior[k]) for k in ('values','measureMutationValues','measureMutationResidual')) and now['rowIndices']==prior['rowIndices'] and now['nodeCount']==prior['nodeCount'],'original full triple baseline join')
        require(ref['layoutSettings'][3]==oldrec['settings'][3],'original triple rules');s=ref['layoutSettings'][3]
        require((s['outerOrder'],tuple(s['innerOrders']),s['sourceNodes'],s['profileNodes'],s['sourceBound'],s['profileBound'],s['regulator'])==(144,(24,24),256,512,48.,14.,.2),'retained native triple setting')
    pins=dict(accepted['sourceFiles'])
    for path in (CHECKPOINT,PLAN,Path(__file__).resolve()):pins[str(path.relative_to(ROOT))]=digest(path)
    for n in pins:
        dest=base/'source'/n;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/n,dest)
    for name in ('bound-wide-adaptive.pickle','wide-adaptive.pickle'):
        shutil.copyfile(previous/name,base/('accepted-'+name));require(digest(base/('accepted-'+name))==accepted['artifacts'][name]['sha256'],'accepted packet byte copy')
    provenance=dict(accepted['provenance'],WIDE_THREE_INPUT_SHA256=digest(CHECKPOINT),WIDE_THREE_INSTRUMENT_SHA256=digest(Path(__file__)),WIDE_THREE_PLAN_SHA256=digest(PLAN))
    data={'r':r,'domains':bp['domains'],'joins':bp['joins'],'refs':refs,'sourceFiles':pins,'provenance':provenance,'previous':previous,'originalDomainRun':original,'smoke':False}
    atomic_pickle(base/'bound-wide-three.pickle',{k:v for k,v in data.items() if k!='r'}|{'dimensionState':dict(vars(dimensions))})
    save(base/'preflight.json',{'sourceFiles':pins,'provenance':provenance,'tasks':TASKS,'acceptedRunDirectory':str(previous),'originalDomainRun':str(original),'methodJoins':METHOD_JOINS,'nativeEngineUnchanged':True})
    return data


def dispatch(base,data):
    old=parallel.evaluate_task;parallel.evaluate_task=evaluate_task
    try:return parallel.dispatch(base,data,TASKS)
    finally:parallel.evaluate_task=old


emit=types.FunctionType(wide.emit.__code__,dict(vars(wide),PREFIX=PREFIX),wide.emit.__name__)
replay_adapter=types.FunctionType(wide.replay_adapter.__code__,dict(vars(wide),PREFIX=PREFIX,emit=emit),wide.replay_adapter.__name__)
SCOPE='Fixed wide boxes: triple quadrature changes with independently checked single/pair Gauss values held. No independent wide-triple adaptive result, uniform or independent-grade/global exceptional coverage, physical tails, Abel limits, scattering or poles.'
old_scope='Fixed wide boxes: single and paired quadrature changes with accepted three-momentum values explicitly held. No independent adaptive check, full wide-box convergence, physical tails, Abel limit, scattering or poles.'
finish,FINISH_JOIN=derived(wide.finish,[(repr('wide-refinement.pickle'),repr('wide-three.pickle'),'eval'),(repr(old_scope),repr(SCOPE),'eval')],{'guard_record':guard_record,'replay_adapter':replay_adapter})
METHOD_JOINS={'record':MAKE_JOIN,'guard':GUARD_JOIN,'worker':WORKER_JOIN,'finish':FINISH_JOIN}


def focused_checks(base,data):
    evidence=[];timings=[];isolations=[]
    class PrefixComplete(Exception):pass
    for task in TASKS:
        ti,di=task;ref=data['refs'][task];bound=data['domains'][di];r=data['r'];variables=next(g['variables'] for g in ref['groups'] if len(g['variables'])==3)
        for stage in (0,3):
            s=settings(ref,stage)[3];worker=engine.BoundedSourceFourierQuadrature.ThreeMomentum(bound['rows'],bound['sources'],r);captured=[];started=time.monotonic()
            def stop(state):
                if state['batchCount']==64:captured.append(state);raise PrefixComplete()
            try:worker.group(ti,variables,s,bound['pairs'],float(bound['abel']['width'].subs(r.regulator,s['regulator'])),bound['positions'],stop)
            except PrefixComplete:pass
            require(len(captured)==1,'complete bounded prefix');actual=captured[0];seconds=time.monotonic()-started
            proof={'task':task,'stage':stage,'actual':actual,'wallSeconds':seconds,'telemetry':parallel.cache_state(worker),'units':[v['unit'] for v in bound['rows'] if len(v['limits'])==3],'massUnit':tuple(3*v for v in bound['momentumUnit'])}
            if stage==0:
                expected_path=data['originalDomainRun']/f'worker-{ti}-{di}'/'partials'/'3-0000064.pickle'
                original_cp=json.loads((M/'S11c_d_momentum_domain_action_checkpoint.json').read_text());rel=str(expected_path.relative_to(data['originalDomainRun']))
                require(digest(expected_path)==next(v['sha256'] for v in original_cp['workerArtifacts'] if v['path']==rel),'original actual partial hash')
                expected=unpickle(expected_path);proof.update(expected=expected,expectedPath=str(expected_path),expectedSha256=digest(expected_path),valueResidual=actual['values']-expected['values'],mutationResidual=actual['measureMutationValues']-expected['measureMutationValues'],massResidual=actual['quadratureMass']-expected['quadratureMass'])
            path=base/f'prefix-{ti}-{di}-{stage}.pickle';atomic_pickle(path,proof);evidence.append(artifact(base,path));timings.append({'task':task,'stage':stage,'nodes':actual['nodeCount'],'wallSeconds':seconds})
            require(actual['nodeCount']==16384 and actual['rowIndices']==[v['index'] for v in bound['rows'] if len(v['limits'])==3] and actual['variables']==variables and actual['test']==ti and actual['setting']==s,'prefix row/node/field/limit join')
            require(np.isfinite(actual['values']).all() and np.isfinite(actual['measureMutationValues']).all() and not np.any((abs(actual['values'])>1e-9)&(abs(actual['measureMutationValues']-actual['values'])<1e-12)),'prefix actual measure response')
            if stage==0:
                require(not np.any(proof['valueResidual']) and not np.any(proof['mutationResidual']) and proof['massResidual']==0 and actual['nodeCount']==expected['nodeCount'] and actual['setting']==expected['setting'] and actual['rowIndices']==expected['rowIndices'] and np.array_equal(actual['positions'],expected['positions']),'exact saved native prefix')
        for stage in (1,2,3):
            a,b=settings(ref,stage-1),settings(ref,stage);changes={(n,k) for n in a for k in set(a[n])|set(b[n]) if a[n].get(k)!=b[n].get(k)}
            require(changes=={(3,'outerOrder' if stage==1 else 'innerOrders')},'isolated triple setting')
            if stage>1:require([i for i,(x,y) in enumerate(zip(a[3]['innerOrders'],b[3]['innerOrders'])) if x!=y]==[stage-2],'one inner coordinate changed')
            isolations.append({'task':task,'stage':stage,'before':a,'after':b})
    atomic_pickle(base/'settings-isolation.pickle',isolations)
    return {'prefixArtifacts':evidence,'prefixTimings':timings,'settingsArtifact':artifact(base,base/'settings-isolation.pickle'),'productionSettings':len(isolations),'methodJoins':METHOD_JOINS}


coarse_checks=wide.coarse_checks


def main():
    p=argparse.ArgumentParser();p.add_argument('--run-directory',type=Path,required=True);p.add_argument('--preflight',action='store_true');p.add_argument('--focused-only',action='store_true');a=p.parse_args()
    base=a.run_directory.resolve();base.relative_to(STORE);base.mkdir(parents=True,exist_ok=False);started=time.monotonic();data=load(base)
    focused=focused_checks(base,data) if a.preflight or a.focused_only else None
    if focused is not None:save(base/'focused-before-workers.json',focused)
    if a.focused_only:print(json.dumps(focused,indent=2));return
    data['smoke']=a.preflight;workers=dispatch(base,data)
    if a.preflight:focused['coarseArtifacts']=coarse_checks(base,data,workers);save(base/'focused-after-workers.json',focused)
    summary=finish(base,data,workers,started,focused);summary['methodJoins']=METHOD_JOINS;save(base/'checks.json',summary);print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
