#!/usr/bin/env python3
"""Independent outer quadrature of native single/pair rows on fixed wide boxes."""
import argparse, ast, copy, hashlib, inspect, json, resource, shutil, textwrap, time, types
from pathlib import Path
import numpy as np
from scipy.integrate import quad_vec
import S11c_d_wide_momentum_refinement as wide
import S11c_d_three_momentum_adaptive as conditional

ROOT, STORE, M = wide.ROOT, wide.STORE, wide.M
CHECKPOINT = M/'S11c_d_wide_momentum_refinement_checkpoint.json'
PLAN = M/'S11c_d_wide_momentum_adaptive_plan.md'
PREFIX = 'WIDE_MOMENTUM_ADAPTIVE_LAB_HELD_RHO4_CONSTANT'
engine, position, parallel, domain = wide.engine, wide.position, wide.parallel, wide.domain
require, digest, save, atomic_pickle, unpickle, artifact = wide.require, wide.digest, wide.save, wide.atomic_pickle, wide.unpickle, wide.artifact
TASKS = wide.TASKS


def settings(reference, stage, smoke=False):
    rules = copy.deepcopy(reference['layoutSettings'])
    if smoke and stage:
        for n in ((1,) if stage == 1 else (1, 2)):
            rules[n].update(sourceNodes=16, profileNodes=16)
        if stage == 3: rules[2].update(panelOrder=1, innerOrders=(1,))
    return rules


native_make_record = types.FunctionType(wide.make_record.__code__, dict(vars(wide), settings=settings), 'make_record')


def load(base):
    accepted, previous = domain.momentum.source.accepted(CHECKPOINT)
    require(accepted['status'] == 'PUBLISHED_ANNEX_VERIFIED', 'accepted wide refinement publication')
    for n,h in accepted['sourceFiles'].items():
        require(digest(ROOT/n) == h == digest(previous/'source'/n), 'accepted current/frozen source')
    for category in ('recordArtifacts', 'workerArtifacts'):
        for a in accepted[category]: require(digest(previous/a['path']) == a['sha256'], 'accepted saved artifact')
    bp = unpickle(previous/'bound-wide-refinement.pickle'); packet = unpickle(previous/'wide-refinement.pickle')
    require(bp['sourceFiles'] == packet['sourceFiles'] == accepted['sourceFiles'] and bp['provenance'] == packet['provenance'] == accepted['provenance'], 'accepted packet joins')
    _, source_base = domain.momentum.source.accepted(domain.momentum.source.SOURCE)
    r, dimensions = domain.momentum.source.native.source.restore_context(unpickle(source_base/'reduced-action.pickle'))
    dimensions.__dict__.update(packet['dimensionState']); refs = {}
    for task in TASKS:
        rec = next(v for v in packet['result']['records'] if (v['test'], v['domain'], v['stage']) == (*task, 6))
        refs[task] = dict(bp['refs'][task], action=rec['action'], groups=rec['groups'],
                         contributions=[{'index':v['index'], 'terms':v['values']} for v in rec['terms']],
                         layoutSettings=rec['settings'], setting=rec['settings'][3])
        rebuilt = position.assemble(r, bp['domains'][task[1]], refs[task], rec['groups'])
        require(np.array_equal(rebuilt['action'], rec['action']) and all(not np.any(v['differences']) for v in rebuilt['terms']), 'accepted full action reconstruction')
        s = rec['settings']; require(s[1]['outerOrder'] == s[2]['outerOrder'] == 432 and tuple(s[2]['innerOrders']) == (64,), 'accepted refined rules')
        require(all((v['sourceNodes'],v['profileNodes'],v['sourceBound'],v['profileBound'],v['regulator']) == (256,512,48.,14.,.2) for v in s.values()), 'accepted fixed source/profile/domain/regulator')
    pins = dict(accepted['sourceFiles'])
    for path in (CHECKPOINT, PLAN, Path(__file__).resolve()): pins[str(path.relative_to(ROOT))] = digest(path)
    for n in pins:
        dest = base/'source'/n; dest.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(ROOT/n,dest)
    for name in ('bound-wide-refinement.pickle','wide-refinement.pickle'):
        shutil.copyfile(previous/name,base/('accepted-'+name)); require(digest(base/('accepted-'+name)) == accepted['artifacts'][name]['sha256'], 'accepted packet byte copy')
    provenance = dict(accepted['provenance'], WIDE_ADAPTIVE_INPUT_SHA256=digest(CHECKPOINT), WIDE_ADAPTIVE_INSTRUMENT_SHA256=digest(Path(__file__)), WIDE_ADAPTIVE_PLAN_SHA256=digest(PLAN))
    data = {'r':r,'domains':bp['domains'],'joins':bp['joins'],'refs':refs,'sourceFiles':pins,'provenance':provenance,'previous':previous,'smoke':False}
    atomic_pickle(base/'bound-wide-adaptive.pickle',{k:v for k,v in data.items() if k != 'r'} | {'dimensionState':dict(vars(dimensions))})
    save(base/'preflight.json',{'sourceFiles':pins,'provenance':provenance,'tasks':TASKS,'methodJoins':conditional.METHOD_JOINS,'acceptedRunDirectory':str(previous),'nativeEngineUnchanged':True})
    return data


def adaptive_group(directory, data, task, worker, variables, setting, artifacts):
    ti,di = task; r=data['r']; bound=data['domains'][di]; layout=len(variables)
    rows=[v for v in bound['rows'] if tuple(l[0] for l in v['limits']) == tuple(variables)]
    positions=bound['positions']; width=float(bound['abel']['width'].subs(r.regulator,setting['regulator']))
    cut=setting['momentumBound']; tolerance=1e-3 if data['smoke'] else 1e-10
    cache={}; calls=hits=0; size=len(rows)*len(positions); point_artifacts=[]
    def packed(g): return np.concatenate((g['values'].ravel(),g['measureMutationValues'].ravel(),[g['quadratureMass']]))
    def function(k):
        nonlocal calls,hits
        calls+=1; key=float(k).hex()
        if key in cache:
            hits+=1; path,h=cache[key]; require(digest(path)==h,'conditional point cache hash'); return packed(unpickle(path)['group'])
        index=len(cache); worker.fixed_outer=float(k)
        def partial(state):
            if state['batchCount']%64==0:
                p=directory/'partials'/f'{layout}-{index:05}-{state["batchCount"]:07}.pickle';p.parent.mkdir(exist_ok=True)
                atomic_pickle(p,{'layout':layout,'outerValue':float(k),'state':state});artifacts.append(artifact(directory,p));save(directory/'partial-inventory.json',artifacts)
        group=worker.group(ti,variables,setting,bound['pairs'],width,positions,partial)
        packet={'task':task,'layout':layout,'outerValue':float(k),'outerUnit':bound['momentumUnit'],
                'conditionalUnits':[tuple(x-y for x,y in zip(v['unit'],bound['momentumUnit'])) for v in rows],
                'conditionalMassUnit':tuple((layout-1)*v for v in bound['momentumUnit']), 'setting':setting,'group':group,
                'sourceFiles':data['sourceFiles'],'provenance':data['provenance'],'methodJoins':conditional.METHOD_JOINS}
        p=directory/'points'/f'{layout}-{index:05}.pickle';p.parent.mkdir(exist_ok=True);atomic_pickle(p,packet)
        a=artifact(directory,p);artifacts.append(a);point_artifacts.append(a);save(directory/f'point-inventory-{layout}.json',point_artifacts)
        require(abs(group['volumeResidual'])<1e-10*(1+abs(group['boxVolume'])) and np.isfinite(group['values']).all() and np.isfinite(group['measureMutationValues']).all(),'finite conditional sum/mass')
        cache[key]=(p,a['sha256']);return packed(group)
    value,error,info=quad_vec(function,-cut,cut,epsabs=tolerance,epsrel=0.,norm='max',quadrature='gk21',workers=1,cache_size=8*1024*1024,limit=1000,full_output=True)
    shape=(len(rows),len(positions));g={'test':ti,'variables':tuple(variables),'rowIndices':[v['index'] for v in rows],
        'values':value[:size].reshape(shape),'measureMutationValues':value[size:2*size].reshape(shape),
        'measureMutationResidual':(value[size:2*size]-value[:size]).reshape(shape),'quadratureMass':value[-1],
        'boxVolume':(2*cut)**layout,'interval':(-cut,cut),'unitFrameErrorEstimate':float(error),'unitFrameTolerance':tolerance,
        'relativeTolerance':0.,'evaluations':int(info.neval),'success':bool(info.success),'status':int(info.status),'message':info.message,
        'intervals':info.intervals,'intervalPackedValues':info.integrals,'unitFrameIntervalErrors':info.errors,
        'calls':calls,'cacheHits':hits,'uniquePoints':len(cache),'conditionalNodeCount':sum(unpickle(directory/a['path'])['group']['nodeCount'] for a in point_artifacts)}
    g['volumeResidual']=g['quadratureMass']-g['boxVolume']
    path=directory/f'group-{layout}.pickle';atomic_pickle(path,{'setting':setting,'group':g,'pointArtifacts':point_artifacts});artifacts.append(artifact(directory,path));save(directory/'group-inventory.json',artifacts)
    require(g['success'] and g['status']==0 and np.isfinite(error) and error<=tolerance and np.isfinite(value).all(),'unfinished adaptive outer quadrature')
    require(len(g['intervals'])==len(g['intervalPackedValues'])==len(g['unitFrameIntervalErrors']) and np.isfinite(g['intervalPackedValues']).all(),'completed adaptive interval evidence')
    intervals=sorted(map(tuple,g['intervals']));require(intervals[0][0]==-cut and intervals[-1][1]==cut and all(a[1]==b[0] for a,b in zip(intervals,intervals[1:])),'adaptive interval partition')
    require(np.max(abs(np.sum(g['intervalPackedValues'],axis=0)-value)) < 1e-12*(1+np.max(abs(value))), 'adaptive saved interval contraction')
    return g


def guard_record(data,task,item,previous=None):
    wide.guard_record(data,task,item,previous)
    require(item['settings']==settings(data['refs'][task],item['stage'],data['smoke']),'adaptive settings identity')
    if item['stage']:
        n=item['evaluatedLayout'];g=next(v for v in item['groups'] if len(v['variables'])==n)
        require(item['method']=='adaptive' and g['success'] and g['status']==0 and g['unitFrameErrorEstimate']<=g['unitFrameTolerance'],'adaptive completion')
        require(g['relativeTolerance']==0. and g['interval']==(-item['settings'][n]['momentumBound'],item['settings'][n]['momentumBound']),'adaptive outer domain')


def evaluate_task(directory,data,task,*,max_batches=None,resume_state=None):
    require(max_batches is None and resume_state is None and task in TASKS,'wide adaptive task')
    started=time.monotonic();ti,di=task;bound=data['domains'][di];ref=data['refs'][task]
    worker=conditional.ConditionalMomentum(bound['rows'],bound['sources'],data['r'],outer_value=0.)
    records=[];artifacts=[]
    for stage in (0,1,3):
        previous=records[-1] if records else None;groups=list(ref['groups'] if previous is None else previous['groups'])
        if stage:
            n=1 if stage==1 else 2;s=settings(ref,stage,data['smoke'])[n];old=next(v for v in groups if len(v['variables'])==n)
            fresh=adaptive_group(directory,data,task,worker,old['variables'],s,artifacts)
            groups=[fresh if len(g['variables'])==n else g for g in groups]
        item=native_make_record(data,task,stage,groups,previous,copy.deepcopy(parallel.cache_state(worker)) if stage else None)
        if stage:item['method']='adaptive'
        p=directory/'records'/f'{stage}.pickle';p.parent.mkdir(exist_ok=True);atomic_pickle(p,item);artifacts.append(artifact(directory,p));save(directory/'record-inventory.json',artifacts)
        guard_record(data,task,item,previous);records.append(item)
    packet={'task':task,'records':records,'artifacts':artifacts,'sourceFiles':data['sourceFiles'],'provenance':data['provenance'],
            'methodJoins':conditional.METHOD_JOINS,'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    atomic_pickle(directory/'result.pickle',packet)
    require(data['sourceFiles']=={n:digest(ROOT/n) for n in data['sourceFiles']},'worker source identity')
    save(directory/'checks.json',{'task':task,'resultSha256':digest(directory/'result.pickle'),'wallSeconds':packet['wallSeconds'],'peakRssKiB':packet['peakRssKiB']})
    return packet


def dispatch(base,data):
    previous=parallel.evaluate_task;parallel.evaluate_task=evaluate_task
    try:return parallel.dispatch(base,data,TASKS)
    finally:parallel.evaluate_task=previous


def derived_emitter():
    old=ast.parse(textwrap.dedent(inspect.getsource(position.emit)));tree=copy.deepcopy(old)
    labels=[n for n in ast.walk(tree) if isinstance(n,ast.Constant) and n.value=='LAYOUT_ORDERS_']
    require(len(labels)==1,'native layout-order label census');labels[0].value='GAUSS_REFERENCE_LAYOUT_ORDERS_'
    restored=copy.deepcopy(tree);next(n for n in ast.walk(restored) if isinstance(n,ast.Constant) and n.value=='GAUSS_REFERENCE_LAYOUT_ORDERS_').value='LAYOUT_ORDERS_'
    require(ast.dump(restored)==ast.dump(old),'native emitter label-only AST join')
    namespace=dict(vars(position),PREFIX=PREFIX);exec(compile(ast.fix_missing_locations(tree),'<wide-adaptive-reference-label>','exec'),namespace)
    proxy=types.SimpleNamespace(**(dict(vars(position)) | {'emit':namespace['emit']}))
    fn=types.FunctionType(wide.emit.__code__,dict(vars(wide),PREFIX=PREFIX,position=proxy),wide.emit.__name__)
    return fn,{'nativeEmitterAstSha256':hashlib.sha256(ast.dump(old).encode()).hexdigest(),'restoredEmitterAstSha256':hashlib.sha256(ast.dump(restored).encode()).hexdigest()}


COMMON_EMIT, EMITTER_JOIN=derived_emitter()


def emit(result,r,bound,rows,finest,provenance):
    COMMON_EMIT(result,r,bound,rows,finest,provenance)
    metadata=engine.FullPencilModes.__new__(engine.FullPencilModes);metadata.r=r;zero=engine.PHYSICAL_METADATA.dimensions.zero;mom=bound['momentumUnit']
    def numeric(name,value,unit):
        body=domain.momentum.number(value);engine.emit(PREFIX+'_'+name,body);engine.emit('METADATA_'+PREFIX+'_'+name,metadata.numeric_metadata(body,unit))
    for item in result['records']:
        if not item['stage']:continue
        n=item['evaluatedLayout'];g=next(v for v in item['groups'] if len(v['variables'])==n);key=f'{item["test"]}_{item["index"]}';size=len(g['rowIndices'])*len(bound['positions'])
        engine.physical(PREFIX+'_OUTER_METHOD_'+key,'adaptive GK21; epsrel=0; displayed Gauss outer order is reference only; inner/source/profile rules held')
        numeric('ADAPTIVE_INTERVALS_'+key,g['intervals'],lambda p:mom)
        numeric('ADAPTIVE_INTERVAL_ERRORS_'+key,g['unitFrameIntervalErrors'],lambda p:zero)
        numeric('ADAPTIVE_ERROR_STATUS_'+key,(g['unitFrameErrorEstimate'],g['unitFrameTolerance'],g['relativeTolerance'],g['status'],int(g['success']),g['evaluations'],g['calls'],g['cacheHits'],g['uniquePoints'],g['conditionalNodeCount']),lambda p:zero)
        numeric('ADAPTIVE_INTERVAL_VALUES_'+key,g['intervalPackedValues'],lambda p:tuple(n*v for v in mom) if p[1]==2*size else bound['rows'][g['rowIndices'][(p[1]%size)//len(bound['positions'])]]['unit'])


def replay_adapter():
    cache=types.FunctionType(position.cache_for_tag.__code__,dict(vars(position),PREFIX=PREFIX),position.cache_for_tag.__name__)
    factory=types.FunctionType(position.replay_adapter.__code__,dict(vars(position),PREFIX=PREFIX,emit=emit,cache_for_tag=cache),position.replay_adapter.__name__)
    return factory()


def finish(base,data,workers,started,focused=None):
    records=[];artifacts=[];inventory=[]
    manifest=json.loads((base/'workers.json').read_text())
    require(manifest['status']=='completed' and len(manifest['outcomes'])==4 and all(v['exitCode']==v['stderrBytes']==0 for v in manifest['outcomes']),'clean adaptive workers')
    for task,w in sorted(workers.items()):
        directory=base/f'worker-{task[0]}-{task[1]}'
        require(w['task']==task and w['sourceFiles']==data['sourceFiles'] and w['provenance']==data['provenance'] and w['methodJoins']==conditional.METHOD_JOINS,'worker packet joins')
        require(digest(directory/'result.pickle')==manifest['resultFiles'][f'{task[0]}-{task[1]}'],'worker manifest hash');artifacts.append(artifact(base,directory/'result.pickle'))
        for a in w['artifacts']:
            require(digest(directory/a['path'])==a['sha256'],'worker record/point/partial hash');artifacts.append(dict(a,path=str((directory/a['path']).relative_to(base))))
            if a['path'].startswith('points/'):
                point=unpickle(directory/a['path']);n=point['layout'];bound=data['domains'][task[1]];g=point['group']
                require(point['task']==task and point['sourceFiles']==data['sourceFiles'] and point['provenance']==data['provenance'] and point['methodJoins']==conditional.METHOD_JOINS,'point provenance')
                require(point['setting']==settings(data['refs'][task],1 if n==1 else 3,data['smoke'])[n] and point['outerUnit']==bound['momentumUnit'],'point setting/unit')
                require(g['test']==task[0] and g['rowIndices']==[v['index'] for v in bound['rows'] if tuple(l[0] for l in v['limits'])==tuple(g['variables'])],'point native rows/field/limits')
                require(point['conditionalUnits']==[tuple(x-y for x,y in zip(bound['rows'][i]['unit'],bound['momentumUnit'])) for i in g['rowIndices']] and point['conditionalMassUnit']==tuple((n-1)*v for v in bound['momentumUnit']),'point full conditional units')
        previous=None
        for item in w['records']:
            guard_record(data,task,item,previous);previous=item;records.append(item)
            p=base/'records'/f'{item["test"]}-{item["index"]}.pickle';p.parent.mkdir(exist_ok=True);atomic_pickle(p,item);inventory.append(artifact(base,p))
    save(base/'record-inventory.json',inventory)
    result={'records':records,'joins':data['joins'],'recordArtifacts':inventory,'workerArtifacts':artifacts}
    atomic_pickle(base/'wide-adaptive.pickle',{'result':result,'sourceFiles':data['sourceFiles'],'provenance':data['provenance'],'dimensionState':dict(vars(engine.PHYSICAL_METADATA.dimensions))})
    before={p.name:digest(p) for p in base.glob('*.pickle')};replay,join=replay_adapter();engine.EMISSION_LINES.clear();engine.PAYLOAD_ENCODER=engine.PayloadEncoder()
    entries,keys,paths=replay(base,result,data['r'],data['domains'][2],data['domains'][2]['rows'],data['refs'],data['provenance'])
    norms=[]
    for item in records:
        if not item['stage']:continue
        g=next(v for v in item['groups'] if len(v['variables'])==item['evaluatedLayout'])
        norms.append({'test':item['test'],'domain':item['domain'],'layout':item['evaluatedLayout'],
                      'rawIntegralChanges':[float(np.max(abs(v))) for v in item['integralChanges']],
                      'actionChange':float(np.max(abs(item['precedingActionDifference']))),
                      'error':g['unitFrameErrorEstimate'],'tolerance':g['unitFrameTolerance'],'evaluations':g['evaluations'],'points':g['uniquePoints']})
    summary={'runDirectory':str(base),'sourceFiles':data['sourceFiles'],'provenance':data['provenance'],'smoke':data['smoke'],'records':len(records),'norms':norms,
             'tagCount':len(entries),'writeKeyCount':len(keys),'metadataPaths':paths,'replayMethodJoin':join,'conditionalMethodJoins':conditional.METHOD_JOINS,'emitterJoin':EMITTER_JOIN,
             'focused':focused,'recordArtifacts':inventory,'workerArtifacts':artifacts,'workerManifest':manifest,
             'packetHashesBeforeEmission':before,'packetHashesAfterEmission':{n:digest(base/n) for n in before},
             'artifacts':{p.name:{'bytes':p.stat().st_size,'sha256':digest(p)} for p in base.iterdir() if p.suffix in ('.pickle','.out')},
             'wallSeconds':time.monotonic()-started,'coordinatorPeakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
             'scope':'Finite independent single/pair outer integration at fixed wide boxes with accepted three-momentum values held; no uniform or independent-grade convergence, global exceptional coverage, physical tails, Abel limit, scattering or poles.'}
    for a in (*inventory,*artifacts): require(digest(base/a['path'])==a['sha256'],'final saved artifact')
    require(before==summary['packetHashesAfterEmission'] and data['sourceFiles']=={n:digest(ROOT/n) for n in data['sourceFiles']} and not engine.PHYSICAL_METADATA.dimensions.constraints,'final source/packet/dimension guard')
    save(base/'checks.json',summary);return summary


def focused_checks(base,data):
    evidence=[]
    for task in TASKS:
        ti,di=task;bound=data['domains'][di];ref=data['refs'][task];r=data['r']
        actual=conditional.ConditionalMomentum(bound['rows'],bound['sources'],r,outer_value=0.)
        native=engine.BoundedSourceFourierQuadrature.CachedMomentum(bound['rows'],bound['sources'],r)
        for n in (1,2):
            setting=ref['layoutSettings'][n];variables=next(g['variables'] for g in ref['groups'] if len(g['variables'])==n);width=float(bound['abel']['width'].subs(r.regulator,setting['regulator']))
            nodes,_,_=native.rule(-setting['momentumBound'],setting['momentumBound'],setting['outerOrder'])
            for pi,k in enumerate((nodes[0],nodes[len(nodes)//2])):
                actual.fixed_outer=float(k);g=actual.group(ti,variables,setting,bound['pairs'],width,bound['positions'])
                expected=native.at_momentum(ti,variables[0],setting,bound['positions'],k) if n==1 else native.inner_at_outer(ti,variables,setting,bound['pairs'],width,bound['positions'],k)
                batches=list(actual.batches(variables,setting,bound['pairs'],width));points=np.concatenate([v[0] for v in batches]);weights=np.concatenate([v[1] for v in batches])
                if n==1:ep=np.asarray([[k]]);ew=np.ones(1)
                else:
                    inner,outer=variables;centers=[k for a,b in bound['pairs'] if inner in (a,b) and (b if a==inner else a)==outer]
                    x,ew,_=native.rule(-setting['momentumBound'],setting['momentumBound'],setting['panelOrder'] if centers else setting['outerOrder'],centers,width);ep=np.column_stack((x,np.full(len(x),k)))
                proof={'task':task,'layout':n,'outerValue':float(k),'setting':setting,'actual':g,'expectedValues':expected,'residual':g['values']-expected,'actualPoints':points,'actualWeights':weights,'expectedPoints':ep,'expectedWeights':ew,
                       'valueUnits':[tuple(x-y for x,y in zip(bound['rows'][i]['unit'],bound['momentumUnit'])) for i in g['rowIndices']],'pointUnit':bound['momentumUnit'],'weightUnit':tuple((n-1)*x for x in bound['momentumUnit'])}
                p=base/f'focused-{ti}-{di}-{n}-{pi}.pickle';atomic_pickle(p,proof);evidence.append(artifact(base,p))
                require(np.array_equal(points,ep) and np.array_equal(weights,ew),'conditional native nodes/weights')
                require(np.max(abs(proof['residual']))<1e-12*(1+np.max(abs(expected))) and abs(g['volumeResidual'])<1e-10*(1+abs(g['boxVolume'])),'conditional native values/mass')
                require(not np.any((abs(g['values'])>1e-9)&(abs(g['measureMutationResidual'])<1e-12)),'conditional actual measure sensitivity')
    return {'conditionalArtifacts':evidence,'conditionalMethodJoins':conditional.METHOD_JOINS,'emitterJoin':EMITTER_JOIN}


def coarse_checks(base,data,workers):
    artifacts=[]
    for task,w in sorted(workers.items()):
        ti,di=task;bound=data['domains'][di];r=data['r']
        for item in w['records']:
            if not item['stage']:continue
            values={i:g['values'][j] for g in item['groups'] for j,i in enumerate(g['rowIndices'])};action=np.zeros_like(item['action']);compiler=engine.BoundedActionQuadrature({})
            for cell in (v for v in bound['cells'] if v['test']==ti):
                for pi,z in enumerate(bound['positions']):
                    env={r.z:z,r.regulator:.2};action[pi,cell['column'],cell['row']]=complex(compiler.evaluate(cell['local'],env))+sum(complex(compiler.evaluate(coef,env))*values[row][pi] for row,coef in cell['terms'])
            proof={'task':task,'stage':item['stage'],'directAction':action,'workerAction':item['action'],'residual':action-item['action'],'equationUnits':bound['equationUnits']}
            p=base/f'coarse-{ti}-{di}-{item["stage"]}.pickle';atomic_pickle(p,proof);artifacts.append(artifact(base,p));require(not np.any(proof['residual']),'independent full native-cell action')
    return artifacts


def main():
    p=argparse.ArgumentParser();p.add_argument('--run-directory',type=Path,required=True);p.add_argument('--preflight',action='store_true');p.add_argument('--focused-only',action='store_true');args=p.parse_args()
    base=args.run_directory.resolve();base.relative_to(STORE);base.mkdir(parents=True,exist_ok=False);started=time.monotonic();data=load(base)
    focused=focused_checks(base,data) if args.preflight or args.focused_only else None
    if focused is not None:save(base/'focused-before-workers.json',focused)
    if args.focused_only:print(json.dumps(focused,indent=2));return
    data['smoke']=args.preflight;workers=dispatch(base,data)
    if args.preflight:focused['coarseArtifacts']=coarse_checks(base,data,workers);save(base/'focused-after-workers.json',focused)
    print(json.dumps(finish(base,data,workers,started,focused),indent=2))


if __name__=='__main__':main()
