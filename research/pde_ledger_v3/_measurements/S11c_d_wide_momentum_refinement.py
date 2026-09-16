#!/usr/bin/env python3
"""Refine native one-/two-momentum rules at fixed accepted wide cutoffs."""
import argparse, copy, json, resource, shutil, time, types
from pathlib import Path
import numpy as np
import sympy as sp
import S11c_d_momentum_domain_action as domain

ROOT, STORE = domain.ROOT, domain.STORE
M = ROOT / '_measurements'
CHECKPOINT = M / 'S11c_d_momentum_domain_action_checkpoint.json'
PLAN = M / 'S11c_d_wide_momentum_refinement_plan.md'
PREFIX = 'WIDE_MOMENTUM_REFINEMENT_LAB_HELD_RHO4_CONSTANT'
engine, position, parallel = domain.engine, domain.position, domain.parallel
require, digest, save, atomic_pickle, unpickle = domain.require, domain.digest, domain.save, domain.atomic_pickle, domain.unpickle
artifact = domain.artifact
TASKS = ((0, 2), (0, 3), (1, 2), (1, 3))


def settings(reference, stage, smoke=False):
    rules = copy.deepcopy(reference['layoutSettings'])
    if stage >= 1:
        rules[1]['outerOrder'] = (3 if stage == 1 else 4) if smoke else (324 if stage == 1 else 432)
    if stage >= 3:
        rules[2]['outerOrder'] = (3 if stage == 3 else 4) if smoke else (324 if stage == 3 else 432)
        inner = (1 if stage < 5 else stage - 3) if smoke else (32 if stage < 5 else 48 if stage == 5 else 64)
        rules[2].update(panelOrder=inner, innerOrders=(inner,))
    if smoke and stage:
        # Only newly evaluated layouts use smoke rules; held groups retain their
        # accepted settings. Source/profile changes here are declared smoke only.
        for n in ((1,) if stage < 3 else (1, 2)):
            rules[n].update(sourceNodes=16, profileNodes=16)
    return rules


def load(base):
    accepted, previous = domain.momentum.source.accepted(CHECKPOINT)
    for n, h in accepted['sourceFiles'].items():
        require(digest(ROOT/n) == h == digest(previous/'source'/n), 'accepted current/frozen source')
    for category in ('recordArtifacts', 'workerArtifacts'):
        for a in accepted[category]:
            require(digest(previous/a['path']) == a['sha256'], 'accepted record/layout/partial hash')
    bp = unpickle(previous/'bound-momentum-domains.pickle')
    packet = unpickle(previous/'momentum-domain-action.pickle')
    require(bp['sourceFiles'] == packet['sourceFiles'] == accepted['sourceFiles'] and bp['provenance'] == packet['provenance'] == accepted['provenance'], 'accepted packet provenance')
    _, sb = domain.momentum.source.accepted(domain.momentum.source.SOURCE)
    r, dimensions = domain.momentum.source.native.source.restore_context(unpickle(sb/'reduced-action.pickle'))
    dimensions.__dict__.update(packet['dimensionState'])
    refs = {}
    for task in TASKS:
        ti, di = task
        rec = next(v for v in packet['result']['records'] if (v['test'], v['index']) == task)
        ref = dict(bp['references'][ti], action=rec['action'], groups=rec['groups'],
                   contributions=[{'index':v['index'], 'terms':v['values']} for v in rec['terms']],
                   layoutSettings=rec['settings'], setting=rec['settings'][3])
        rebuilt = position.assemble(r, bp['domains'][di], ref, ref['groups'])
        require(np.array_equal(rebuilt['action'], rec['action']) and all(not np.any(v['differences']) for v in rebuilt['terms']), 'accepted full baseline join')
        refs[task] = ref
    pins = dict(accepted['sourceFiles'])
    for path in (CHECKPOINT, PLAN, Path(__file__).resolve()):
        pins[str(path.relative_to(ROOT))] = digest(path)
    for n in pins:
        dest = base/'source'/n; dest.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(ROOT/n, dest)
    for name in ('bound-momentum-domains.pickle', 'momentum-domain-action.pickle'):
        shutil.copyfile(previous/name, base/('accepted-'+name))
        require(digest(base/('accepted-'+name)) == accepted['artifacts'][name]['sha256'], 'accepted packet byte copy')
    provenance = dict(accepted['provenance'], WIDE_REFINEMENT_INPUT_SHA256=digest(CHECKPOINT),
                      WIDE_REFINEMENT_INSTRUMENT_SHA256=digest(Path(__file__)), WIDE_REFINEMENT_PLAN_SHA256=digest(PLAN))
    data = {'r':r, 'domains':{i:bp['domains'][i] for i in (2,3)}, 'joins':{i:bp['joins'][i] for i in (2,3)},
            'refs':refs, 'sourceFiles':pins, 'provenance':provenance, 'previous':previous, 'smoke':False}
    atomic_pickle(base/'bound-wide-refinement.pickle', {k:v for k,v in data.items() if k != 'r'} | {'dimensionState':dict(vars(dimensions))})
    save(base/'preflight.json', {'sourceFiles':pins, 'provenance':provenance, 'tasks':TASKS,
                               'acceptedRunDirectory':str(previous), 'nativeEngineUnchanged':True})
    return data


def make_record(data, task, stage, groups, previous=None, telemetry=None):
    ti, di = task; ref = data['refs'][task]; bound = data['domains'][di]
    n = 0 if stage == 0 else 1 if stage <= 2 else 2
    item = {'test':ti, 'index':100*di+stage, 'domain':di, 'stage':stage,
            'method':'retained' if stage == 0 else 'gauss', 'evaluatedLayout':n,
            'groups':groups, 'settings':settings(ref, stage, data['smoke']),
            **position.assemble(data['r'], bound, ref, groups)}
    if previous is not None:
        item['precedingActionDifference'] = item['action'] - previous['action']
        item['integralChanges'] = [a['values']-b['values'] for a,b in zip(groups, previous['groups'])]
        item['termChanges'] = [np.asarray(a['values'])-np.asarray(b['values']) for a,b in zip(item['terms'], previous['terms'])]
        item['heldTermMask'] = [[len(bound['rows'][row]['limits']) != n for row in term['rows']] for term in item['terms']]
        item['heldTermDifferences'] = [change[np.asarray(mask, dtype=bool)] for change,mask in zip(item['termChanges'], item['heldTermMask'])]
    if telemetry is not None: item['telemetry'] = telemetry
    return item


def guard_record(data, task, item, previous=None):
    ti, di = task; bound = data['domains'][di]; n = item['evaluatedLayout']
    require({i for g in item['groups'] for i in g['rowIndices']} == set(range(80)), 'full action row census')
    for g in item['groups']:
        size = len(g['variables']); s = item['settings'][size]
        require(g['test'] == ti and g['rowIndices'] == [row['index'] for row in bound['rows'] if tuple(v[0] for v in row['limits']) == tuple(g['variables'])], 'native row/field/limit join')
        require(all(tuple(map(float, v[1:])) == (-s['momentumBound'], s['momentumBound']) for row in bound['rows'] if len(row['limits']) == size for v in row['limits']), 'native finite limits')
        require(np.isfinite(g['values']).all() and np.isfinite(g['measureMutationValues']).all() and abs(g['volumeResidual']) < 1e-10*(1+abs(g['boxVolume'])), 'finite group/mass')
        require(np.array_equal(g['measureMutationResidual'], g['measureMutationValues']-g['values']) and not np.any((abs(g['values']) > 1e-9) & (abs(g['measureMutationResidual']) < 1e-12)), 'actual measure control')
        if previous is not None and size != n:
            old = next(v for v in previous['groups'] if len(v['variables']) == size)
            for key in ('values','measureMutationValues','measureMutationResidual'):
                require(np.array_equal(g[key], old[key]), 'held group identity')
            require(s == previous['settings'][size], 'held setting identity')
    if previous is not None:
        require(all(not np.any(v) for v in item['heldTermDifferences']), 'held native terms')
    if 'telemetry' in item:
        t = item['telemetry']; rows = [v for v in bound['rows'] if len(v['limits']) <= max(n,1)]
        require(set(t['sourceFrequencyCensus']) == {(ti,f['sourceIndex']) for row in rows for f in row['factors']}, 'evaluated source census')
        profiles = set().union(*(f['coefficient'].atoms(sp.Integral) for row in rows for f in row['factors']))
        require(t['evaluatedProfileIntegrals'] == profiles and profiles <= set(bound['profileUnits']), 'evaluated profile census')


def evaluate_task(directory, data, task, *, max_batches=None, resume_state=None):
    require(max_batches is None and resume_state is None and task in TASKS, 'wide refinement task')
    started = time.monotonic(); ti, di = task; bound = data['domains'][di]; ref = data['refs'][task]
    worker = engine.BoundedSourceFourierQuadrature.ThreeMomentum(bound['rows'], bound['sources'], data['r'])
    records = []; artifacts = []
    for stage in range(7):
        previous = records[-1] if records else None
        groups = list(ref['groups'] if previous is None else previous['groups'])
        if stage:
            n = 1 if stage <= 2 else 2; s = settings(ref, stage, data['smoke'])[n]
            old = next(v for v in groups if len(v['variables']) == n)
            width = float(bound['abel']['width'].subs(data['r'].regulator, s['regulator']))
            def partial(state):
                if state['batchCount'] % 64 == 0:
                    path = directory/'partials'/f'{stage}-{state["batchCount"]:07}.pickle'; path.parent.mkdir(exist_ok=True)
                    atomic_pickle(path, state); artifacts.append(artifact(directory,path)); save(directory/'partial-inventory.json',artifacts)
            fresh = worker.group(ti, old['variables'], s, bound['pairs'], width, bound['positions'], partial)
            path = directory/f'group-{stage}.pickle'; atomic_pickle(path, {'setting':s, 'group':fresh}); artifacts.append(artifact(directory,path))
            groups = [fresh if len(g['variables']) == n else g for g in groups]
        item = make_record(data, task, stage, groups, previous, copy.deepcopy(parallel.cache_state(worker)) if stage else None)
        path = directory/'records'/f'{stage}.pickle'; path.parent.mkdir(exist_ok=True); atomic_pickle(path,item)
        artifacts.append(artifact(directory,path)); save(directory/'record-inventory.json',artifacts)
        guard_record(data, task, item, previous); records.append(item)
    packet = {'task':task, 'records':records, 'artifacts':artifacts, 'sourceFiles':data['sourceFiles'],
              'provenance':data['provenance'], 'wallSeconds':time.monotonic()-started,
              'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    atomic_pickle(directory/'result.pickle',packet)
    require(data['sourceFiles'] == {n:digest(ROOT/n) for n in data['sourceFiles']}, 'worker source identity')
    save(directory/'checks.json', {'task':task,'resultSha256':digest(directory/'result.pickle'),
                                  'wallSeconds':packet['wallSeconds'],'peakRssKiB':packet['peakRssKiB']})
    return packet


def dispatch(base,data):
    original = parallel.evaluate_task; parallel.evaluate_task = evaluate_task
    try: return parallel.dispatch(base,data,TASKS)
    finally: parallel.evaluate_task = original


def emit(result,r,bound,rows,finest,provenance):
    fn = types.FunctionType(position.emit.__code__,dict(vars(position),PREFIX=PREFIX),position.emit.__name__)
    fn(dict(result,joins={}),r,bound,rows,finest,provenance)
    metadata = engine.FullPencilModes.__new__(engine.FullPencilModes); metadata.r = r
    zero = engine.PHYSICAL_METADATA.dimensions.zero
    def numeric(name,value,unit):
        body = domain.momentum.number(value); engine.emit(PREFIX+'_'+name,body)
        engine.emit('METADATA_'+PREFIX+'_'+name,metadata.numeric_metadata(body,unit))
    for di,join in result['joins'].items():
        for row in join['rows']:
            numeric(f'NATIVE_LIMITS_{di}_{row["row"]}',row['new'],lambda p:bound['momentumUnit'])
    for item in result['records']:
        key=f'{item["test"]}_{item["index"]}'
        numeric('RECORD_DOMAIN_STAGE_LAYOUT_'+key,(item['domain'],item['stage'],item['evaluatedLayout']),lambda p:zero)
        if 'heldTermDifferences' in item:
            for i,(mask,difference) in enumerate(zip(item['heldTermMask'],item['heldTermDifferences'])):
                numeric(f'HELD_TERM_MASK_{key}_{i}',list(map(int,mask)),lambda p:zero)
                numeric(f'HELD_TERM_DIFFERENCE_{key}_{i}',difference,lambda p:bound['equationUnits'][item['terms'][i]['index'][2]])


def replay_adapter():
    cache = types.FunctionType(position.cache_for_tag.__code__,dict(vars(position),PREFIX=PREFIX),position.cache_for_tag.__name__)
    factory = types.FunctionType(position.replay_adapter.__code__,dict(vars(position),PREFIX=PREFIX,emit=emit,cache_for_tag=cache),position.replay_adapter.__name__)
    return factory()


def finish(base,data,workers,started,focused=None):
    records=[]; artifacts=[]; inventory=[]
    for task,w in sorted(workers.items()):
        directory=base/f'worker-{task[0]}-{task[1]}'
        require(w['task']==task and w['sourceFiles']==data['sourceFiles'] and w['provenance']==data['provenance'],'worker packet joins')
        artifacts.append(artifact(base,directory/'result.pickle'))
        for a in w['artifacts']:
            require(digest(directory/a['path'])==a['sha256'],'worker record/group/partial hash')
            artifacts.append(dict(a,path=str((directory/a['path']).relative_to(base))))
        previous=None
        for item in w['records']:
            guard_record(data,task,item,previous); previous=item; records.append(item)
            path=base/'records'/f'{item["test"]}-{item["index"]}.pickle';path.parent.mkdir(exist_ok=True);atomic_pickle(path,item);inventory.append(artifact(base,path))
    save(base/'record-inventory.json',inventory)
    result={'records':records,'joins':data['joins'],'recordArtifacts':inventory,'workerArtifacts':artifacts}
    atomic_pickle(base/'wide-refinement.pickle',{'result':result,'sourceFiles':data['sourceFiles'],'provenance':data['provenance'],'dimensionState':dict(vars(engine.PHYSICAL_METADATA.dimensions))})
    before={p.name:digest(p) for p in base.glob('*.pickle')}; replay,join=replay_adapter()
    engine.EMISSION_LINES.clear();engine.PAYLOAD_ENCODER=engine.PayloadEncoder()
    entries,keys,paths=replay(base,result,data['r'],data['domains'][2],data['domains'][2]['rows'],data['refs'],data['provenance'])
    norms=[{'test':v['test'],'domain':v['domain'],'stage':v['stage'],'evaluatedLayout':v['evaluatedLayout'],
            'rawIntegralChanges':[float(np.max(abs(a))) for a in v['integralChanges']],
            'actionChange':float(np.max(abs(v['precedingActionDifference']))),
            'rawTermChange':max(float(np.max(abs(a))) if len(a) else 0. for a in v['termChanges'])} for v in records if v['stage']]
    summary={'runDirectory':str(base),'sourceFiles':data['sourceFiles'],'provenance':data['provenance'],'smoke':data['smoke'],
             'records':len(records),'norms':norms,'tagCount':len(entries),'writeKeyCount':len(keys),'metadataPaths':paths,
             'replayMethodJoin':join,'focused':focused,'recordArtifacts':inventory,'workerArtifacts':artifacts,
             'workerManifest':json.loads((base/'workers.json').read_text()),'packetHashesBeforeEmission':before,
             'packetHashesAfterEmission':{n:digest(base/n) for n in before},
             'artifacts':{p.name:{'bytes':p.stat().st_size,'sha256':digest(p)} for p in base.iterdir() if p.suffix in ('.pickle','.out')},
             'wallSeconds':time.monotonic()-started,'coordinatorPeakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
             'scope':'Fixed wide boxes: single and paired quadrature changes with accepted three-momentum values explicitly held. No independent adaptive check, full wide-box convergence, physical tails, Abel limit, scattering or poles.'}
    for a in (*inventory,*artifacts):require(digest(base/a['path'])==a['sha256'],'final saved artifact')
    require(before==summary['packetHashesAfterEmission'] and data['sourceFiles']=={n:digest(ROOT/n) for n in data['sourceFiles']} and not engine.PHYSICAL_METADATA.dimensions.constraints,'final source/packet/dimension guard')
    save(base/'checks.json',summary);return summary


def focused_checks(base,data):
    evidence=[]
    class PrefixComplete(Exception): pass
    for task in TASKS:
        ti,di=task;ref=data['refs'][task];bound=data['domains'][di]
        for n in (1,2):
            s=ref['layoutSettings'][n];variables=next(g['variables'] for g in ref['groups'] if len(g['variables'])==n)
            worker=engine.BoundedSourceFourierQuadrature.ThreeMomentum(bound['rows'],bound['sources'],data['r']);captured=[]
            def stop(state):
                if state['batchCount']==64:captured.append(state);raise PrefixComplete()
            width=float(bound['abel']['width'].subs(data['r'].regulator,s['regulator']))
            try:
                group=worker.group(ti,variables,s,bound['pairs'],width,bound['positions'],stop)
                expected=next(g for g in ref['groups'] if len(g['variables'])==n); actual=group
            except PrefixComplete:
                actual=captured[0];expected=unpickle(data['previous']/f'worker-{ti}-{di}'/'partials'/f'{n}-0000064.pickle')
            proof={'task':task,'layout':n,'actual':actual,'expected':expected,'valueResidual':actual['values']-expected['values'],
                   'mutationResidual':actual['measureMutationValues']-expected['measureMutationValues'],'massResidual':actual['quadratureMass']-expected['quadratureMass']}
            path=base/f'prefix-{ti}-{di}-{n}.pickle';atomic_pickle(path,proof);evidence.append(artifact(base,path))
            require(not np.any(proof['valueResidual']) and not np.any(proof['mutationResidual']) and proof['massResidual']==0 and actual['nodeCount']==expected['nodeCount'] and actual['rowIndices']==expected['rowIndices'],'actual accepted native prefix')
        for stage in range(1,7):
            before,after=settings(ref,stage-1),settings(ref,stage)
            changes={(n,k) for n in before for k in set(before[n])|set(after[n]) if before[n].get(k)!=after[n].get(k)}
            expected={(1 if stage<=2 else 2,'outerOrder')} if stage<=4 else {(2,'innerOrders'),(2,'panelOrder')}
            require(changes==expected,'isolated physical quadrature coordinate')
    return {'prefixArtifacts':evidence,'productionSettingIsolation':True}


def coarse_checks(base,data,workers):
    evidence=[]
    for task,w in sorted(workers.items()):
        ti,di=task;bound=data['domains'][di];r=data['r']
        worker=engine.BoundedSourceFourierQuadrature.ThreeMomentum(bound['rows'],bound['sources'],r)
        for item in w['records']:
            if not item['stage']:continue
            n=item['evaluatedLayout'];g=next(v for v in item['groups'] if len(v['variables'])==n);s=item['settings'][n]
            fresh=worker.group(ti,g['variables'],s,bound['pairs'],float(bound['abel']['width'].subs(r.regulator,s['regulator'])),bound['positions'])
            require(np.array_equal(fresh['values'],g['values']) and np.array_equal(fresh['measureMutationValues'],g['measureMutationValues']) and fresh['quadratureMass']==g['quadratureMass'] and fresh['nodeCount']==g['nodeCount'],'coarse worker/serial group')
            values={i:g['values'][j] for g in item['groups'] for j,i in enumerate(g['rowIndices'])};action=np.zeros_like(item['action']);compiler=engine.BoundedActionQuadrature({})
            for cell in (v for v in bound['cells'] if v['test']==ti):
                for pi,z in enumerate(bound['positions']):
                    env={r.z:z,r.regulator:.2};action[pi,cell['column'],cell['row']]=complex(compiler.evaluate(cell['local'],env))+sum(complex(compiler.evaluate(coef,env))*values[row][pi] for row,coef in cell['terms'])
            proof={'task':task,'stage':item['stage'],'serialGroup':fresh,'workerGroup':g,'directAction':action,'workerAction':item['action'],'residual':action-item['action']}
            path=base/f'coarse-{ti}-{di}-{item["stage"]}.pickle';atomic_pickle(path,proof);evidence.append(artifact(base,path));require(not np.any(proof['residual']),'independent native-cell contraction')
    return evidence


def main():
    p=argparse.ArgumentParser();p.add_argument('--run-directory',type=Path,required=True);p.add_argument('--preflight',action='store_true');p.add_argument('--focused-only',action='store_true');a=p.parse_args()
    base=a.run_directory.resolve();base.relative_to(STORE);base.mkdir(parents=True,exist_ok=False);started=time.monotonic();data=load(base)
    focused=focused_checks(base,data) if a.preflight or a.focused_only else None
    if focused is not None:save(base/'focused-before-workers.json',focused)
    if a.focused_only:print(json.dumps(focused,indent=2));return
    data['smoke']=a.preflight;workers=dispatch(base,data)
    if a.preflight:focused['coarseArtifacts']=coarse_checks(base,data,workers);save(base/'focused-after-workers.json',focused)
    summary=finish(base,data,workers,started,focused);print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
