#!/usr/bin/env python3
"""Independent outer quadrature with accepted wide-triple inner rules."""
import argparse, copy, hashlib, json, shutil, time, types
from pathlib import Path
import numpy as np
import S11c_d_wide_momentum_adaptive as adaptive
import S11c_d_wide_three_momentum as triple

ROOT, STORE, M = triple.ROOT, triple.STORE, triple.M
engine, position, parallel, domain = triple.engine, triple.position, triple.parallel, triple.domain
require, digest, save, atomic_pickle, unpickle, artifact = (
    triple.require, triple.digest, triple.save, triple.atomic_pickle, triple.unpickle, triple.artifact)
conditional = adaptive.conditional
TASKS = triple.TASKS
CHECKPOINT = M/'S11c_d_wide_three_momentum_checkpoint.json'
PLAN = M/'S11c_d_wide_three_adaptive_plan.md'
AUTHORITY = ROOT/'directives/S11c_d_NONLINEAR_POLE_CONTRACT.md'
REPAIR = M/'S11c_d_nonlinear_pole_repair_checkpoint.json'
PREFIX = 'WIDE_THREE_ADAPTIVE_LAB_HELD_RHO4_CONSTANT'
DISPOSITION = ('The shared v10 baseline and nonlinearPoleV2 addendum jointly govern new work. '
               'This fixed-box action quadrature calls no pole constructor; all physical '
               'operands, numerical inputs, native engine and prior frozen artifacts remain unchanged.')


def settings(reference, stage, smoke=False):
    require(stage in (0,3), 'wide triple adaptive stage')
    rules = copy.deepcopy(reference['layoutSettings'])
    if smoke and stage:
        rules[3].update(innerOrders=(1,1), sourceNodes=16, profileNodes=16)
    return rules


make_record, MAKE_JOIN = triple.derived(adaptive.wide.make_record, [
    ('n = 0 if stage == 0 else 1 if stage <= 2 else 2',
     'n = 0 if stage == 0 else 3', 'exec')], {'settings':settings})
guard_record = types.FunctionType(adaptive.guard_record.__code__,
    dict(vars(adaptive), wide=types.SimpleNamespace(guard_record=triple.guard_record), settings=settings),
    adaptive.guard_record.__name__)
evaluate_task, WORKER_JOIN = triple.derived(adaptive.evaluate_task, [
    ('(0,1,3)', '(0,3)', 'eval'), ('n = 1 if stage == 1 else 2', 'n = 3', 'exec')],
    dict(vars(adaptive), settings=settings, native_make_record=make_record, guard_record=guard_record))
derived_emitter = types.FunctionType(adaptive.derived_emitter.__code__,
                                    dict(vars(adaptive), PREFIX=PREFIX), 'derived_emitter')
COMMON_EMIT, EMITTER_JOIN = derived_emitter()
emit = types.FunctionType(adaptive.emit.__code__,
                         dict(vars(adaptive), PREFIX=PREFIX, COMMON_EMIT=COMMON_EMIT), 'emit')
replay_adapter = types.FunctionType(adaptive.replay_adapter.__code__,
                                   dict(vars(adaptive), PREFIX=PREFIX, emit=emit), 'replay_adapter')
OLD_SCOPE = ('Finite independent single/pair outer integration at fixed wide boxes with accepted '
             'three-momentum values held; no uniform or independent-grade convergence, global '
             'exceptional coverage, physical tails, Abel limit, scattering or poles.')
SCOPE = ('Finite independent triple outer quadrature at fixed cutoffs3/4 with accepted '
         'inner32/32, source/profile256/512 and single/pair values held. No uniform or '
         'independent-grade/global exceptional coverage, infinite tails, Abel limit, scattering or poles.')
finish, FINISH_JOIN = triple.derived(adaptive.finish, [
    (repr('wide-adaptive.pickle'), repr('wide-three-adaptive.pickle'), 'eval'),
    (repr(OLD_SCOPE), repr(SCOPE), 'eval')],
    dict(vars(adaptive), settings=settings, guard_record=guard_record, replay_adapter=replay_adapter,
         EMITTER_JOIN=EMITTER_JOIN))
METHOD_JOINS = {'record':MAKE_JOIN, 'worker':WORKER_JOIN, 'finish':FINISH_JOIN,
                'tripleGuard':triple.GUARD_JOIN, 'conditional':conditional.METHOD_JOINS,
                'emitter':EMITTER_JOIN}


def load(base):
    accepted, previous = domain.momentum.source.accepted(CHECKPOINT)
    require(accepted['status']=='PUBLISHED_ANNEX_VERIFIED', 'accepted triple refinement')
    for name,h in accepted['sourceFiles'].items():
        require(digest(ROOT/name)==h==digest(previous/'source'/name), 'accepted current/frozen source')
    for category in ('recordArtifacts','workerArtifacts'):
        for a in accepted[category]:
            require(digest(previous/a['path'])==a['sha256'], 'accepted saved artifact')
    bp = unpickle(previous/'bound-wide-three.pickle')
    packet = unpickle(previous/'wide-three.pickle')
    require(bp['sourceFiles']==packet['sourceFiles']==accepted['sourceFiles'] and
            bp['provenance']==packet['provenance']==accepted['provenance'], 'accepted packet provenance')
    _, source_base = domain.momentum.source.accepted(domain.momentum.source.SOURCE)
    r, dimensions = domain.momentum.source.native.source.restore_context(unpickle(source_base/'reduced-action.pickle'))
    dimensions.__dict__.update(packet['dimensionState'])
    refs = {}
    for task in TASKS:
        rec = next(v for v in packet['result']['records'] if (v['test'],v['domain'],v['stage'])==(*task,3))
        ref = dict(bp['refs'][task], action=rec['action'], groups=rec['groups'],
                   contributions=[{'index':v['index'],'terms':v['values']} for v in rec['terms']],
                   layoutSettings=rec['settings'], setting=rec['settings'][3])
        rebuilt = position.assemble(r,bp['domains'][task[1]],ref,rec['groups'])
        require(np.array_equal(rebuilt['action'],rec['action']) and
                all(not np.any(v['differences']) for v in rebuilt['terms']), 'accepted native action')
        s = rec['settings']
        require(s[1]['outerOrder']==s[2]['outerOrder']==432 and tuple(s[2]['innerOrders'])==(64,),
                'accepted single/pair rules')
        require(s[3]['outerOrder']==216 and tuple(s[3]['innerOrders'])==(32,32), 'accepted triple rule')
        require(all((v['sourceNodes'],v['profileNodes'],v['sourceBound'],v['profileBound'],v['regulator'])
                    ==(256,512,48.,14.,.2) for v in s.values()), 'fixed source/profile/domain/regulator')
        refs[task] = ref
    correction = json.loads(REPAIR.read_text())
    require(correction['status']=='ACCEPTED_ADDITIVE_CONTRACT_REPAIR', 'effective pole correction')
    require(digest(AUTHORITY)=='2b6018faf0d5ca772bf9c980174d8ebe7545127a3bf43ba01d1b3c0b6a5be423',
            'approved pole authority identity')
    pins = dict(accepted['sourceFiles'])
    for path in (CHECKPOINT,PLAN,Path(__file__).resolve(),AUTHORITY,REPAIR):
        pins[str(path.relative_to(ROOT))] = digest(path)
    for name in pins:
        dest=base/'source'/name; dest.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(ROOT/name,dest)
    for name in ('bound-wide-three.pickle','wide-three.pickle'):
        shutil.copyfile(previous/name,base/('accepted-'+name))
        require(digest(base/('accepted-'+name))==accepted['artifacts'][name]['sha256'], 'accepted byte copy')
    provenance = dict(accepted['provenance'], WIDE_THREE_ADAPTIVE_INPUT_SHA256=digest(CHECKPOINT),
        WIDE_THREE_ADAPTIVE_INSTRUMENT_SHA256=digest(Path(__file__)), WIDE_THREE_ADAPTIVE_PLAN_SHA256=digest(PLAN),
        NONLINEAR_POLE_V2_AUTHORITY_SHA256=digest(AUTHORITY), NONLINEAR_POLE_REPAIR_SHA256=digest(REPAIR),
        POLE_DEPENDENCY_DISPOSITION=DISPOSITION)
    data = {'r':r,'domains':bp['domains'],'joins':bp['joins'],'refs':refs,'sourceFiles':pins,
            'provenance':provenance,'previous':previous,'smoke':False}
    atomic_pickle(base/'bound-wide-three-adaptive.pickle',
                  {k:v for k,v in data.items() if k!='r'} | {'dimensionState':dict(vars(dimensions))})
    save(base/'preflight.json', {'sourceFiles':pins,'provenance':provenance,'tasks':TASKS,
         'methodJoins':METHOD_JOINS,'acceptedRunDirectory':str(previous),'nativeEngineUnchanged':True})
    return data


class NativeOneOuter(engine.BoundedSourceFourierQuadrature.ThreeMomentum):
    """Retain native recursion, replacing only its first rule by one unit node."""
    def __init__(self,*args,outer_value,**kwargs):
        super().__init__(*args,**kwargs); self.fixed_outer=float(outer_value); self.first=True; self.captured=[]

    def rule(self,lower,upper,order,centers=(),width=None):
        if self.first:
            self.first=False
            require(not centers and lower<self.fixed_outer<upper, 'native first outer rule')
            return np.array([self.fixed_outer]),np.ones(1),[lower,upper]
        return super().rule(lower,upper,order,centers,width)

    def batches(self,*args):
        self.first=True
        for points,weights in super().batches(*args):
            self.captured.append((points.copy(),weights.copy()))
            yield points,weights


class CapturedConditional(conditional.ConditionalMomentum):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs); self.captured=[]

    def batches(self,*args):
        for points,weights in super().batches(*args):
            self.captured.append((points.copy(),weights.copy()))
            yield points,weights


def focused_checks(base,data):
    artifacts=[]; timings=[]
    class PrefixComplete(Exception): pass
    def run(worker,ti,variables,s,bound,width,bounded):
        saved=[]; started=time.monotonic()
        def stop(state):
            if state['batchCount']==64:
                saved.append(state); raise PrefixComplete()
        try:
            group=worker.group(ti,variables,s,bound['pairs'],width,bound['positions'],stop if bounded else None)
        except PrefixComplete:
            group=saved[0]
        if bounded: require(len(saved)==1 and group['nodeCount']==16384, 'complete conditional prefix')
        return group,time.monotonic()-started
    for task in TASKS:
        ti,di=task; bound=data['domains'][di]; ref=data['refs'][task]; r=data['r']; s=ref['layoutSettings'][3]
        variables=next(g['variables'] for g in ref['groups'] if len(g['variables'])==3)
        width=float(bound['abel']['width'].subs(r.regulator,s['regulator']))
        probe=engine.BoundedSourceFourierQuadrature.ThreeMomentum(bound['rows'],bound['sources'],r)
        nodes,_,_=probe.rule(-s['momentumBound'],s['momentumBound'],s['outerOrder'])
        for index,k in enumerate((nodes[0],nodes[len(nodes)//2],0.)):
            bounded=index<2; setting=s if bounded else settings(ref,3,True)[3]
            a=CapturedConditional(bound['rows'],bound['sources'],r,outer_value=k)
            b=NativeOneOuter(bound['rows'],bound['sources'],r,outer_value=k)
            actual,ta=run(a,ti,variables,setting,bound,width,bounded)
            expected,tb=run(b,ti,variables,setting,bound,width,bounded)
            points,weights=[np.concatenate([v[i] for v in a.captured]) for i in (0,1)]
            ep,ew=[np.concatenate([v[i] for v in b.captured]) for i in (0,1)]
            proof={'task':task,'index':index,'bounded':bounded,'outerValue':float(k),'setting':setting,
                'actual':actual,'expected':expected,'actualPoints':points,'expectedPoints':ep,
                'actualWeights':weights,'expectedWeights':ew,
                'valueResidual':actual['values']-expected['values'],
                'mutationResidual':actual['measureMutationValues']-expected['measureMutationValues'],
                'massResidual':actual['quadratureMass']-expected['quadratureMass'],
                'conditionalUnits':[tuple(x-y for x,y in zip(bound['rows'][i]['unit'],bound['momentumUnit']))
                                    for i in actual['rowIndices']],
                'pointUnit':bound['momentumUnit'],'weightUnit':tuple(2*x for x in bound['momentumUnit']),
                'actualTelemetry':parallel.cache_state(a),'expectedTelemetry':parallel.cache_state(b),
                'wallSeconds':{'conditional':ta,'native':tb}}
            path=base/f'conditional-{ti}-{di}-{index}.pickle';atomic_pickle(path,proof);artifacts.append(artifact(base,path))
            require(np.array_equal(points,ep) and np.array_equal(weights,ew), 'native conditional nodes/weights')
            require(not np.any(proof['valueResidual']) and not np.any(proof['mutationResidual'])
                    and proof['massResidual']==0, 'native conditional values/mutations/mass')
            require(actual['rowIndices']==expected['rowIndices']==[v['index'] for v in bound['rows'] if len(v['limits'])==3]
                    and actual['variables']==expected['variables']==variables and actual['test']==expected['test']==ti,
                    'conditional row/field/limit identity')
            require(not np.any((abs(actual['values'])>1e-9)&
                    (abs(actual['measureMutationValues']-actual['values'])<1e-12)), 'actual measure sensitivity')
            if not bounded:
                require(abs(actual['volumeResidual'])<1e-10*(1+abs(actual['boxVolume'])) and
                        actual['boxVolume']==(2*setting['momentumBound'])**2, 'conditional finite area')
            timings.append({'task':task,'index':index,'bounded':bounded,'nodes':actual['nodeCount'],
                            'conditionalSeconds':ta,'nativeSeconds':tb})
    return {'conditionalArtifacts':artifacts,'conditionalPrefixes':8,'completeCoarseConditionals':4,
            'timings':timings,'methodJoins':METHOD_JOINS}


coarse_checks = adaptive.coarse_checks


def dispatch(base,data):
    old=parallel.evaluate_task; parallel.evaluate_task=evaluate_task
    try: return parallel.dispatch(base,data,TASKS)
    finally: parallel.evaluate_task=old


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--run-directory',type=Path,required=True)
    parser.add_argument('--preflight',action='store_true'); args=parser.parse_args()
    base=args.run_directory.resolve(); base.relative_to(STORE); base.mkdir(parents=True,exist_ok=False)
    started=time.monotonic(); data=load(base)
    focused=focused_checks(base,data) if args.preflight else None
    if focused is not None: save(base/'focused-before-workers.json',focused)
    data['smoke']=args.preflight; workers=dispatch(base,data)
    if args.preflight:
        focused['coarseArtifacts']=coarse_checks(base,data,workers); save(base/'focused-after-workers.json',focused)
    summary=finish(base,data,workers,started,focused); summary['methodJoins']=METHOD_JOINS
    save(base/'checks.json',summary); print(json.dumps(summary,indent=2))


if __name__=='__main__': main()
