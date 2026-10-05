#!/usr/bin/env python3
"""Exact resume checks and serial/parallel comparisons on actual numerical operands."""
import argparse,ast,copy,inspect,json,shutil,textwrap,time
from pathlib import Path
import numpy as np
import S11c_d_parallel_momentum as p


def same_group(a,b):
    keys=('values','measureMutationValues','quadratureMass','nodeCount','batchCount','rowIndices','variables','test')
    residuals={k:float(np.max(abs(np.asarray(a[k])-np.asarray(b[k])))) for k in ('values','measureMutationValues','quadratureMass','nodeCount','batchCount')}
    p.require(all(np.array_equal(a[k],b[k]) for k in keys),'serial/parallel group mismatch')
    return residuals


def reference_prefix(directory,data,task,stop=256):
    r,bound,rows,variables,finest,provenance,pins=data;ti,index=task
    setting=dict(finest[ti]['setting'],sourceNodes=256,profileNodes=128 if index==1 else 256)
    width=float(bound['abel']['width'].subs(r.regulator,setting['regulator']))
    worker=p.engine.BoundedSourceFourierQuadrature.ThreeMomentum(bound['rows'],bound['sources'],r);states={}
    def checkpoint(s):
        if s['batchCount'] in (64,stop):states[s['batchCount']]=s;p.atomic_pickle(directory/f'prefix-{s["batchCount"]}.pickle',s)
        if s['batchCount']==stop:raise p.PrefixComplete()
    started=time.monotonic()
    try:worker.group(ti,variables,setting,bound['pairs'],width,bound['positions'],checkpoint)
    except p.PrefixComplete:pass
    p.require(stop in states,'serial prefix completion')
    return {'states':states,'telemetry':p.cache_state(worker),'seconds':time.monotonic()-started}


def focused_constructor():
    original=ast.parse(textwrap.dedent(inspect.getsource(p.native.construct)));tree=copy.deepcopy(original)
    assignment=next(n for n in tree.body[0].body if isinstance(n,ast.Assign) and isinstance(n.targets[0],ast.Name) and n.targets[0].id=='choices')
    prior=copy.deepcopy(assignment.value);assignment.value=ast.parse('((144,24,24,128,128),(2,1,1,256,128),(2,1,1,256,256))',mode='eval').body
    proof=copy.deepcopy(tree);next(n for n in proof.body[0].body if isinstance(n,ast.Assign) and isinstance(n.targets[0],ast.Name) and n.targets[0].id=='choices').value=prior
    p.require(ast.dump(proof)==ast.dump(original),'focused quadrature-order-only AST join')
    ns=dict(vars(p.native));exec(compile(ast.fix_missing_locations(tree),'<focused-native-constructor>','exec'),ns)
    return ns['construct']


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--directory',type=Path,required=True);parser.add_argument('--mode',choices=('resume','parallel'),required=True);a=parser.parse_args()
    base=a.directory.resolve();base.relative_to(p.STORE);base.mkdir(parents=True,exist_ok=False);started=time.monotonic();data=p.prepare(base)
    if a.mode=='resume':
        d=base/'serial';d.mkdir();serial=reference_prefix(d,data,(0,1))
        d=base/'resumed';d.mkdir();resumed=p.evaluate_task(d,data,(0,1),max_batches=256,resume_state=serial['states'][64])
        residuals=same_group(serial['states'][256],resumed['group'])
        p.require(serial['telemetry']['sourceFrequencyCensus']==resumed['telemetry']['sourceFrequencyCensus'],'resume frequency census')
        bad=copy.deepcopy(serial['states'][64]);bad['nodeCount']+=1;d=base/'bad-count';d.mkdir();rejected=False
        try:p.evaluate_task(d,data,(0,1),max_batches=256,resume_state=bad)
        except ValueError as e:
            p.require(str(e)=='resume generated node count','unexpected mutation error');rejected=True
        p.require(rejected,'bad resume count accepted')
        result={'mode':a.mode,'residuals':residuals,'frequencyCensusExact':True,'badCountRejected':rejected,'resumedPrefixNodes':resumed['resumePrefixNodes'],'serialSeconds':serial['seconds'],'resumedSeconds':resumed['wallSeconds']}
    else:
        tasks=[(ti,index) for ti in range(2) for index in (1,2)];serial={};start=time.monotonic()
        for task in tasks:
            d=base/f'serial-prefix-{task[0]}-{task[1]}';d.mkdir();serial[task]=reference_prefix(d,data,task)
        serial_seconds=time.monotonic()-start;d=base/'parallel-prefix';d.mkdir();start=time.monotonic();parallel=p.dispatch(d,data,tasks,max_batches=256);parallel_seconds=time.monotonic()-start
        residuals={str(t):same_group(serial[t]['states'][256],parallel[t]['group']) for t in tasks}
        p.require(all(serial[t]['telemetry']['sourceFrequencyCensus']==parallel[t]['telemetry']['sourceFrequencyCensus'] for t in tasks),'parallel frequency census')
        coarse=list(data);coarse[4]={ti:dict(v,setting=dict(v['setting'],outerOrder=2,innerOrders=(1,1))) for ti,v in data[4].items()};coarse=tuple(coarse)
        constructor=focused_constructor();sd=base/'serial-full';pd=base/'parallel-full';sd.mkdir();pd.mkdir()
        for d in (sd,pd):shutil.copyfile(base/'accepted-bound-momentum.pickle',d/'accepted-bound-momentum.pickle')
        r,bound,rows,variables,finest,provenance,pins=coarse
        serial_full=constructor(sd,r,bound,rows,variables,finest)
        jobs=p.dispatch(pd,coarse,tasks);parallel_full=p.construct_from_workers(pd,coarse,jobs,constructor)
        full_residuals=[]
        for x,y in zip(serial_full['records'],parallel_full['records'],strict=True):
            same_group(x['group'],y['group']);p.require(np.array_equal(x['action'],y['action']),'parallel full action')
            for xx,yy in zip(x['terms'],y['terms'],strict=True):p.require(np.array_equal(xx['values'],yy['values']) and np.array_equal(xx['differences'],yy['differences']),'parallel full terms')
            full_residuals.append(float(np.max(abs(x['action']-y['action']))))
        p.engine.PAYLOAD_ENCODER=p.engine.PayloadEncoder()
        p.engine.EMISSION_LINES.clear();left,keys,paths=p.native.emit_and_replay(sd,serial_full,r,bound,rows,finest,provenance)
        p.engine.PAYLOAD_ENCODER=p.engine.PayloadEncoder()
        p.engine.EMISSION_LINES.clear();right,rkeys,rpaths=p.native.emit_and_replay(pd,parallel_full,r,bound,rows,finest,provenance)
        differences=[k for k in left if left[k]!=right.get(k)];p.save(base/'emission-differences.json',differences)
        p.require(left==right and keys==rkeys and paths==rpaths,'parallel emission/metadata payloads')
        result={'mode':a.mode,'prefixResiduals':residuals,'frequencyCensusExact':True,'serialPrefixWallSeconds':serial_seconds,'fourWorkerPrefixWallSeconds':parallel_seconds,'measuredPrefixSpeedup':serial_seconds/parallel_seconds,
            'fullActionResiduals':full_residuals,'fullRecordCount':len(serial_full['records']),'emissionTags':len(left),'writeKeys':len(keys),'metadataPaths':paths,'emissionDifferences':len(differences),'serialTranscriptSha256':p.digest(sd/'full.out'),'parallelTranscriptSha256':p.digest(pd/'full.out'),
            'workerPeakRssKiB':{str(t):v['peakRssKiB'] for t,v in parallel.items()},'focusedRuleOrders':{'outer':2,'inner':[1,1]},'scope':'Underresolved finite full quadratures test execution equivalence only; no physical convergence result.'}
    result.update(sourceFiles=data[-1],provenance=data[-2],methodJoins=p.METHOD_JOINS,wallSeconds=time.monotonic()-started)
    p.save(base/'checks.json',result);p.require(data[-1]=={n:p.digest(p.ROOT/n) for n in data[-1]},'focused sources changed');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
