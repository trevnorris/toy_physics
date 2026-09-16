#!/usr/bin/env python3
"""Check conditional native operands and run an explicitly underresolved adaptive smoke test."""
import argparse,itertools,json,time
from pathlib import Path
import numpy as np
import S11c_d_three_momentum_adaptive as a


class PrefixComplete(Exception):pass


def prefix(worker,ti,variables,setting,bound,width):
    states=[]
    def checkpoint(s):
        if s['batchCount']==64:states.append(s);raise PrefixComplete()
    try:worker.group(ti,variables,setting,bound['pairs'],width,bound['positions'],checkpoint)
    except PrefixComplete:pass
    a.require(len(states)==1,'conditional prefix complete')
    return states[0]


def relative(left,right):
    return float(np.max(abs(left-right)))/max(float(np.max(abs(left))),float(np.max(abs(right))),1e-300)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--directory',type=Path,required=True);options=parser.parse_args()
    base=options.directory.resolve();base.relative_to(a.STORE);base.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    data=a.load(base);r,bound,rows,variables,finest,provenance,pins=data
    accepted,previous=a.native.momentum.source.accepted(a.CHECKPOINT)
    evidence=[];bounds=[]
    for ti in range(2):
        setting=finest[ti]['setting'];width=float(bound['abel']['width'].subs(r.regulator,setting['regulator']))
        original=a.engine.BoundedSourceFourierQuadrature.ThreeMomentum(bound['rows'],bound['sources'],r)
        outer_nodes,outer_weights,_=original.rule(-setting['momentumBound'],setting['momentumBound'],setting['outerOrder'])
        conditioned=a.ConditionalMomentum(bound['rows'],bound['sources'],r,outer_value=outer_nodes[0])
        node_diff=weight_diff=0.
        for x,y in zip(itertools.islice(original.batches(variables,setting,bound['pairs'],width),64),itertools.islice(conditioned.batches(variables,setting,bound['pairs'],width),64),strict=True):
            node_diff=max(node_diff,float(np.max(abs(x[0]-y[0]))));weight_diff=max(weight_diff,float(np.max(abs(x[1]-y[1]*outer_weights[0]))))
        path=previous/f'worker-{ti}-2/partials/00000064.pickle'
        artifact=next(v for v in accepted['workerArtifacts'] if v['path']==str(path.relative_to(previous)))
        a.require(a.digest(path)==artifact['sha256'],'saved original prefix hash')
        saved=a.unpickle(path);fresh=prefix(conditioned,ti,variables,setting,bound,width)
        a.require(saved['setting']==setting and saved['rowIndices']==fresh['rowIndices'] and saved['variables']==fresh['variables'],'prefix source/layout join')
        value=fresh['values']*outer_weights[0];mutation=fresh['measureMutationValues']*outer_weights[0]
        record={'test':ti,'sourceArtifact':artifact,'outerValue':outer_nodes[0],'outerWeight':outer_weights[0],
            'native':saved,'conditioned':fresh,'valueResidual':value-saved['values'],'mutationResidual':mutation-saved['measureMutationValues'],
            'relativeValueResidual':relative(value,saved['values']),'relativeMutationResidual':relative(mutation,saved['measureMutationValues']),
            'massResidual':fresh['quadratureMass']*outer_weights[0]-saved['quadratureMass'],'nodeResidual':node_diff,'weightResidual':weight_diff,
            'measureMutationResponse':float(np.max(abs(mutation-value)))}
        evidence.append(record);a.atomic_pickle(base/f'prefix-{ti}.pickle',record)
        a.require(node_diff==0 and weight_diff<1e-18 and record['relativeValueResidual']<1e-10 and record['relativeMutationResidual']<1e-10 and abs(record['massResidual'])<1e-12,'native conditioned prefix')
        a.require(record['measureMutationResponse']>0,'prefix actual measure response')
        # Cheap complete grouping comparison changes only numerical rule orders.
        setting=dict(setting,outerOrder=2,innerOrders=(1,1),sourceNodes=32,profileNodes=32)
        original=a.engine.BoundedSourceFourierQuadrature.ThreeMomentum(bound['rows'],bound['sources'],r)
        group=original.group(ti,variables,setting,bound['pairs'],width,bound['positions'])
        nodes,weights,_=original.rule(-setting['momentumBound'],setting['momentumBound'],setting['outerOrder'])
        points=[]
        for k in nodes:
            worker=a.ConditionalMomentum(bound['rows'],bound['sources'],r,outer_value=k)
            points.append(worker.group(ti,variables,setting,bound['pairs'],width,bound['positions']))
        combined=dict(group,values=sum(g['values']*w for g,w in zip(points,weights)),measureMutationValues=sum(g['measureMutationValues']*w for g,w in zip(points,weights)))
        combined['quadratureMass']=sum(g['quadratureMass']*w for g,w in zip(points,weights))
        direct=a.native.full_action(r,bound,finest[ti],group);assembled=a.native.full_action(r,bound,finest[ti],combined)
        full={'test':ti,'setting':setting,'native':group,'conditionedPoints':points,'combined':combined,'relativeValueResidual':relative(group['values'],combined['values']),
            'relativeMutationResidual':relative(group['measureMutationValues'],combined['measureMutationValues']),
            'massResidual':group['quadratureMass']-combined['quadratureMass'],'actionResidual':direct['action']-assembled['action'],
            'termsResiduals':[np.asarray(x['values'])-np.asarray(y['values']) for x,y in zip(direct['terms'],assembled['terms'],strict=True)]}
        a.atomic_pickle(base/f'full-grouping-{ti}.pickle',full);bounds.append(full)
        a.require(full['relativeValueResidual']<1e-10 and full['relativeMutationResidual']<1e-10 and abs(full['massResidual'])<1e-12 and np.max(abs(full['actionResidual']))<1e-12,'full native conditional grouping')
    coarse=list(data);coarse[4]={ti:dict(v,adaptiveOverrides={'innerOrders':(1,1),'sourceNodes':16,'profileNodes':16,'adaptiveTolerance':1e-3}) for ti,v in finest.items()};coarse=tuple(coarse)
    a.save(base/'focused-before-adaptive.json',{'prefixes':[{'test':v['test'],'relativeValueResidual':v['relativeValueResidual'],'relativeMutationResidual':v['relativeMutationResidual'],'nodeResidual':v['nodeResidual'],'weightResidual':v['weightResidual'],'massResidual':v['massResidual']} for v in evidence],
        'fullGrouping':[{'test':v['test'],'relativeValueResidual':v['relativeValueResidual'],'relativeMutationResidual':v['relativeMutationResidual'],'massResidual':v['massResidual'],'maxActionResidual':float(np.max(abs(v['actionResidual'])))} for v in bounds],
        'scope':'Production-order prefixes and complete underresolved grouping checks; subsequent adaptive smoke is also underresolved.'})
    packets=a.dispatch(base,coarse);checks=a.finish(base,coarse,packets)
    checks.update(focused=json.loads((base/'focused-before-adaptive.json').read_text()),wallSeconds=time.monotonic()-started,
        scope='Instrument preflight: accepted production-order conditional prefixes, underresolved complete grouping and adaptive smoke only. No new physical quadrature result or convergence claim.')
    checks['focusedArtifacts']={p.name:a.digest(p) for p in base.glob('prefix-*.pickle')}
    checks['focusedArtifacts'].update({p.name:a.digest(p) for p in base.glob('full-grouping-*.pickle')})
    a.save(base/'checks.json',checks);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
