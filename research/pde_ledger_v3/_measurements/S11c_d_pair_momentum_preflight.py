#!/usr/bin/env python3
"""Check the paired inner-at-outer evaluator and cached rules on actual rows."""
import argparse,json,time
from pathlib import Path
import numpy as np
import S11c_d_pair_momentum_check as s


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--directory',type=Path,required=True);args=parser.parse_args()
    base=args.directory.resolve();base.relative_to(s.STORE);base.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    r,bound,rows,variables,finest,provenance,pins=s.load(base)
    worker=s.engine.BoundedSourceFourierQuadrature.CachedMomentum(bound['rows'],bound['sources'],r)
    records=[]
    for ti in range(2):
        setting=finest[ti]['setting'];width=float(bound['abel']['width'].subs(r.regulator,setting['regulator']))
        value=worker.group(ti,variables,setting,bound['pairs'],width,bound['positions'])
        reference=next(g for g in finest[ti]['groups'] if g['variables']==variables)
        nodes,weights,_=worker.rule(-setting['momentumBound'],setting['momentumBound'],setting['outerOrder'])
        points=np.asarray([worker.inner_at_outer(ti,variables,setting,bound['pairs'],width,bound['positions'],x) for x in nodes])
        independent=np.sum(points*weights[:,None,None],axis=0)
        mutated=np.sum(points*(weights*1.001)[:,None,None],axis=0)
        mixed=s.full_action(r,bound,finest[ti],value)
        records.append({'test':ti,'group':value,'nodes':nodes,'weights':weights,'pointValues':points,'pointSum':independent,
            'pointSumResidual':independent-value['values'],'referenceResidual':value['values']-reference['values'],
            'weightMutation':mutated-independent,'mixedAction':mixed})
    rules=[]
    for (points,order),(nodes,weights) in worker.fixed_rules.items():
        actual=s.engine.BoundedSourceFourierQuadrature.rule(points,order)
        rules.append({'points':points,'order':order,'nodesEqual':np.array_equal(nodes,actual[0]),'weightsEqual':np.array_equal(weights,actual[1]),
            'readOnly':not nodes.flags.writeable and not weights.flags.writeable})
    s.atomic_pickle(base/'focused.pickle',{'records':records,'rules':rules})
    checks={'runDirectory':str(base),'provenance':provenance,'sourceFiles':pins,'nativeEngineAstJoin':True,'profileEvaluatorAstJoin':True,
        'rows':len(rows),'tests':2,'positions':len(bound['positions']),'pointEvaluations':sum(len(v['nodes']) for v in records),
        'maxPointSumResidual':max(float(np.max(abs(v['pointSumResidual']))) for v in records),
        'maxReferenceResidual':max(float(np.max(abs(v['referenceResidual']))) for v in records),
        'maxFullReferenceResidual':max(float(np.max(abs(v['mixedAction']['differenceFromReference']))) for v in records),
        'minTestMaxMutation':min(float(np.max(abs(v['weightMutation']))) for v in records),
        'ruleChecks':rules,'focusedPacketSha256':s.digest(base/'focused.pickle'),'wallSeconds':time.monotonic()-started}
    s.save(base/'checks.json',checks)
    if any(checks[k]>1e-10 for k in ('maxPointSumResidual','maxReferenceResidual','maxFullReferenceResidual')) or checks['minTestMaxMutation']<1e-10 or any(not v[k] for v in rules for k in ('nodesEqual','weightsEqual','readOnly')):
        raise ValueError('focused two-momentum guard')
    if pins!={n:s.digest(s.ROOT/n) for n in pins}:raise ValueError('focused source changed')
    print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
