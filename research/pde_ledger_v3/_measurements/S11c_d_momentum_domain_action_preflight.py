#!/usr/bin/env python3
"""Bounded native-worker, transform, rule, action and emission joins."""
import argparse,json,time
from pathlib import Path
import numpy as np
import S11c_d_momentum_domain_action as domain
import S11c_d_position_domain_preflight as previous
from S11c_d_momentum_domain_action import require,digest,save,atomic_pickle,unpickle,engine


def checks_before_workers(base,data):
    r=data['r'];old=unpickle(base/'accepted-position-domain.pickle');oldbound=unpickle(base/'accepted-bound-position-domains.pickle')['domains'][2];evidence=[];prefixes=[]
    # Actual production transform methods join the already accepted independent
    # Fourier records, so the chosen rule is tested at the newly bound ranges.
    prep=unpickle(base/'accepted-momentum-domain-preparation.pickle')
    for w in prep['result']['workers']:
        ti,oldindex=w['task'];index=oldindex+1;bound=data['domains'][index];worker=engine.BoundedSourceFourierQuadrature.ThreeMomentum(bound['rows'],bound['sources'],r);setting=domain.settings(data['finest'][ti],index)[3]
        momenta=data['joins'][index]['momenta']
        for record in w['records']:
            environment={k:record['assignments'][:,j] for j,k in enumerate(momenta)}
            if record['kind']=='source':
                v=worker.source_value(bound['sources'][(ti,record['sourceIndex'])],environment,setting,{});expected=record['gaussValues'][1];identity=record['sourceIndex']
            else:
                v=worker.profile_value(record['original'],environment,setting,{});expected=record['gaussValues'][-1];identity=record['profileIndex']
            residual=v-expected;adaptive=v-record['adaptiveValues'];scale=1+abs(record['adaptiveValues'])
            path=base/'transforms'/f'{ti}-{index}-{record["kind"]}-{identity:02}.pickle';path.parent.mkdir(exist_ok=True);atomic_pickle(path,{'record':record,'nativeValues':v,'ruleResidual':residual,'adaptiveResidual':adaptive,'setting':setting})
            require(np.max(abs(residual)/scale)<1e-9 and np.max(abs(adaptive)/scale)<1e-9,'selected native transform route/accepted adaptive join')
            evidence.append(dict(domain.artifact(base,path),kind=record['kind'],scaledRuleResidual=float(np.max(abs(residual)/scale)),scaledAdaptiveResidual=float(np.max(abs(adaptive)/scale))))
    for ti in range(2):
        for group in data['finest'][ti]['groups']:
            variables=group['variables'];n=len(variables);setting=domain.settings(data['finest'][ti],0)[n]
            a=previous.first_batch(engine.BoundedSourceFourierQuadrature.ThreeMomentum(oldbound['rows'],oldbound['sources'],r),ti,variables,setting,oldbound,r)
            b=previous.first_batch(engine.BoundedSourceFourierQuadrature.ThreeMomentum(data['bound']['rows'],data['bound']['sources'],r),ti,variables,setting,data['bound'],r)
            proof={'test':ti,'layout':n,'old':a,'bound':b,'valueResidual':a['values']-b['values'],'mutationResidual':a['measureMutationValues']-b['measureMutationValues'],'massResidual':a['quadratureMass']-b['quadratureMass']}
            path=base/f'prefix-{ti}-{n}.pickle';atomic_pickle(path,proof);require(not np.any(proof['valueResidual']) and not np.any(proof['mutationResidual']) and proof['massResidual']==0 and a['nodeCount']==b['nodeCount'] and a['rowIndices']==b['rowIndices'],'retained native prefix identity');prefixes.append(domain.artifact(base,path))
    # Explicit production-setting one-axis-change controls and native layout
    # cutoffs; these are instrument joins, not a convergence claim.
    for ti in range(2):
        for index in range(1,4):
            current=domain.settings(data['finest'][ti],index);prev=domain.settings(data['finest'][ti],index-1)
            for n,s in current.items():
                require({k for k in set(s)|set(prev[n]) if s.get(k)!=prev[n].get(k)}==({'profileNodes'} if index==1 else {'momentumBound'}),'isolated rule/domain parameter')
                require(all(tuple(map(float,limit[1:]))==(-s['momentumBound'],s['momentumBound']) for row in data['domains'][index]['rows'] if len(row['limits'])==n for limit in row['limits']),'actual native momentum limits/rule')
    return evidence,prefixes


def check_coarse(base,data,workers):
    artifacts=[];r=data['r']
    for (ti,index),w in sorted(workers.items()):
        bound=data['domains'][index];worker=engine.BoundedSourceFourierQuadrature.ThreeMomentum(bound['rows'],bound['sources'],r);groups=[]
        for g in w['groups']:
            s=w['settings'][len(g['variables'])];fresh=worker.group(ti,g['variables'],s,bound['pairs'],float(bound['abel']['width'].subs(r.regulator,s['regulator'])),bound['positions'])
            require(np.array_equal(fresh['values'],g['values']) and np.array_equal(fresh['measureMutationValues'],g['measureMutationValues']) and fresh['quadratureMass']==g['quadratureMass'] and fresh['nodeCount']==g['nodeCount'],'coarse native serial/worker identity');groups.append(fresh)
        values={i:g['values'][j] for g in groups for j,i in enumerate(g['rowIndices'])};action=np.zeros_like(w['action']);compiler=engine.BoundedActionQuadrature({})
        for cell in (v for v in bound['cells'] if v['test']==ti):
            for pi,z in enumerate(bound['positions']):
                env={r.z:z,r.regulator:.2};action[pi,cell['column'],cell['row']]=complex(compiler.evaluate(cell['local'],env))+sum(complex(compiler.evaluate(coef,env))*values[row][pi] for row,coef in cell['terms'])
        proof={'directAction':action,'workerAction':w['action'],'residual':action-w['action'],'groups':groups};path=base/f'coarse-{ti}-{index}.pickle';atomic_pickle(path,proof);require(not np.any(proof['residual']),'independent native full-action contraction');artifacts.append(domain.artifact(base,path))
    return artifacts


def main():
    p=argparse.ArgumentParser();p.add_argument('--run-directory',type=Path,required=True);a=p.parse_args();base=a.run_directory.resolve();base.relative_to(domain.STORE);base.mkdir(parents=True,exist_ok=False);started=time.monotonic();data=domain.load(base)
    transforms,prefixes=checks_before_workers(base,data);save(base/'focused-before-smoke.json',{'transformArtifacts':transforms,'prefixArtifacts':prefixes,'workerMethodJoin':domain.METHOD_JOIN})
    data['smoke']=True;workers=domain.dispatch(base,data);coarse=check_coarse(base,data,workers);summary=domain.finish(base,data,workers,started)
    summary.update(transformArtifacts=transforms,prefixArtifacts=prefixes,coarseArtifacts=coarse)
    for a in (*transforms,*prefixes,*coarse):require(digest(base/a['path'])==a['sha256'],'final focused evidence')
    save(base/'checks.json',summary);print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
