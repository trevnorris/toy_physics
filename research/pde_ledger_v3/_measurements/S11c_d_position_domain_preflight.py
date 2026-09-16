#!/usr/bin/env python3
"""Bounded checks of position-cutoff binding, complete actions and four workers."""
import argparse,copy,json,time
from pathlib import Path
import numpy as np
import sympy as sp
import S11c_d_position_domain as domain
from S11c_d_position_domain import engine,require,digest,save,atomic_pickle,unpickle


def check_transforms(base,data):
    r=data['r'];momenta=tuple(r.normal_map[g[2]] for g in r.momentum_groups)
    inventory=[];summaries=[]
    for index in (1,2):
        bound=data['domains'][index];setting=domain.settings(data['finest'][0],domain.CHOICES[index])[3]
        worker=engine.BoundedSourceFourierQuadrature.ThreeMomentum(bound['rows'],bound['sources'],r)
        if index==1: # Source cutoff 48 is identical in the two new domain cases.
            for (ti,si),record in sorted(bound['sources'].items()):
                env={k:record['assignments'][:,j] for j,k in enumerate(momenta)}
                value=worker.source_value(record,env,setting,{})
                refined=worker.source_value(record,env,dict(setting,sourceNodes=384),{})
                cut=setting['sourceBound'];width=float(record['testWidth']);ell=record['profileWidth']
                points=sorted({-cut,cut,0.,-width,width,-ell,ell})
                x,w=engine.BoundedSourceFourierQuadrature.rule(points,256)
                original=sp.lambdify((r.zp,*momenta),record['boundSource'],'numpy',cse=True)
                direct=np.asarray([np.dot(np.broadcast_to(np.asarray(original(x,*assignment),dtype=complex),x.shape),w) for assignment in record['assignments']])
                # Mutation acts on the actual source integration weights.
                q=engine.BoundedSourceFourierQuadrature(r.zp,record['boundAmplitude'])
                mutated=q.gauss(record['frequencies'],points,256,weight_scale=1.001)
                residual=value-direct;change=refined-value;response=mutated-value
                result={'test':ti,'sourceIndex':si,'setting':setting,'frequencies':record['frequencies'],'assignments':record['assignments'],
                    'original':record['boundSource'],'values':value,'direct':direct,'refined':refined,'sourceResidual':residual,'orderChange':change,'mutationValues':mutated,'mutationResidual':response,
                    'unit':record['integralUnit']}
                path=base/'transforms'/f'source-{ti}-{si:02}.pickle';path.parent.mkdir(exist_ok=True);atomic_pickle(path,result);inventory.append(domain.artifact(base,path));save(base/'transform-inventory.json',inventory)
                require(np.max(abs(residual)/(1+abs(direct)))<1e-9 and np.max(abs(change)/(1+abs(refined)))<1e-9,'finite source/order residual')
                require(float(np.max(abs(response)))>0,'source weight response')
                summaries.append({'kind':'source','test':ti,'sourceIndex':si,'scaledOriginalResidual':float(np.max(abs(residual)/(1+abs(direct)))),'scaledOrderChange':float(np.max(abs(change)/(1+abs(refined))))})
        for j,proof in enumerate(data['joins'][index]['profiles']):
            integral=proof['new'];original=next(v for v in bound['profiles'] if v['bound'].function==integral.function)
            info=engine.BoundedSourceFourierQuadrature.affine_range(original['transfer'],{k:(-sp.Integer(2),sp.Integer(2)) for k in momenta})
            fractions=np.linspace(0,1,65);env={k:(-2+4*fractions if info['coefficients'][k]>0 else 2-4*fractions if info['coefficients'][k]<0 else np.zeros(65)) for k in momenta}
            value=worker.profile_value(integral,env,setting,{})
            refined=worker.profile_value(integral,env,dict(setting,profileNodes=384),{})
            x,w=engine.BoundedSourceFourierQuadrature.rule([-setting['profileBound'],0,setting['profileBound']],256)
            fn=sp.lambdify((r.xi,*momenta),integral.function,'numpy',cse=True)
            direct=np.asarray([np.dot(np.asarray(fn(x,*[env[k][i] for k in momenta]),complex),w) for i in range(65)])
            mutated=np.asarray([np.dot(np.asarray(fn(x,*[env[k][i] for k in momenta]),complex),w*1.001) for i in range(65)])
            result={'index':index,'profile':j,'integral':integral,'transfer':original['transfer'],'range':info,'environment':env,'setting':setting,'values':value,'direct':direct,'refined':refined,
                'sourceResidual':value-direct,'orderChange':refined-value,'mutationValues':mutated,'mutationResidual':mutated-value,'unit':proof['unit']}
            path=base/'transforms'/f'profile-{index}-{j}.pickle';path.parent.mkdir(exist_ok=True);atomic_pickle(path,result);inventory.append(domain.artifact(base,path));save(base/'transform-inventory.json',inventory)
            require(np.max(abs(value-direct))<1e-9 and np.max(abs(refined-value))<1e-9,'profile/order residual')
            require(np.max(abs(mutated-value))>0,'profile weight response')
            summaries.append({'kind':'profile','index':index,'profile':j,'originalResidual':float(np.max(abs(value-direct))),'orderChange':float(np.max(abs(refined-value)))})
    return inventory,summaries


def first_batch(worker,ti,variables,setting,bound,r):
    saved=[]
    def stop(state):saved.append(state);raise domain.parallel.PrefixComplete()
    try:worker.group(ti,variables,setting,bound['pairs'],float(bound['abel']['width'].subs(r.regulator,setting['regulator'])),bound['positions'],stop)
    except domain.parallel.PrefixComplete:pass
    require(len(saved)==1,'bounded prefix completion');return saved[0]


def check_prefixes(base,data):
    r=data['r'];original=data['bound'];baseline=data['domains'][0];records=[]
    for ti in range(2):
        for group in data['finest'][ti]['groups']:
            variables=tuple(group['variables']);n=len(variables);setting=domain.settings(data['finest'][ti],domain.CHOICES[0])[n]
            a=first_batch(engine.BoundedSourceFourierQuadrature.ThreeMomentum(original['rows'],original['sources'],r),ti,variables,setting,original,r)
            b=first_batch(engine.BoundedSourceFourierQuadrature.ThreeMomentum(baseline['rows'],baseline['sources'],r),ti,variables,setting,baseline,r)
            proof={'test':ti,'layout':n,'original':a,'rebound':b,'valueResidual':a['values']-b['values'],'mutationResidual':a['measureMutationValues']-b['measureMutationValues'],'massResidual':a['quadratureMass']-b['quadratureMass']}
            path=base/f'prefix-{ti}-{n}.pickle';atomic_pickle(path,proof)
            require(not np.any(proof['valueResidual']) and not np.any(proof['mutationResidual']) and proof['massResidual']==0,'baseline prefix identity')
            require(a['nodeCount']==b['nodeCount'] and a['rowIndices']==b['rowIndices'] and a['variables']==b['variables'],'prefix assignment/census')
            records.append(domain.artifact(base,path))
    return records


def check_coarse(base,data,workers):
    evidence=[];r=data['r']
    for (ti,index),w in sorted(workers.items()):
        bound=data['domains'][index];worker=engine.BoundedSourceFourierQuadrature.ThreeMomentum(bound['rows'],bound['sources'],r);groups=[]
        for g in w['groups']:
            setting=w['settings'][len(g['variables'])]
            fresh=worker.group(ti,g['variables'],setting,bound['pairs'],float(bound['abel']['width'].subs(r.regulator,setting['regulator'])),bound['positions'])
            require(np.array_equal(fresh['values'],g['values']) and np.array_equal(fresh['measureMutationValues'],g['measureMutationValues']) and fresh['quadratureMass']==g['quadratureMass'],'coarse serial/worker join');groups.append(fresh)
        # Independent literal assembly using every bound cell and all integral rows.
        lookup={idx:g['values'][j] for g in groups for j,idx in enumerate(g['rowIndices'])}
        action=np.zeros_like(w['action']);compiler=engine.BoundedActionQuadrature({})
        for cell in (c for c in bound['cells'] if c['test']==ti):
            for pi,z in enumerate(bound['positions']):
                env={r.z:z,r.regulator:.2}
                action[pi,cell['column'],cell['row']]=complex(compiler.evaluate(cell['local'],env))+sum(complex(compiler.evaluate(coef,env))*lookup[row][pi] for row,coef in cell['terms'])
        residual=action-w['action'];path=base/f'coarse-{ti}-{index}.pickle';atomic_pickle(path,{'directAction':action,'workerAction':w['action'],'residual':residual,'groups':groups})
        require(not np.any(residual),'all-cell direct/worker action join');evidence.append(domain.artifact(base,path))
    # A wrong profile limit must trigger the unchanged native domain guard.
    wrong=next(iter(data['domains'][0]['profileUnits']));worker=engine.BoundedSourceFourierQuadrature.ThreeMomentum(data['bound']['rows'],data['bound']['sources'],r)
    rejected=False
    try:worker.profile_value(wrong,{k:np.zeros(1) for k in (r.normal_map[g[2]] for g in r.momentum_groups)},domain.settings(data['finest'][0],domain.CHOICES[2])[3],{})
    except ValueError as e:rejected=str(e)=='nested profile cutoff mismatch'
    require(rejected,'incorrect profile-limit negative control');return evidence


def main():
    p=argparse.ArgumentParser();p.add_argument('--run-directory',type=Path,required=True);a=p.parse_args()
    base=a.run_directory.resolve();base.relative_to(domain.STORE);base.mkdir(parents=True,exist_ok=False);started=time.monotonic();data=domain.load(base)
    inventory,transforms=check_transforms(base,data);prefixes=check_prefixes(base,data)
    data['smoke']=True;save(base/'focused-before-smoke.json',{'transformChecks':transforms,'transformArtifacts':inventory,'prefixArtifacts':prefixes})
    workers=domain.dispatch(base,data);coarse=check_coarse(base,data,workers);summary=domain.finish(base,data,workers,started)
    summary.update(transformChecks=transforms,transformArtifacts=inventory,prefixArtifacts=prefixes,coarseArtifacts=coarse,profileLimitMutationRejected=True)
    for v in (*inventory,*prefixes,*coarse):require(digest(base/v['path'])==v['sha256'],'focused artifact identity')
    save(base/'checks.json',summary);print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
