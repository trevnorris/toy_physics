#!/usr/bin/env python3
"""Focused finite momentum rule, source and nested-profile instrument checks."""
import argparse
import ast
import json
from pathlib import Path
import shutil
import time
import numpy as np
import sympy as sp
import S11c_d_momentum_action_check as m


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--directory',type=Path,required=True);args=parser.parse_args()
    base=args.directory.resolve();base.relative_to(m.STORE);base.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    paths=(m.engine.HERE,Path(m.__file__).resolve(),Path(__file__).resolve())
    pins={str(p.relative_to(m.ROOT)):m.digest(p) for p in paths}
    for p in paths:
        target=base/'source'/p.relative_to(m.ROOT);target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
    accepted,sb=m.source.accepted(m.SOURCE_CHECKPOINT);dc,db=m.source.accepted(m.DOMAIN_CHECKPOINT)
    nc,nb=m.source.accepted(m.source.NUMERICAL);sc,rb=m.source.accepted(m.source.SOURCE)
    packet=m.unpickle(sb/'bound-sources.pickle');bound=packet['result']
    domains=m.unpickle(db/'domains.pickle')['result'];r,dimensions=m.source.native.source.restore_context(m.unpickle(rb/'reduced-action.pickle'))
    contractor=m.engine.BoundedSourceFourierQuadrature.FiniteMomentum((),{},r)
    old=ast.parse((sb/'source'/m.engine.HERE.relative_to(m.ROOT)).read_text());new=ast.parse(m.engine.HERE.read_text())
    cls=next(n for n in new.body if getattr(n,'name',None)=='BoundedSourceFourierQuadrature')
    additions=[n for n in cls.body if getattr(n,'name',None)=='FiniteMomentum']
    if len(additions)!=1:raise ValueError('momentum helper census')
    cls.body.remove(additions[0])
    if ast.dump(old)!=ast.dump(new):raise ValueError('whole accepted engine AST join')
    setting=dict(m.unpickle(nb/'actions.pickle')['result']['results'][0]['settings'],kind='legacy')
    abel=domains['abel'];width=float(abel['width'].subs(r.regulator,setting['regulator']))
    density=sp.lambdify((abel['momentum'],abel['center'],r.regulator),abel['density'],'numpy')
    primitive=sp.lambdify((abel['momentum'],abel['center'],r.regulator),abel['primitive'],'numpy')
    rules=[]
    for center in (0.,0.13,1.99):
        x,w,points=contractor.rule(-2.,2.,16,(center,),width)
        actual=np.dot(w,density(x,center,setting['regulator']));expected=primitive(2.,center,setting['regulator'])-primitive(-2.,center,setting['regulator'])
        unsplit_x,unsplit_w,_=contractor.rule(-2.,2.,16)
        ordinary=np.dot(unsplit_w,density(unsplit_x,center,setting['regulator']))
        rules.append({'center':center,'panels':points,'nodes':len(x),'massResidual':float(w.sum()-4.),
            'secondMomentResidual':float(np.dot(w,x*x)-16/3),'abelResidual':[float((actual-expected).real),float((actual-expected).imag)],
            'omittedPanelsDifference':float(abs(actual-ordinary)),'weightMutation':float(abs(np.dot(w*1.001,density(x,center,setting['regulator']))-actual))})
    environment={v:np.array([-1.7,-.3,.1,1.2])*(1-i/5) for i,v in enumerate(bound['momenta'])}
    environment[r.regulator]=np.full(4,setting['regulator'])
    profiles=[]
    compiler=m.engine.BoundedActionQuadrature({r.xi:(-setting['profileBound'],setting['profileBound'],setting['profileNodes'])})
    for index,p in enumerate(domains['profiles']):
        finite=sp.Integral(p['bound'].function,(r.xi,-sp.Rational(str(setting['profileBound'])),sp.Rational(str(setting['profileBound']))))
        actual=contractor.profile_value(finite,environment,setting,{})
        expected=compiler.integrate(p['bound'],environment)
        profiles.append({'index':index,'actual':actual,'expected':expected,'residual':actual-expected})
    sources=[];nodes,weights=m.engine.BoundedSourceFourierQuadrature.rule((-setting['sourceBound'],setting['sourceBound']),setting['sourceNodes'])
    for record in bound['records']:
        value=contractor.source_value(dict(record,testWidth=bound['testWidthsMomenta'][record['test']][0],profileWidth=bound['profileWidth']),environment,setting,{})
        function=sp.lambdify((r.zp,*bound['momenta']),record['boundSource'],'numpy',cse=True,docstring_limit=0)
        literal=np.asarray([np.dot(weights,np.broadcast_to(np.asarray(function(nodes,*(environment[v][j] for v in bound['momenta'])),dtype=complex),nodes.shape)) for j in range(4)])
        sources.append({'test':record['test'],'sourceIndex':record['sourceIndex'],'actual':value,'expected':literal,'residual':value-literal})
    m.atomic_pickle(base/'operands.pickle',{'rules':rules,'profiles':profiles,'sources':sources,'environment':environment,'setting':setting,'width':width})
    summary={'nativeEngineAstJoin':True,'sourceFiles':pins,'sourceCheckpointSha256':m.digest(m.SOURCE_CHECKPOINT),'domainCheckpointSha256':m.digest(m.DOMAIN_CHECKPOINT),
        'engineSha256':m.digest(m.engine.HERE),'checkerSha256':m.digest(Path(m.__file__)),'preflightSha256':m.digest(Path(__file__)),
        'profileCount':len(profiles),'boundSourceCount':len(sources),'rules':rules,
        'profileMaxResidual':max(float(np.max(abs(v['residual']))) for v in profiles),
        'sourceMaxScaledResidual':max(float(np.max(abs(v['residual'])/(1+abs(v['expected'])))) for v in sources),
        'operandsSha256':m.digest(base/'operands.pickle'),'wallSeconds':time.monotonic()-start}
    m.save(base/'checks.json',summary)
    if pins!={str(p.relative_to(m.ROOT)):m.digest(p) for p in paths}:raise ValueError('focused source changed during checks')
    if (summary['profileMaxResidual']>1e-10 or summary['sourceMaxScaledResidual']>1e-10 or
        any(max(map(abs,v['abelResidual']))>1e-10 or abs(v['massResidual'])>1e-12 or abs(v['secondMomentResidual'])>1e-12 or
            min(v['omittedPanelsDifference'],v['weightMutation'])<=1e-10 for v in rules)):
        raise ValueError('finite momentum focused instrument guard')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
