#!/usr/bin/env python3
"""A selected finite boundary/regulator comparison with explicit modal phases."""
import argparse
import ast
import copy
import inspect
import json
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import textwrap
import time

import numpy as np
import sympy as sp
import S11c_d_finite_scattering_resolution as study
import S11c_d_position_domain as position

f=study.f
PLAN=f.M/'S11c_d_finite_scattering_domain_plan.md'
CHECKPOINT=f.M/'S11c_d_finite_scattering_resolution_checkpoint.json'
CASES=(('matching_basis',48.,0.2),('boundary_source',64.,0.2),('regulator',64.,0.1))


def tree(function):
    return ast.parse(textwrap.dedent(inspect.getsource(function))).body[0]


def compiled(node,namespace):
    module=ast.fix_missing_locations(ast.Module(body=[node],type_ignores=[]))
    scope=dict(namespace);exec(compile(module,__file__,'exec'),scope)
    return scope[node.name]


def adapters():
    old=tree(f.construct);new=copy.deepcopy(old);new.name='variable_construct'
    f.require(ast.unparse(new.body[1])=='interval = 48.0','original interval assignment')
    interval_assignment=new.body.pop(1)
    new.args.args.extend([ast.arg(arg='interval'),ast.arg(arg='regulator')])
    new.args.defaults.extend([ast.Constant(48.),ast.Constant(0.2)])
    class Regulator(ast.NodeTransformer):
        count=0
        def visit_Constant(self,node):
            if type(node.value) is float and node.value==0.2:
                self.count+=1;return ast.copy_location(ast.Name(id='regulator',ctx=ast.Load()),node)
            return node
    replace=Regulator()
    # Replace only body bindings, leaving default arguments literal.
    new.body=[replace.visit(n) for n in new.body]
    f.require(replace.count==2,'both actual regulator bindings')
    restored=copy.deepcopy(new);restored.name=old.name;restored.args.args=restored.args.args[:-2];restored.args.defaults=restored.args.defaults[:-2]
    class Reverse(ast.NodeTransformer):
        def visit_Name(self,node):
            return ast.copy_location(ast.Constant(0.2),node) if node.id=='regulator' else node
    restored.body=[Reverse().visit(n) for n in restored.body];restored.body.insert(1,interval_assignment)
    f.require(ast.dump(restored)==ast.dump(old),'whole finite constructor reverse AST join')
    construct=compiled(new,vars(f))

    old_rebind=tree(position.rebind);new_rebind=copy.deepcopy(old_rebind);new_rebind.name='variable_rebind'
    class Cutoffs(ast.NodeTransformer):
        def __init__(self,mapping):self.mapping=mapping;self.counts={n:0 for n in mapping}
        def visit_Constant(self,node):
            if type(node.value) is int and node.value in self.mapping:
                self.counts[node.value]+=1;return ast.copy_location(ast.Constant(self.mapping[node.value]),node)
            return node
    cut=Cutoffs({32:48,10:14});new_rebind=cut.visit(new_rebind)
    f.require(cut.counts=={32:3,10:3},'accepted source/profile cutoff sites')
    restored=Cutoffs({48:32,14:10}).visit(copy.deepcopy(new_rebind));restored.name=old_rebind.name
    f.require(ast.dump(restored)==ast.dump(old_rebind),'whole cutoff adapter reverse AST join')
    rebind=compiled(new_rebind,vars(position))

    old_inspect=tree(study.inspect);new_inspect=copy.deepcopy(old_inspect);new_inspect.name='variable_inspect'
    calls=[n for n in ast.walk(new_inspect) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr=='polynomial_basis']
    f.require(len(calls)==1 and ast.unparse(calls[0].args[1])=='48.0','common-grid basis interval site')
    calls[0].args[1]=ast.parse("system['settings']['sourceBound']",mode='eval').body
    restored=copy.deepcopy(new_inspect);restored.name=old_inspect.name
    next(n for n in ast.walk(restored) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr=='polynomial_basis').args[1]=ast.Constant(48.)
    f.require(ast.dump(restored)==ast.dump(old_inspect),'whole inspector reverse AST join')
    inspect_case=compiled(new_inspect,vars(study))
    return construct,rebind,inspect_case,{'constructorReverseJoin':True,'regulatorSites':replace.count,
        'cutoffReverseJoin':True,'cutoffSites':cut.counts,'inspectorReverseJoin':True}


def phases(system,solution):
    ends={'LEFT':float(system['nodes'][0]),'RIGHT':float(system['nodes'][-1])}
    incoming=[np.exp(1j*v['k']*ends[end]) for end,c in system['channels'].items() for v in c['incoming']]
    outgoing=[np.exp(-1j*v['k']*ends[end]) for end,v in solution['outgoingLabels']]
    f.require(all(v['k'].imag==0 for _,v in solution['outgoingLabels']) and
        all(v['k'].imag==0 for c in system['channels'].values() for v in c['incoming']),'open-mode phase frame')
    return np.asarray(incoming),np.asarray(outgoing)


def observable(directory,inspect_case):
    record=inspect_case(directory)
    system=f.unpickle(directory/'finite-system.pickle');solution=f.unpickle(directory/'finite-solution.pickle')
    inc,out=phases(system,solution)
    origin=out[:,None]*record['scattering']*inc[None,:]
    inverse=np.diag(1/out);current=inverse.conj().T@record['current']@inverse
    residual=np.diag(origin.conj().T@current@origin).real/record['incomingFlux']-record['totalCurrentRatio']
    f.require(np.max(abs(residual))<1e-12,'actual rephased full-current contraction')
    record.update(originScattering=origin,originFields=record['fields']*inc[None,None,:],
        originCurrent=current,originCurrentResidual=residual,incomingPhase=inc,outgoingPhase=out,
        interval=system['settings']['sourceBound'],regulator=system['settings']['regulator'])
    return record


def compare(previous,current):
    a=dict(previous,scattering=previous['originScattering'],fields=previous['originFields'])
    b=dict(current,scattering=current['originScattering'],fields=current['originFields'])
    difference=study.compare(a,b)
    difference['boundaryAnchoredModalDifferences']=difference.pop('modalDifferences')
    difference['phaseFrame']='Open amplitudes and incident-normalized common-grid fields at origin; evanescent modes remain at their actual boundaries.'
    return difference


def case(base,index,seconds):
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    def timeout(*_):raise TimeoutError('selected finite-domain budget; preserve all completed work')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(seconds);started=time.monotonic()
    name,interval,regulator=CASES[index];f.require(interval>0 and regulator>0,'finite positive domain/regulator')
    pilot=json.loads(study.CHECKPOINT.read_text());data=f.load(base,Path(pilot['bindingReuse']['directory']))
    construct,rebind,_,proofs=adapters();r,bound,local,jets,channels,pins,operands=data
    rebound,cutoff_proof=rebind(r,bound,(interval,14.))
    if interval==48.:f.require(rebound['rows']==bound['rows'] and rebound['cutoffBindings']==bound['cutoffBindings'],'unchanged baseline cutoff data')
    f.require(rebound['profiles']==bound['profiles'] and rebound['sources']==bound['sources'],'same profile/source operands')
    bounded_sources={i:sp.Integral(v['originalBoundAmplitude']*bound['sources'][(0,i)]['boundCharacter'],
        (r.zp,-sp.Rational(str(interval)),sp.Rational(str(interval)))) for i,v in jets.items()}
    for path in (Path(__file__),Path(position.__file__),Path(study.__file__),PLAN,CHECKPOINT):
        key=str(path.resolve().relative_to(f.ROOT));pins[key]=f.digest(path)
        target=base/'source'/key;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(path.read_bytes())
    data=(r,rebound,local,jets,channels,pins,operands)
    f.atomic_pickle(base/'domain-binding.pickle',{'bound':rebound,'cutoffProof':cutoff_proof,
        'genericBoundedSources':bounded_sources,'adapterProofs':proofs})
    inputs=json.loads((base/'inputs.json').read_text());inputs.update(sourceFiles=pins,
        numericalInterval=interval,numericalRegulator=regulator,adapterProofs=proofs,
        numericalUnits={'interval':list(map(str,f.engine.PHYSICAL_METADATA.dimensions.measure(r.zp))),
                        'regulator':list(map(str,f.engine.PHYSICAL_METADATA.dimensions.measure(r.regulator)))})
    f.save(base/'inputs.json',inputs)
    result,setting=construct(base,data,129,16,4,256,512,interval=interval,regulator=regulator)
    f.require(all(f.digest(f.ROOT/n)==h for n,h in pins.items()),'case sources unchanged')
    f.require(all(f.digest(Path(n))==h for n,h in operands.items()),'case inputs unchanged')
    summary={'runDirectory':str(base),'case':name,'unknowns':645,'rank':result['rank'],
        'balancedCondition':result['balancedCondition'],'settings':setting,
        'maximumEquationResidual':float(np.max(abs(result['equationResidual']))),
        'maximumScaledEquationResidual':float(np.max(abs(result['scaledEquationResidual']))),
        'outgoingFluxRatio':result['outgoingFluxRatio'].tolist(),
        'maximumMatrixActionResidual':max(float(np.max(abs(g['actionResidual']))) for g in result['groups']),
        'newMomentumNodes':sum(g['nodes'] for g in result['groups']),
        'wallSeconds':time.monotonic()-started,'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'sourceFiles':pins,'operandHashes':operands,
        'artifacts':{p.name:{'bytes':p.stat().st_size,'sha256':f.digest(p)} for p in base.glob('*.pickle')},
        'scope':'Finite domain/regulator response with approximate modal boundaries; continuum expansion and small-signal resolution pending.'}
    f.save(base/'checks.json',summary);signal.alarm(0);print(json.dumps(summary,indent=2))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True)
    parser.add_argument('--case',type=int,choices=range(3));parser.add_argument('--seconds',type=int,default=900)
    args=parser.parse_args();base=args.run_directory.resolve();base.relative_to(f.STORE);base.mkdir(parents=True,exist_ok=False)
    if args.case is not None:return case(base,args.case,args.seconds)
    started=time.monotonic();checkpoint=json.loads(CHECKPOINT.read_text())
    f.require(checkpoint['status']=='VALIDATED_FINITE_RESPONSE_RESOLUTION','accepted finite resolution study')
    baseline=Path(checkpoint['checks']['records'][-1]['checks']['runDirectory'])
    f.require(f.digest(Path(checkpoint['runDirectory'])/'checks.json')==checkpoint['checksSha256'],'accepted resolution payload')
    _,_,inspect_case,proofs=adapters();previous=observable(baseline,inspect_case)
    pins={str(p.resolve().relative_to(f.ROOT)):f.digest(p) for p in
        (Path(__file__),Path(f.__file__),Path(study.__file__),Path(position.__file__),PLAN,CHECKPOINT,f.ACCEPTANCE)}
    for name in pins:
        target=base/'source'/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes((f.ROOT/name).read_bytes())
    f.save(base/'inputs.json',{'sourceFiles':pins,'cases':CASES,'adapterProofs':proofs,
        'reportingTargets':{'resolvedRelative':0.01,'amplitudeAbsolute':1e-4,'normalizedCurrentAbsolute':1e-6}})
    environment=os.environ.copy();environment.update({k:'1' for k in
        ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')})
    records=[]
    for index,(name,_,_) in enumerate(CASES):
        remaining=int(900-(time.monotonic()-started));f.require(remaining>0,'remaining comparison budget')
        child=base/name;child.mkdir();directory=child/'complete'
        command=[sys.executable,'-u',str(Path(__file__)),'--run-directory',str(directory),'--case',str(index),'--seconds',str(remaining)]
        f.save(child/'invocation.json',{'command':command});then=time.monotonic()
        with (child/'stdout.txt').open('xb') as out,(child/'stderr.txt').open('xb') as err:
            run=subprocess.run(command,cwd=f.ROOT,env=environment,stdin=subprocess.DEVNULL,stdout=out,stderr=err)
        f.save(child/'outcome.json',{'exitCode':run.returncode,'wallSeconds':time.monotonic()-then,'stderrBytes':(child/'stderr.txt').stat().st_size})
        f.require(run.returncode==0 and (child/'stderr.txt').stat().st_size==0,('finite-domain case completed',name))
        current=observable(directory,inspect_case)
        f.require(current['checks']==json.loads((child/'stdout.txt').read_text()),'case stdout identity')
        difference=compare(previous,current)
        f.atomic_pickle(child/'comparison.pickle',{'previous':previous,'current':current,'difference':difference})
        records.append({'name':name,'checks':current['checks'],'difference':difference['summary'],
            'boundaryResidualMaxima':current['boundaryResidualMaxima'],
            'independentSolveCoefficientDifference':current['independentSolveCoefficientDifference'],
            'comparisonSha256':f.digest(child/'comparison.pickle')})
        f.save(base/'inventory.json',records);previous=current
    f.require(all(f.digest(f.ROOT/n)==h for n,h in pins.items()),'study sources unchanged')
    summary={'sourceFiles':pins,'records':records,'wallSeconds':time.monotonic()-started,
        'scope':'Selected finite boundary/source and regulator sensitivity; no transparent-boundary, infinite-limit or continuum-expansion claim.'}
    f.save(base/'checks.json',summary);print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
