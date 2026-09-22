#!/usr/bin/env python3
"""One new complete numerical row, with saved source nodes and bounded GK checks."""
import argparse
import ast
import gc
import inspect
import json
from pathlib import Path
import pickle
import resource
import signal
import time

import numpy as np
import scipy.integrate as integrate
import sympy as sp
import S11c_d_remaining_case_frequency_end_continuation_inputs as io
import S11c_d_remaining_case_frequency_end_pilot as storage

M,REPO=io.M,io.REPO
F=REPO/'_scratch/s11c/s11c-remaining-case-frequency-20260921'
READY=F/'row-inputs';READY_SHA='6b0c5edbb2567a9cf81a95c518f80098073b3ad8b3867debc6ceaa1289d4cc39'
INVENTORY_SHA='2d025399874fe772b8e3199916ccd97978db0ca33f8ef1cb0bbc73a2f1340a4f'
CP=M/'S11c_d_remaining_case_frequency_source_coefficients_checkpoint.json'
CP_SHA='4330b4ed945b3cf6f64f74e08e850ee33676fdce163945a899e832c441370782'
CASE='LAB_HELD__RHOBR_CONSTANT';ROW=28
require,same=io.require,io.same


def evaluate(expression, environment, memo=None):
    """Literal new numerical expression evaluation, without symbolic compilation."""
    memo={} if memo is None else memo
    if expression in memo:return memo[expression]
    if expression in environment:value=environment[expression]
    elif expression is sp.I:value=1j
    elif expression is sp.pi:value=np.pi
    elif isinstance(expression,sp.Rational):value=int(expression.p)/int(expression.q)
    elif isinstance(expression,sp.Float):value=float(expression)
    elif isinstance(expression,sp.Add):value=sum(evaluate(v,environment,memo) for v in expression.args)
    elif isinstance(expression,sp.Mul):
        value=1
        for v in expression.args:value=value*evaluate(v,environment,memo)
    elif isinstance(expression,sp.Pow) and isinstance(expression.exp,sp.Rational):
        exponent=int(expression.exp) if isinstance(expression.exp,sp.Integer) else int(expression.exp.p)/int(expression.exp.q)
        operand=np.asarray(evaluate(expression.base,environment,memo),dtype=complex)
        if not isinstance(expression.exp,sp.Integer):
            require(np.all(operand.imag==0) and np.all(operand.real>0),'noninteger numerical power requires positive-real saved branch')
        value=operand**exponent
    elif expression.func in (sp.exp,sp.tanh,sp.sin,sp.cos):
        function={sp.exp:np.exp,sp.tanh:np.tanh,sp.sin:np.sin,sp.cos:np.cos}[expression.func]
        value=function(evaluate(expression.args[0],environment,memo))
    else:raise TypeError(('unsupported numerical row node',type(expression),repr(expression)))
    memo[expression]=value
    return value


def source_matrix(coefficients,nodes,weights,bound,size):
    """New full source action via simultaneous Chebyshev value/jet recurrence."""
    degree=len(coefficients)-1
    require(degree<=2,'declared at-most-second-order source pilot')
    x=nodes/bound
    previous=[np.zeros_like(nodes) for _ in range(degree+1)]
    current=[np.ones_like(nodes)]+[np.zeros_like(nodes) for _ in range(degree)]
    value=np.empty((len(nodes),size),complex)
    for column in range(size):
        value[:,column]=weights*sum(a*b for a,b in zip(coefficients,current))
        if column==0:
            following=[x]+([np.ones_like(nodes)/bound] if degree else [])+([np.zeros_like(nodes)] if degree==2 else [])
        else:
            following=[2*x*current[0]-previous[0]]
            for n in range(1,degree+1):following.append(2*n/bound*current[n-1]+2*x*current[n]-previous[n])
        previous,current=current,following
    return value


def independent_source(coefficients,nodes,weights,bound,size):
    theta=np.arccos(nodes/bound);n=np.arange(size)[None,:];angle=theta[:,None]*n
    functions=[np.cos(angle)]
    if len(coefficients)>1:functions.append(n*np.sin(angle)/(bound*np.sin(theta)[:,None]))
    if len(coefficients)>2:functions.append(n*(np.sin(angle)*np.cos(theta)[:,None]-n*np.cos(angle)*np.sin(theta)[:,None])/(bound**2*np.sin(theta)[:,None]**3))
    return weights[:,None]*sum(a[:,None]*b for a,b in zip(coefficients,functions))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-directory',required=True,type=Path)
    base=ap.parse_args().run_directory.resolve();base.relative_to(REPO/'_scratch/s11c');base.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));signal.alarm(900);start=time.monotonic()
    reader,journal=io.Reader(),storage.Journal(base)
    cp=reader.json(CP,CP_SHA);require(cp['status']=='ACCEPTED_CASE_FREQUENCY_SOURCE_COEFFICIENTS','accepted actual source coefficients')
    ready=reader.json(READY/'complete/checks.json',READY_SHA);require(ready['newScientificCalls']==0,'saved row inspection')
    inventory=reader.json(READY/'completed-input-artifact-inventory.json',INVENTORY_SHA)
    def metadata(name):return reader.json(READY/'complete'/name,inventory[name]['sha256'])
    def load(record):return reader.packet(record.get('logical',record.get('path')),record['sha256'])
    def address(record):
        value=load(record['packet'])
        for k in record['keys']:value=value[k]
        return value
    selected=metadata('cases/'+CASE+'/row-'+str(ROW)+'.json')
    require(selected['status']=='UNMATCHED_FULL_ROW_INPUT' and not selected['savedCompleteMatches'] and selected['layoutDimension']==1,'one genuinely missing whole 1D row')
    old=load(selected['sourceInputs']['packet']);scalars=load(selected['scalarInputs']);system=load(selected['basis'])
    physical=reader.json(selected['ownPhysicalRoutes']['packet']['logical'],selected['ownPhysicalRoutes']['packet']['sha256'])
    common=load(physical['context']);context=common['contextPair'][0]
    require(same(*common['contextPair']) and same(*common['basisPair']),'full accepted own/shared physical context and units')
    for v in physical.values():
        if isinstance(v,dict) and 'logical' in v:reader.retain(v['logical'],v['sha256'])
    row=old['bound']['rows'][ROW];require(row['index']==ROW and len(row['factors'])==1,'actual single-factor row pilot')
    factor=row['factors'][0];si=factor['sourceIndex'];jet=address(selected['jets'][str(si)]);source=old['bound']['sources'][0,si]
    basis_candidates=next(v for v in metadata('cases/'+CASE+'/source-basis-inputs.json') if v['sourceIndex']==si)
    require(not basis_candidates['savedCoefficientBasisCandidates'] and same(basis_candidates['coefficientRoute'],selected['jets'][str(si)]),'no full saved source-basis candidate across inspected original outputs')
    coefficient=scalars['actual']['factor',ROW,0];variable=row['limits'][0][0];settings=old['settings'];positions=system['nodes'];size=len(positions)
    require(same(settings,system['settings']) and size==129 and jet['probe'].args==(context['zp'],),'actual full source/trial/coordinate inputs')
    require(same(jet['originalBoundAmplitude'],scalars['actual']['source',si]) and same(jet['amplitudeUnit'],source['amplitudeUnit']) and same(jet['integralUnit'],source['integralUnit']),'exact own new source amplitude and units')
    source_meta=metadata('saved-prepared-basis/'+CASE+'.json');prepared=load(source_meta['packet']);saved=prepared['original']
    nodes=saved['source_nodes'];weights=saved['source_weights'];source_bound=settings['sourceBound']
    require(saved['size']==size and context['settings']==settings and same(source_meta['sourceSettings'],json.loads(json.dumps(settings))),'whole native saved source-rule caller settings')
    for v in saved['jet_data'].values():
        require(not same((jet['probe'],jet['coefficients'],jet['amplitudeUnit'],jet['integralUnit']),
                         (v['probe'],v['coefficients'],v['amplitudeUnit'],v['integralUnit'])),'no completed saved own source-basis input is reconstructed')
    require(nodes.shape==weights.shape==(4*settings['sourceNodes'],) and np.all(weights>0) and np.all(abs(nodes)<source_bound),'complete retained source nodes and weights')
    source_route={'packet':source_meta['packet'],'keys':['original','source_nodes']};weight_route={'packet':source_meta['packet'],'keys':['original','source_weights']}
    del prepared,saved;gc.collect()
    scope={'case':CASE,'rowIndex':ROW,'sourceIndex':si,'frequency':str(scalars['frequency']),'physicalInput':selected,'sourceNodes':source_route,'sourceWeights':weight_route,
        'trialSize':size,'positions':selected['basis'],'method':'new numerical tree evaluator + new full source-action recurrence + adaptive Gauss-Kronrod',
        'firstRule':'gk21','secondRule':'gk15','epsabs':[1e-8,5e-9],'epsrel':[1e-6,5e-7],'intervalLimit':256,'maximumDistinctMomentumPoints':4096,
        'comparisonScaledTolerance':1e-6,'sourceBasisTolerance':2e-10,'remainingSecondsRule':'second route only if3* measured first cost +60 fits remaining900s',
        'scope':'One new full129x129 integral-row pilot, not a scattering/pole/domain result or a replacement for the approved native16/4/4 quadrature results. Source256/profile512 settings and physical finite domains remain recorded; no profile integration enabled in this single-row pilot.'}
    journal.json('pilot-scope.json',scope)
    journal.write('row-input.pickle',{'row':row,'coefficient':coefficient,'source':source,'jet':jet,'context':context,'settings':settings,'physicalRoutes':physical,'scalarInputs':selected['scalarInputs'],'positions':positions,'fieldUnits':old['fieldUnits'],'equationUnits':old['equationUnits']})
    require(not coefficient.has(sp.Integral) and len(row['limits'])==1,'single-row pilot requires no uninspected old profile transform')
    require(coefficient.free_symbols<={variable,context['z'],context['regulator']} and source['frequency'].free_symbols<={variable} and all(a.free_symbols<={context['zp']} for a in jet['coefficients']),'all actual coefficient and source-frequency variables covered')
    # Pin whole native caller and the explicitly new numerical implementation.
    native_meta=metadata('native-callers.json')
    for record in native_meta.values():reader.retain(record['file']['logical'],record['file']['sha256'])
    paths=(Path(__file__).resolve(),M/'S11c_d_remaining_case_frequency_row_pilot_plan.md',Path(io.__file__),Path(storage.__file__),M.parent/'directives/S11c_d_EXPLORATORY_ACCEPTANCE.md',Path(inspect.getsourcefile(integrate.quad_vec)))
    for p in paths:reader.retain(p)
    journal.json('native-and-new-numerical-callers.json',{'native':native_meta,'new':{v.__name__:inspect.getsource(v) for v in (evaluate,source_matrix,independent_source)},'quadVecSource':reader.retain(inspect.getsourcefile(integrate.quad_vec)),'nativeBasisOrRuleConstructorCalled':False})
    def forbidden(*args,**kwargs):raise RuntimeError('completed native science disabled in new row pilot')
    for name in ('diff','lambdify','cancel','expand','factor','solve','gcd','resultant','integrate'):setattr(sp,name,forbidden)
    io.native.f.source_jets=io.native.f.polynomial_basis=io.native.f.BasisMomentum.prepare_basis=forbidden
    io.native.Pair.__init__=io.native.maps=io.native.continue_pair=forbidden
    # Actual new coefficient values and new full source-action call.
    journal.write('source-action/input.pickle',{'jet':jet,'nodeRoute':source_route,'weightRoute':weight_route,'nodes':nodes,'weights':weights,'bound':source_bound,'size':size,'context':context})
    values=[]
    for n,a in enumerate(jet['coefficients']):
        prefix='source-action/coefficients/'+str(n)
        journal.write(prefix+'/input.pickle',{'expression':a,'variable':context['zp'],'nodes':source_route,'jet':selected['jets'][str(si)]})
        with np.errstate(over='raise',invalid='raise',divide='raise',under='ignore'):
            v=np.broadcast_to(np.asarray(evaluate(a,{context['zp']:nodes}),complex),nodes.shape)
        output=journal.write(prefix+'/value.pickle',v);journal.json(prefix+'/completed.json',{'value':output})
        values.append(v)
    journal.write('source-action/coefficient-values.pickle',values)
    matrix=source_matrix(values,nodes,weights,source_bound,size)
    source_output=journal.write('source-action/value.pickle',matrix);journal.json('source-action/completed.json',{'value':source_output})
    journal.write('source-action/independent-input.pickle',{'coefficientValues':journal.artifacts['source-action/coefficient-values.pickle'],'sourceInput':journal.artifacts['source-action/input.pickle'],'method':'closed trigonometric Chebyshev value/derivatives'})
    other=independent_source(values,nodes,weights,source_bound,size);journal.write('source-action/independent-value.pickle',other)
    source_spread=float(np.max(abs(matrix-other))/(1+np.max(abs(matrix))))
    journal.json('source-action/comparison.json',{'maximumScaledDifference':source_spread,'tolerance':2e-10})
    require(np.isfinite(matrix).all() and source_spread<2e-10,'new source-action recurrence/trigonometric comparison')
    del other
    cache={};factor_cache={};number=0;new_evaluations=0;reused_evaluations=0
    def integrand(momentum):
        nonlocal number,new_evaluations,reused_evaluations
        key=float(momentum).hex()
        if key in cache:
            reused_evaluations+=1;record=cache[key]
            with Path(record['path']).open('rb') as stream:return pickle.load(stream)
        require(new_evaluations<4096,'bounded unique numerical pilot evaluations')
        prefix='integrand/'+str(number);number+=1
        journal.write(prefix+'/input.pickle',{'momentum':float(momentum),'variable':variable,'positions':selected['basis'],'rowInput':journal.artifacts['row-input.pickle'],'sourceAction':source_output})
        environment={variable:float(momentum),context['z']:positions,context['regulator']:settings['regulator']}
        with np.errstate(over='raise',invalid='raise',divide='raise',under='ignore'):
            c=np.broadcast_to(np.asarray(evaluate(coefficient,environment),complex),positions.shape)
            journal.write(prefix+'/coefficient-value.pickle',c)
            frequency=complex(evaluate(source['frequency'],{variable:float(momentum)}))
            journal.write(prefix+'/frequency-value.pickle',frequency)
            require(abs(frequency.imag)<1e-14,'actual source Fourier frequency remains real')
            source_value=np.exp(-1j*frequency.real*nodes)@matrix
            value=c[:,None]*source_value[None,:]
        factor_cache[key]=journal.write(prefix+'/factors.pickle',{'coefficient':c,'sourceFrequency':frequency,'sourceFourier':source_value})
        output=journal.write(prefix+'/value.pickle',value)
        journal.json(prefix+'/completed.json',{'input':journal.artifacts[prefix+'/input.pickle'],'value':output})
        require(value.shape==(129,129) and np.isfinite(value).all(),'finite full own row integrand')
        cache[key]=output;new_evaluations+=1
        return value
    # New actual sign/measure controls at a declared nonzero probe momentum.
    control_momentum=.731
    control_value=integrand(control_momentum)
    with Path(factor_cache[control_momentum.hex()]['path']).open('rb') as stream:control_factors=pickle.load(stream)
    c,freq=control_factors['coefficient'],control_factors['sourceFrequency']
    journal.write('controls/input.pickle',{'momentum':control_momentum,'actualIntegrand':cache[control_momentum.hex()],'actualFactors':factor_cache[control_momentum.hex()],'sourceAction':source_output,'wrongPhaseSign':1,'weightScale':1.001})
    wrong=c[:,None]*(np.exp(1j*freq.real*nodes)@matrix)[None,:];changed=control_value*1.001
    journal.write('controls/value.pickle',{'wrongPhaseSign':wrong,'changedMeasure':changed})
    responses={'sourcePhaseSignResponse':float(np.max(abs(wrong-control_value))),'measureResponse':float(np.max(abs(changed-control_value)))}
    journal.json('controls/checks.json',responses);require(all(v>1e-12 for v in responses.values()),'actual source phase and integration measure controls respond')
    results=[];costs=[];limit=float(settings['momentumBound'])
    for index,(rule,ea,er) in enumerate((('gk21',1e-8,1e-6),('gk15',5e-9,5e-7))):
        if index:
            cost={'firstRouteSeconds':costs[0],'remainingSeconds':900-(time.monotonic()-start),'requiredReserveSeconds':3*costs[0]+60}
            journal.json('measured-pilot-cost.json',cost);require(cost['requiredReserveSeconds']<cost['remainingSeconds'],'measured second-route cost fits whole-job budget')
        prefix='quadrature/'+str(index);journal.write(prefix+'/input.pickle',{'interval':(-limit,limit),'rule':rule,'epsabs':ea,'epsrel':er,'norm':'max','workers':1,'limit':256,'cacheSize':32*1024**2,'rowInput':journal.artifacts['row-input.pickle'],'sourceAction':source_output})
        tick=time.monotonic();value,error,info=integrate.quad_vec(integrand,-limit,limit,epsabs=ea,epsrel=er,norm='max',quadrature=rule,workers=1,cache_size=32*1024**2,limit=256,full_output=True)
        cost=time.monotonic()-tick;costs.append(cost)
        output=journal.write(prefix+'/value.pickle',{'row':value,'errorEstimate':error,'success':info.success,'status':info.status,'message':info.message,'evaluations':info.neval,'intervals':info.intervals,'intervalValues':info.integrals,'intervalErrors':info.errors,'wallSeconds':cost})
        journal.json(prefix+'/completed.json',{'input':journal.artifacts[prefix+'/input.pickle'],'value':output,'newIntegrandEvaluations':new_evaluations,'reusedIntegrandEvaluations':reused_evaluations})
        require(info.success and np.isfinite(value).all(),'bounded adaptive row convergence')
        results.append(value)
    difference=results[1]-results[0];absolute_difference=float(np.max(abs(difference)));row_norm=float(np.max(abs(results[1])));comparison=absolute_difference/(1+row_norm)
    journal.write('row-comparison.pickle',{'difference':difference,'maximumAbsoluteDifference':absolute_difference,'rowMaxNorm':row_norm,'maximumScaledDifference':comparison,'first':journal.artifacts['quadrature/0/value.pickle'],'second':journal.artifacts['quadrature/1/value.pickle']})
    require(comparison<1e-6,'selected independent GK-rule row comparison')
    reader.postcheck();journal.json('inputs.json',{'acceptedSourceCoefficients':reader.retain(CP,CP_SHA),'savedRowInspection':reader.retain(READY/'complete/checks.json',READY_SHA),'consumedRoutes':reader.routes})
    checks={'status':'COMPLETED_BOUNDED_NEW_FREQUENCY_ROW_PILOT','case':CASE,'rowIndex':ROW,'sourceIndex':si,'rowShape':[129,129],
        'sourceBasisScaledDifference':source_spread,'rowScaledDifference':comparison,'rowMaxNorm':row_norm,'rowAbsoluteDifference':absolute_difference,'newIntegrandEvaluations':new_evaluations,'reusedIntegrandEvaluations':reused_evaluations,'quadratureSeconds':costs,'controlResponses':responses,
        'allConsumedHashesUnchanged':True,'artifacts':dict(journal.artifacts),'wallSeconds':time.monotonic()-start,
        'scope':'One genuinely missing full1D row at1-0.01i via two bounded adaptive rules. Observed numerical spread only, not a scattering error bound, pole/domain or full four-case response. No old native source basis/rule/evaluator/row/map recomputed.'}
    journal.json('checks.json',checks);signal.alarm(0);print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
