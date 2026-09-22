#!/usr/bin/env python3
"""Finish only the incomplete GK15 row comparison using complete saved points."""
import ast
import copy
import inspect
import json
from pathlib import Path
import S11c_d_remaining_case_frequency_row_pilot_recover as recovery

prior=recovery.prior
M,F,require=prior.M,prior.F,prior.require
OLD=F/'row-pilot-recovery-01'
INVENTORY_SHA='911817bd3382466caeb85e195f747be2995478964304376457f1fd3f81717171'
RECOVERY_SHA='7d9d446165f6e3301dabd5dc225f1585930f5879a1feb65fa58e96ba24f82a40'
LOCALS=('coefficient','source','variable','positions','context','settings','source_output','selected','matrix','nodes','si','source_spread','cache','factor_cache','number','new_evaluations','reused_evaluations','responses','results','costs','route_input')


def adapted(text):
    module,previous_join=recovery.adapted(text)
    evaluate=copy.deepcopy(module.body[0])
    original=copy.deepcopy(module.body[1]);fn=copy.deepcopy(original)
    restoration=next(i for i,n in enumerate(fn.body) if isinstance(n,ast.Assign) and isinstance(n.value,ast.Call) and isinstance(n.value.func,ast.Name) and n.value.func.id=='restore')
    integrand_index=next(i for i,n in enumerate(fn.body) if isinstance(n,ast.FunctionDef) and n.name=='integrand')
    loop_index=next(i for i,n in enumerate(fn.body) if isinstance(n,ast.For) and 'gk21' in ast.unparse(n.iter))
    limit_index=loop_index-1
    require(ast.unparse(fn.body[limit_index]).startswith("limit = float(settings['momentumBound'])"),'native quadrature interval route')
    restore_node=ast.parse(','.join(LOCALS)+' = restore_remaining(base,reader,journal)').body[0]
    integrand=copy.deepcopy(fn.body[integrand_index])
    bound=next(n for n in ast.walk(integrand) if isinstance(n,ast.Compare) and ast.unparse(n.left)=='new_evaluations')
    require(ast.unparse(bound)=='new_evaluations < 4096','original declared point cap');bound.comparators=[ast.Constant(value=8192)]
    loop=copy.deepcopy(fn.body[loop_index]);old_loop=copy.deepcopy(loop)
    specs=loop.iter.args[0];require(len(specs.elts)==2,'exact two original rules')
    loop.iter=ast.List(elts=[ast.Tuple(elts=[ast.Constant(value=1),copy.deepcopy(specs.elts[1])],ctx=ast.Load())],ctx=ast.Load())
    for n in ast.walk(loop):
        if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and isinstance(n.func.value,ast.Name) and n.func.value.id=='journal':
            if n.func.attr=='json' and isinstance(n.args[0],ast.Constant) and n.args[0].value=='measured-pilot-cost.json':n.args[0]=ast.Constant(value='continuation-measured-cost.json')
            elif n.func.attr=='write' and isinstance(n.args[0],ast.BinOp) and isinstance(n.args[0].right,ast.Constant) and n.args[0].right.value=='/input.pickle':n.func=ast.Name(id='route_input',ctx=ast.Load())
    fn.body=fn.body[:restoration]+[restore_node,integrand,copy.deepcopy(fn.body[limit_index]),loop]+fn.body[loop_index+1:]
    extra=ast.parse("checks.update(priorCompletedIntegrands=4096,newIntegrandsThisContinuation=new_evaluations-4096,savedPointLoadsThisContinuation=reused_evaluations,priorUnfinishedRuleReuseCountNotReconstructed=True,completedGK21Reused=True,completedControlsReused=True,cumulativeDistinctPointBudget=8192,unfinishedGK15ControllerResumedFromSavedPointFunction=True)").body[0]
    output=next(i for i,n in enumerate(fn.body) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and ast.unparse(n.value.func)=='journal.json' and isinstance(n.value.args[0],ast.Constant) and n.value.args[0].value=='checks.json')
    fn.body.insert(output,extra)
    # Reverse every allowed change, including the complete unexecuted prior main.
    reverse=copy.deepcopy(fn);reverse.body.pop(output)
    i=next(i for i,n in enumerate(reverse.body) if isinstance(n,ast.FunctionDef) and n.name=='integrand')
    next(n for n in ast.walk(reverse.body[i]) if isinstance(n,ast.Compare) and ast.unparse(n.left)=='new_evaluations').comparators=[ast.Constant(value=4096)]
    require(ast.dump(reverse.body[i])==ast.dump(original.body[integrand_index]),'whole integrand unchanged except explicit measured point budget')
    actual_loop=next(n for n in reverse.body if isinstance(n,ast.For) and 'gk15' in ast.unparse(n.iter))
    actual_loop.iter=copy.deepcopy(old_loop.iter)
    for n in ast.walk(actual_loop):
        if isinstance(n,ast.Call):
            if isinstance(n.func,ast.Name) and n.func.id=='route_input':n.func=ast.Attribute(value=ast.Name(id='journal',ctx=ast.Load()),attr='write',ctx=ast.Load())
            elif ast.unparse(n.func)=='journal.json' and isinstance(n.args[0],ast.Constant) and n.args[0].value=='continuation-measured-cost.json':n.args[0]=ast.Constant(value='measured-pilot-cost.json')
    require(ast.dump(actual_loop)==ast.dump(old_loop),'whole original quadrature loop body and unchanged GK15 settings')
    new_loop_index=reverse.body.index(actual_loop)
    reverse.body[restoration:new_loop_index+1]=copy.deepcopy(original.body[restoration:loop_index+1])
    require(ast.dump(reverse)==ast.dump(original),'whole previous helper main reverse identity')
    return ast.fix_missing_locations(ast.Module(body=[evaluate,fn],type_ignores=[])),{'previous':previous_join,'wholePreviousMainReverseAST':True,'wholeIntegrandReverseAST':True,'wholeQuadratureLoopReverseAST':True,'restoredLocals':list(LOCALS),'priorCompletedPointCount':4096,'newPointAllowance':4096,'completedGK21Executed':False,'controlsExecuted':False}


def restore_remaining(base,reader,journal):
    inventory=reader.json(OLD/'completed-pilot-file-inventory.json',INVENTORY_SHA)
    failure=reader.json(OLD/'failed-point-cap-outcome.json')
    require(failure['actualExits']==[1,1,1] and failure['guardInterruption'] is None and failure['completedIntegrandValues']==4096 and not failure['secondRuleComplete'],'actual point-budget stop, completed prefix and unfinished comparison')
    # Keep a reference inventory; do not copy or relink thousands of old values.
    for name,record in inventory.items():
        route=reader.retain(record['logical'],record['sha256'])
        require(route==record,'unchanged actual logical/canonical point route')
        journal.artifacts[name]={'path':record['logical'],'sha256':record['sha256'],'bytes':record['bytes']}
    def packet(name):
        record=inventory[name];return reader.packet(record['logical'],record['sha256'])
    def metadata(name):
        record=inventory[name];return reader.json(record['logical'],record['sha256'])
    row=packet('row-input.pickle');source_input=packet('source-action/input.pickle');matrix=packet('source-action/value.pickle')
    scope=metadata('pilot-scope.json');selected=scope['physicalInput'];second_input=packet('quadrature/1/input.pickle')
    first=packet('quadrature/0/value.pickle');first_receipt=metadata('quadrature/0/completed.json')
    require(first['success'] and first['status']==0 and first_receipt['newIntegrandEvaluations']==2668,'actual successful saved GK21 result')
    require(second_input['rule']=='gk15' and second_input['epsabs']==5e-9 and second_input['epsrel']==5e-7 and second_input['workers']==1 and second_input['limit']==256,'unchanged unfinished comparison settings')
    source_output=second_input['sourceAction'];journal.artifacts['row-input.pickle']=second_input['rowInput']
    cache={};factor_cache={};points=[]
    completed=[int(name.split('/')[1]) for name in inventory if name.startswith('integrand/') and name.endswith('/completed.json')]
    require(sorted(completed)==list(range(1,4097)),'exact contiguous saved full point receipts')
    variable=row['row']['limits'][0][0]
    for number in sorted(completed):
        prefix='integrand/'+str(number);inp=packet(prefix+'/input.pickle');receipt=metadata(prefix+'/completed.json')
        require(receipt['input']==journal.artifacts[prefix+'/input.pickle'] and receipt['value']==journal.artifacts[prefix+'/value.pickle'],'actual immediate point receipt')
        require(prior.same(inp,{'momentum':inp['momentum'],'variable':variable,'positions':selected['basis'],'rowInput':second_input['rowInput'],'sourceAction':source_output}),'every full saved coefficient/source action/position input joins unchanged callback')
        key=float(inp['momentum']).hex();require(key not in cache,'unique complete original point owner')
        cache[key]=receipt['value'];factor_cache[key]=journal.artifacts[prefix+'/factors.pickle'];points.append({'pointIndex':number,'key':key,'input':receipt['input'],'value':receipt['value']})
    controls=metadata('controls/checks.json');require(all(v>1e-12 for v in controls.values()),'saved completed sign and measure responses')
    first_summary={k:v for k,v in first.items() if k not in ('row','intervalValues','intervals','intervalErrors')}
    first_summary.update(rowShape=list(first['row'].shape),intervalCount=len(first['intervals']),rowMaxNorm=float(prior.np.max(abs(first['row']))))
    journal.json('saved-first-rule-summary.json',first_summary)
    journal.json('saved-point-return-catalogue.json',points)
    journal.json('unfinished-comparison-scope.json',{'firstRule':first_receipt,'priorPointCount':4096,'newPointAllowance':4096,'cumulativeDistinctPointLimit':8192,'previousGuardSeconds':failure['guardWallSeconds'],'measuredFirstRuleSeconds':first['wallSeconds'],'remainingRule':{k:v for k,v in second_input.items() if k not in ('rowInput','sourceAction')},'budgetBasis':'One missing comparison only. First full rule used2667 quadrature values in15.935s; another4096-point allowance is well within the900s job budget. No frequency/domain/grid sweep or tolerance change.','partialControllerState':'The failed quad_vec call has no final return or serialized adaptive accumulator. Its unfinished controller runs against the exact saved point function; prior function evaluations and completed GK21 are not repeated. No prior controller state is claimed restored.','priorUnfinishedRuleReuseCount':'not recorded; not inferred','numericalAcceptance':False})
    # Pin every consumed caller and full physical route; preserve old logs by hash.
    def references(value):
        if isinstance(value,dict):
            if 'sha256' in value and ('logical' in value or 'path' in value):reader.retain(value.get('logical',value.get('path')),value['sha256'])
            else:
                for child in value.values():references(child)
        elif isinstance(value,(list,tuple)):
            for child in value:references(child)
    references(scope);references(row['physicalRoutes']);references(metadata('native-and-new-numerical-callers.json'));references(metadata('recovery-native-power-join.json'));references(metadata('recovery-saved-prefix-and-failure.json'))
    for p in OLD.rglob('*'):
        if p.is_file() and 'complete' not in p.relative_to(OLD).parts:reader.retain(p)
    reader.json(prior.CP,prior.CP_SHA);reader.retain(prior.READY/'complete/checks.json',prior.READY_SHA)
    for p in (Path(__file__).resolve(),M/'S11c_d_remaining_case_frequency_row_pilot_finish_plan.md',Path(recovery.__file__),M/'S11c_d_remaining_case_frequency_row_pilot_recovery_plan.md',M/'S11c_d_remaining_case_frequency_row_pilot_power_repair.md',Path(prior.__file__),M/'S11c_d_remaining_case_frequency_row_pilot_plan.md',Path(prior.io.__file__),Path(prior.storage.__file__),M.parent/'directives/S11c_d_EXPLORATORY_ACCEPTANCE.md',Path(inspect.getsourcefile(prior.integrate.quad_vec))):reader.retain(p)
    journal.json('finish-whole-helper-joins.json',ADAPTER_JOIN)
    journal.json('completed-point-and-failure-provenance.json',{'inventory':reader.retain(OLD/'completed-pilot-file-inventory.json',INVENTORY_SHA),'failure':failure,'sourceActionRecomputed':False,'completedFirstRuleRecomputed':False,'controlsRecomputed':False})
    def forbidden(*args,**kwargs):raise RuntimeError('completed source/first-rule/native science disabled in row finish')
    prior.source_matrix=prior.independent_source=forbidden
    for name in ('diff','lambdify','cancel','expand','factor','solve','gcd','resultant','integrate'):setattr(prior.sp,name,forbidden)
    prior.io.native.f.source_jets=prior.io.native.f.polynomial_basis=prior.io.native.f.BasisMomentum.prepare_basis=forbidden
    prior.io.native.Pair.__init__=prior.io.native.maps=prior.io.native.continue_pair=forbidden
    def route_input(name,requested):
        require(name=='quadrature/1/input.pickle' and prior.same(requested,second_input),'exact unchanged full unfinished quad_vec input')
        return journal.artifacts[name]
    saved={'coefficient':row['coefficient'],'source':row['source'],'variable':variable,'positions':row['positions'],'context':row['context'],'settings':row['settings'],'source_output':source_output,'selected':selected,'matrix':matrix,'nodes':source_input['nodes'],'si':scope['sourceIndex'],'source_spread':metadata('source-action/comparison.json')['maximumScaledDifference'],'cache':cache,'factor_cache':factor_cache,'number':4097,'new_evaluations':4096,'reused_evaluations':0,'responses':controls,'results':[first['row']],'costs':[first['wallSeconds']],'route_input':route_input}
    return tuple(saved[k] for k in LOCALS)


if __name__=='__main__':
    require(prior.io.digest(Path(prior.__file__))==recovery.HELPER_SHA and prior.io.digest(Path(recovery.__file__))==RECOVERY_SHA,'immutable original and principal-power continuation helpers')
    module,ADAPTER_JOIN=adapted(Path(prior.__file__).read_text())
    namespace=vars(prior);namespace['restore_remaining']=restore_remaining
    exec(compile(module,str(Path(prior.__file__))+'[saved-points-GK15-finish]','exec'),namespace)
    namespace['main']()
