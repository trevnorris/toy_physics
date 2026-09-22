#!/usr/bin/env python3
"""Resume the unfinished row evaluator with its accepted principal-power rule."""
import ast
import copy
import inspect
import json
from pathlib import Path
import S11c_d_remaining_case_frequency_row_pilot as prior

M,F,require=prior.M,prior.F,prior.require
OLD=F/'row-pilot'
HELPER_SHA='c579b8881e890f0257fb7cebb5b50ca30dce356b3b2141f622f0b35deb7ad2ec'
INVENTORY_SHA='4ac6401a41ad9865d62773865b724a74c2b292dd0b06a3d34fc99be176f6a5b3'
BASELINE_SHA='0c350c5f9c554c2e585920745157be826daf51775e05676a0554d8670a750663'
ENGINE_SHA='5bf77729682960bb59f10c672a5fe3f67a98b0a9533da7c6189a37183596050f'
LOCALS=('coefficient','source','variable','positions','context','settings','source_output','selected','matrix','nodes','si','source_spread')


def adapted(text):
    module=ast.parse(text)
    original_evaluate=next(n for n in module.body if isinstance(n,ast.FunctionDef) and n.name=='evaluate')
    evaluate=copy.deepcopy(original_evaluate)
    power=next(n for n in ast.walk(evaluate) if isinstance(n,ast.If) and ast.unparse(n.test)=='not isinstance(expression.exp, sp.Integer)')
    parent=next(n for n in ast.walk(evaluate) if isinstance(getattr(n,'body',None),list) and power in n.body)
    offset=parent.body.index(power);parent.body.pop(offset)
    reverse=copy.deepcopy(evaluate)
    reverse_parent=next(n for n in ast.walk(reverse) if isinstance(n,ast.If) and 'operand = np.asarray' in ast.unparse(n) and n.test.__class__ is ast.BoolOp)
    reverse_parent.body.insert(offset,copy.deepcopy(power))
    require(ast.dump(reverse)==ast.dump(original_evaluate),'whole evaluator differs only by removed nonnative positive-real restriction')
    original=next(n for n in module.body if isinstance(n,ast.FunctionDef) and n.name=='main')
    fn=copy.deepcopy(original)
    begin=next(i for i,n in enumerate(fn.body) if isinstance(n,ast.Assign) and ast.unparse(n.targets[0])=='cp')
    end=next(i for i,n in enumerate(fn.body) if isinstance(n,ast.Assign) and ast.unparse(n.targets[0])=='cache')
    restoration=ast.parse(','.join(LOCALS)+' = restore(base,reader,journal)').body[0]
    fn.body[begin:end]=[restoration]
    count=next(n for n in fn.body if isinstance(n,ast.Assign) and ast.unparse(n.targets[0])=='number')
    require(ast.literal_eval(count.value)==0,'original first point address');count.value=ast.Constant(value=1)
    checks=next(i for i,n in enumerate(fn.body) if isinstance(n,ast.Assign) and ast.unparse(n.targets[0])=='checks')
    extra=ast.parse("checks.update(completedSourceActionReused=True,newSourceActions=0,acceptedPrincipalComplexPowerRule=True,originalFailurePreserved=True)").body[0]
    fn.body.insert(checks+1,extra)
    restored=copy.deepcopy(fn);restored.body.pop(checks+1)
    next(n for n in restored.body if isinstance(n,ast.Assign) and ast.unparse(n.targets[0])=='number').value=ast.Constant(value=0)
    restored.body[begin:begin+1]=copy.deepcopy(original.body[begin:end])
    require(ast.dump(restored)==ast.dump(original),'whole original main reverse identity')
    return ast.fix_missing_locations(ast.Module(body=[evaluate,fn],type_ignores=[])),{'prefixStatementRange':[begin,end],'restoredLocals':list(LOCALS),'wholeOriginalMainReverseAST':True,'wholeEvaluatorReverseAST':True,'completedSourcePrefixReexecuted':False,'pointJournalStartsAt':1}


def restore(base,reader,journal):
    inventory=reader.json(OLD/'completed-file-inventory.json',INVENTORY_SHA)
    failure=reader.json(OLD/'failed-outcome-inspection.json')
    require(failure['actualExits']==[1,1,1] and failure['guardInterruption'] is None and failure['sourceActionCompleted'] and not failure['integrandCompleted'] and not failure['quadratureStarted'],'actual completed and unfinished slots')
    for name,record in inventory.items():
        original=reader.retain(record['path'],record['sha256']);target=base/name
        require(not target.exists() and not target.is_symlink(),'fresh saved-prefix reference')
        target.parent.mkdir(parents=True,exist_ok=True);target.symlink_to(original['canonical'])
        reader.retain(target,record['sha256']);journal.artifacts[name]=record
    logs={}
    for p in OLD.rglob('*'):
        if p.is_file() and 'complete' not in p.relative_to(OLD).parts:logs[str(p.relative_to(OLD))]=reader.retain(p)
    def packet(name):return reader.packet(inventory[name]['path'],inventory[name]['sha256'])
    row=packet('row-input.pickle');source_input=packet('source-action/input.pickle');matrix=packet('source-action/value.pickle')
    scope=reader.json(inventory['pilot-scope.json']['path'],inventory['pilot-scope.json']['sha256'])
    selected=scope['physicalInput'];source_spread=failure['sourceComparison']['maximumScaledDifference']
    def references(value):
        if isinstance(value,dict):
            if 'sha256' in value and ('logical' in value or 'path' in value):reader.retain(value.get('logical',value.get('path')),value['sha256'])
            else:
                for child in value.values():references(child)
        elif isinstance(value,(list,tuple)):
            for child in value:references(child)
    references(scope);references(row['physicalRoutes'])
    cp=reader.json(prior.CP,prior.CP_SHA);reader.retain(prior.READY/'complete/checks.json',prior.READY_SHA)
    reader.retain(prior.READY/'completed-input-artifact-inventory.json',prior.INVENTORY_SHA)
    original_callers=reader.json(inventory['native-and-new-numerical-callers.json']['path'],inventory['native-and-new-numerical-callers.json']['sha256']);references(original_callers)
    # Full native source proves the existing principal-complex convention; never
    # invoke the old cached compiler, its substitutions or lambdify.
    baseline=reader.json(M/'S11c_d_frequency_matrix_checkpoint.json',BASELINE_SHA)
    native_path=M.parent/'scripts/S11c_d_mixing_scattering_sympy_audit.py'
    require(baseline['sourceFiles']['scripts/S11c_d_mixing_scattering_sympy_audit.py']==ENGINE_SHA,'accepted native numerical compiler pin')
    native_file=reader.retain(native_path,ENGINE_SHA)
    frozen_file=reader.retain(Path(baseline['runDirectory'])/'source/scripts/S11c_d_mixing_scattering_sympy_audit.py',ENGINE_SHA)
    tree=ast.parse(native_path.read_text());cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='BoundedActionQuadrature')
    compiled=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='compiled')
    numerical_power=next(n for n in ast.walk(compiled) if isinstance(n,ast.Lambda))
    require(ast.unparse(numerical_power)=='lambda x, p: np.asarray(x, dtype=complex) ** p','literal accepted principal-complex numerical power')
    for p in (Path(__file__).resolve(),M/'S11c_d_remaining_case_frequency_row_pilot_recovery_plan.md',M/'S11c_d_remaining_case_frequency_row_pilot_power_repair.md',Path(prior.__file__),M/'S11c_d_remaining_case_frequency_row_pilot_plan.md',Path(prior.io.__file__),Path(prior.storage.__file__),M.parent/'directives/S11c_d_EXPLORATORY_ACCEPTANCE.md',Path(inspect.getsourcefile(prior.integrate.quad_vec))):reader.retain(p)
    # Structural node inventory only. No coefficient, root or branch calculation.
    seen=set();powers=[];unsupported=[]
    def walk(expr):
        if id(expr) in seen:return
        seen.add(id(expr))
        sp=prior.sp
        if isinstance(expr,sp.Pow):
            if not isinstance(expr.exp,sp.Rational):unsupported.append(repr(expr))
            elif not isinstance(expr.exp,sp.Integer):powers.append({'base':expr.base,'exponent':expr.exp})
        elif not (isinstance(expr,(sp.Add,sp.Mul,sp.Rational,sp.Float,sp.Symbol)) or expr is sp.I or expr is sp.pi or expr.func in (sp.exp,sp.tanh,sp.sin,sp.cos)):
            unsupported.append(repr(expr))
        for child in expr.args:walk(child)
    walk(row['coefficient']);walk(row['source']['frequency'])
    journal.write('recovery-principal-power-inputs.pickle',{'coefficient':row['coefficient'],'sourceFrequency':row['source']['frequency'],'nonintegerPowers':powers,'fullInput':inventory['row-input.pickle']})
    journal.json('recovery-native-power-join.json',{'current':native_file,'frozen':frozen_file,'wholeNativeCompiler':ast.unparse(compiled),'literalPower':ast.unparse(numerical_power),'nonintegerPowerCount':len(powers),'unsupportedNodes':unsupported,'nativeCompilerCalled':False,'scope':'Existing principal branch of each actual expression. No new branch selection or outgoing-domain claim.'})
    require(not unsupported,'all actual unfinished evaluator nodes supported before numerical work')
    journal.json('recovery-saved-prefix-and-failure.json',{'inventory':reader.retain(OLD/'completed-file-inventory.json',INVENTORY_SHA),'failure':failure,'logs':logs,'newSourceActions':0,'newSourceCoefficientEvaluations':0,'sourceProofRecomputation':False,'originalSourceComparisonReused':source_spread,'unfinishedFullCoefficientCallOnly':True})
    journal.json('recovery-whole-helper-joins.json',ADAPTER_JOIN)
    def forbidden(*args,**kwargs):raise RuntimeError('completed source/native science disabled in saved row continuation')
    prior.source_matrix=prior.independent_source=forbidden
    for name in ('diff','lambdify','cancel','expand','factor','solve','gcd','resultant','integrate'):setattr(prior.sp,name,forbidden)
    prior.io.native.f.source_jets=prior.io.native.f.polynomial_basis=prior.io.native.f.BasisMomentum.prepare_basis=forbidden
    prior.io.native.Pair.__init__=prior.io.native.maps=prior.io.native.continue_pair=forbidden
    saved={'coefficient':row['coefficient'],'source':row['source'],'variable':row['row']['limits'][0][0],'positions':row['positions'],'context':row['context'],'settings':row['settings'],'source_output':inventory['source-action/value.pickle'],'selected':selected,'matrix':matrix,'nodes':source_input['nodes'],'si':scope['sourceIndex'],'source_spread':source_spread}
    return tuple(saved[k] for k in LOCALS)


if __name__=='__main__':
    require(prior.io.digest(Path(prior.__file__))==HELPER_SHA,'immutable original row pilot helper')
    module,ADAPTER_JOIN=adapted(Path(prior.__file__).read_text())
    namespace=vars(prior);namespace['restore']=restore
    exec(compile(module,str(Path(prior.__file__))+'[saved-prefix-principal-power-continuation]','exec'),namespace)
    namespace['main']()
