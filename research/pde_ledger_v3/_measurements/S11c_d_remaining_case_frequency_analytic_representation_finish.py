#!/usr/bin/env python3
"""Finish saved analytic validation while retaining original live representations."""
import ast
import collections
import copy
import hashlib
import importlib.util
import inspect
import json
from pathlib import Path
import sys

REPO=Path('/var/projects/toy_physics');M=REPO/'research/pde_ledger_v3/_measurements'
FAMILY=REPO/'_scratch/s11c/s11c-remaining-case-frequency-20260921'
PREVIOUS=FAMILY/'analytic-acceptance-recovery-01'
DIAGNOSTIC=FAMILY/'analytic-representation-diagnostic-recovery-01'
ROOT=FAMILY/'analytic-acceptance-recovery-02'
PLAN=M/'S11c_d_remaining_case_frequency_analytic_representation_finish_plan.md'
REPAIR=M/'S11c_d_remaining_case_frequency_analytic_representation_repair.json'
sys.dont_write_bytecode=True
spec=importlib.util.spec_from_file_location('prior_analytic_validation_finish',M/'S11c_d_remaining_case_frequency_analytic_validation_finish.py')
prior=importlib.util.module_from_spec(spec);spec.loader.exec_module(prior)
v=prior.v
STATE={'views':{},'viewCalls':0,'controlCases':set()}


def reference(source,destination,expected):
    v.check(source,expected['sha256']);assert source.stat().st_size==expected['bytes']
    assert str(source.resolve())==expected['resolved']
    assert (str(source.readlink()) if source.is_symlink() else None)==expected['rawLink']
    destination.parent.mkdir(parents=True,exist_ok=True)
    assert not destination.exists() and not destination.is_symlink()
    destination.symlink_to(source);v.check(destination,expected['sha256'])


def completed_summary(name):
    path=PREVIOUS/name
    v.check(path,STATE['repair']['previousFiles'][name]['sha256'])
    return v.read(path)


def restore_main_prefix(roots,den,chart,baseline,checks):
    # These exact inputs, unit guards, atlas prefix and controls completed in
    # the first validator and were hash-joined by the prior continuation.
    for name,expected in prior.STATE['repair']['completedPrefixInputs'].items():v.check(name,expected)
    completed_summary('completed-initial-unit-controls-reuse.json')
    completed_summary('completed-initial-derivative-atlas-reuse.json')
    v.save('completed-initial-validation-reuse.json',{'previous':str(PREVIOUS),
        'originalValidatorSha256':prior.STATE['repair']['originalValidatorSha256'],
        'previousHelperSha256':STATE['repair']['previousHelperSha256'],'prefixChecksRepeated':False})
    return roots['frequency']


def restore_routes(roots,baseline,den,labels):
    summary=completed_summary('validated-derivative-lift-routes.json')
    directory=v.P/'analytic-proof-reuse'
    assert not v.derivative_values and not v.lift_values and not v.certificate_uses
    for i,_ in enumerate(v.read(directory/'derivative-call-routes.json')):
        path=directory/'derivative-calls'/str(i)
        op=v.packet(path/'input.pickle');entry=v.packet(path/'value.pickle')
        v.derivative_values.setdefault(op['active'],[]).append((op,entry))
    for i,_ in enumerate(v.read(directory/'lift-call-routes.json')):
        path=directory/'lift-calls'/str(i)
        op=v.packet(path/'input.pickle');entry=v.packet(path/'value.pickle')
        v.lift_values[op['active']]=(op,entry)
        for proof in entry['proofs']:
            if proof['certificate'] is not None:v.certificate_uses.append((op['active'],proof['certificate']))
    # Restore completed lookup locals directly from the already validated
    # input/value packets; no route, branch or certificate guard is repeated.
    v.save('completed-route-validation-reuse.json',{'summary':summary,
        'previousSummarySha256':v.check(PREVIOUS/'validated-derivative-lift-routes.json'),
        'derivativeValues':sum(map(len,v.derivative_values.values())),'liftValues':len(v.lift_values),
        'savedBranchConsumers':len(v.certificate_uses),'routeValidationRepeated':False,
        'scientificConstructionRepeated':False})
    return collections.Counter({k:summary[k] for k in ('newDerivatives','reusedDerivatives','newLifts','reusedLifts') if k in summary})


def restore_denominators(roots,den):
    summary=completed_summary('validated-new-denominators.json')
    result=v.packet(v.P/'remaining-case-denominator-chart.pickle')
    v.save('completed-domain-validation-reuse.json',{'summary':summary,
        'packetSha256':v.check(v.P/'remaining-case-denominator-chart.pickle'),
        'previousSummarySha256':v.check(PREVIOUS/'validated-new-denominators.json'),
        'domainValidationRepeated':False})
    return result


def restore_baseline(label,summary):
    old=completed_summary('validated-case-'+label+'.json')
    v.same({k:old[k] for k in summary},summary)
    complete=completed_summary('validated-new-records-'+label+'.json')
    assert complete=={'records':0,'seedControlComparisons':0,'proofScalarsSoFar':0}
    v.save('completed-baseline-case-reuse.json',{'case':label,'summary':old,
        'previousSummarySha256':v.check(PREVIOUS/('validated-case-'+label+'.json')),
        'caseValidationRepeated':False,'uncheckpointedNewCaseValidationClaimedComplete':False})
    return dict(summary)


def view_matches(op,path,view):
    return (str(path)==view['sourcePath'] and
        tuple(op['representationStrings'])==tuple(x['liveString'] for x in view['fields']) and
        tuple(hashlib.sha256(s.encode()).hexdigest() for s in op['representationStrings'])==tuple(op['representationSha256']) and
        op['exactEqualLive']==view['savedExactEqualLive'] and
        tuple(v.sp.srepr(x) for x in op['unit'])==tuple(view['unit']))


def representations_saved(op,path):
    index,view=STATE['views'][str(path)]
    v.check(path,view['sourceSha256']);v.check(view['recordPath'],view['recordSha256'])
    v.save('live-restored-view-'+str(index)+'.json',{'sourcePath':str(path),'sourceSha256':view['sourceSha256'],
        'diagnosticView':str(DIAGNOSTIC/('representation-pair-'+str(index)+'.json')),
        'diagnosticSha256':v.check(DIAGNOSTIC/('representation-pair-'+str(index)+'.json')),
        'savedExactEqualLive':op['exactEqualLive'],'restoredExactEqual':op['left']==op['right'],
        'liveStringsPreserved':True,'physicalTypedOperandComparisonsUnchanged':True})
    assert view_matches(op,path,view),'exact saved live string/hash/flag/unit/address join'
    # Preserve the actual restored residual equation. The original live flag
    # is checked against its immutable live record, never relabelled restored.
    v.same(op['rawResidual'],op['left']-op['right'])
    assert (op['left']==op['right'])==view['restoredExactEqual']
    if view['savedExactEqualLive']!=view['restoredExactEqual']:
        normalized=path.parent if path.parent.name=='normalized' else path.parent/'normalized'
        cert=v.packet(normalized/'certificate.pickle');mutation=v.packet(normalized/'mutation.pickle')
        v.same((cert['LEFT'],cert['RIGHT']),(op['left'],op['right']))
        assert cert['RESIDUAL']==0 and mutation['RESIDUAL']!=0
    label=view['case']
    if label not in STATE['controlCases']:
        altered={}
        for name in ('liveString','liveFlag','unit','address'):
            changed=dict(op);address=path
            if name=='liveString':changed['representationStrings']=('changed-physical-coefficient',)+tuple(op['representationStrings'][1:])
            elif name=='liveFlag':changed['exactEqualLive']=not op['exactEqualLive']
            elif name=='unit':changed['unit']=(op['unit'][0]+1,)+tuple(op['unit'][1:])
            else:address=path.parent/'wrong-owner-operands.pickle'
            altered[name]={'path':str(address),'liveStrings':changed['representationStrings'],
                'liveFlag':changed['exactEqualLive'],'unit':[v.sp.srepr(x) for x in changed['unit']]}
        v.save('representation-mutation-operands-'+label+'.json',{'originalPath':str(path),'changed':altered})
        controls={}
        for name,item in altered.items():
            changed=dict(op,representationStrings=tuple(item['liveStrings']),exactEqualLive=item['liveFlag'])
            if name=='unit':changed['unit']=(op['unit'][0]+1,)+tuple(op['unit'][1:])
            controls[name]=not view_matches(changed,Path(item['path']),view)
        v.save('representation-mutation-controls-'+label+'.json',controls);assert all(controls.values())
        STATE['controlCases'].add(label)
    STATE['viewCalls']+=1


def finish_provenance():
    assert STATE['viewCalls']==252 and len(STATE['controlCases'])==3
    return {'previousValidation':str(PREVIOUS),'previousHelperSha256':STATE['repair']['previousHelperSha256'],
        'previousFailure':STATE['repair']['previousFailure'],'originalFailure':prior.STATE['repair']['failedOutcome'],
        'diagnosticChecksSha256':STATE['repair']['diagnosticChecksSha256'],
        'helperSha256':STATE['repair']['helperSha256'],'adapterJoins':STATE['repair']['adapterJoins'],
        'liveRestoredViews':STATE['viewCalls'],'representationMutationControls':12,
        'completedRoutesDomainsBaselineReused':True,'rawLiveHistoryRetained':True,'scientificConstructionRepeated':False}


def comparison_adapter():
    original=ast.parse(inspect.getsource(v.comparison));changed=copy.deepcopy(original)
    class Change(ast.NodeTransformer):
        def __init__(self):self.count=0
        def visit_Call(self,node):
            if isinstance(node.func,ast.Name) and node.func.id=='representations':
                operand=ast.unparse(node.args[0]);assert operand in ('op','nop');self.count+=1
                node.func.id='representations_saved'
                node.args.append(ast.parse("directory/'operands.pickle'" if operand=='op' else "directory/'normalized/operands.pickle'",mode='eval').body)
            return self.generic_visit(node)
    class Reverse(ast.NodeTransformer):
        def visit_Call(self,node):
            if isinstance(node.func,ast.Name) and node.func.id=='representations_saved':node.func.id='representations';node.args.pop()
            return self.generic_visit(node)
    transform=Change();changed=transform.visit(changed);assert transform.count==2
    assert ast.dump(Reverse().visit(copy.deepcopy(changed)))==ast.dump(original)
    ns=dict(vars(v),representations_saved=representations_saved)
    exec(compile(ast.fix_missing_locations(changed),'<saved-live-restored-comparison>','exec'),ns)
    return ns['comparison'],{'wholeFunctionReverseAST':True,'representationReaderCalls':2,
        'originalAST':hashlib.sha256(ast.dump(original).encode()).hexdigest()}


def main_adapter():
    original=ast.parse(inspect.getsource(v.main));changed=copy.deepcopy(original);fn=changed.body[0]
    start=next(i for i,x in enumerate(fn.body) if isinstance(x,ast.Expr) and ast.unparse(x.value)=="same(chart['chart'], roots)")
    stop=next(i for i,x in enumerate(fn.body) if isinstance(x,ast.Assign) and ast.unparse(x.targets[0])=='counts')
    prefix=copy.deepcopy(fn.body[start:stop]);fn.body[start:stop]=ast.parse('w=restore_main_prefix(roots,den,chart,baseline,c)').body
    case=next(x for x in fn.body if isinstance(x,ast.For) and ast.unparse(x.target)=='(label, summary)')
    case.body[:0]=ast.parse("if label==n.BASELINE:\n validated[label]=restore_baseline(label,summary)\n continue").body
    at=next(i for i,x in enumerate(fn.body) if isinstance(x,ast.Expr) and isinstance(x.value,ast.Call) and ast.unparse(x.value.func)=='save'
        and isinstance(x.value.args[0],ast.Constant) and x.value.args[0].value=='checks.json')
    fn.body[at:at]=ast.parse("result['validationContinuation']=finish_provenance()").body
    reverse=copy.deepcopy(changed);rf=reverse.body[0]
    rf.body=[x for x in rf.body if not (isinstance(x,ast.Assign) and ast.unparse(x.targets[0])=="result['validationContinuation']")]
    rcase=next(x for x in rf.body if isinstance(x,ast.For) and ast.unparse(x.target)=='(label, summary)');rcase.body.pop(0)
    rf.body[start:start+1]=prefix;assert ast.dump(reverse)==ast.dump(original)
    ns=dict(vars(v),restore_main_prefix=restore_main_prefix,restore_baseline=restore_baseline,finish_provenance=finish_provenance)
    exec(compile(ast.fix_missing_locations(changed),'<saved-analytic-representation-continuation>','exec'),ns)
    return ns['main'],{'wholeFunctionReverseAST':True,'initialPrefixRestored':True,'completedBaselineCaseReused':True,
        'originalAST':hashlib.sha256(ast.dump(original).encode()).hexdigest(),
        'removedPrefixAST':hashlib.sha256(ast.dump(ast.Module(body=prefix,type_ignores=[])).encode()).hexdigest()}


def main():
    repair=json.loads(REPAIR.read_text());STATE['repair']=repair;v.B=ROOT
    v.check(Path(__file__),repair['helperSha256']);v.check(PLAN,repair['planSha256']);v.check(REPAIR)
    v.check(Path(prior.__file__),repair['previousHelperSha256'])
    v.check(Path(v.__file__),repair['originalValidatorSha256'])
    prior.STATE.update(repair=json.loads(prior.REPAIR.read_text()),ownerCalls=191,controls=True)
    v.check(prior.REPAIR,repair['previousRepairSha256'])
    v.check(prior.PLAN,prior.STATE['repair']['newPlanSha256'])
    refs={}
    for field,origin,label in (('previousFiles',PREVIOUS,'previous-validation'),('diagnosticFiles',DIAGNOSTIC,'accepted-representation-diagnostic')):
        for name,item in repair[field].items():
            reference(origin/name,ROOT/label/name,item);refs[label+'/'+name]=dict(item,original=str(origin/name))
    for name in ('validated-derivative-lift-routes.json','validated-new-denominators.json',
            'validated-new-records-LAB_HELD__RHO4_CONSTANT.json','validated-case-LAB_HELD__RHO4_CONSTANT.json'):
        reference(PREVIOUS/name,ROOT/name,repair['previousFiles'][name])
    v.save('continuation-reference-inventory.json',refs)
    assert v.check(DIAGNOSTIC/'checks.json')==repair['diagnosticChecksSha256']
    v.h.source.receipts.inspect_guard(DIAGNOSTIC,'diagnose')
    assert (DIAGNOSTIC/'checks.json').read_bytes()==(DIAGNOSTIC/'diagnose.stdout').read_bytes()
    for index,item in enumerate(v.read(DIAGNOSTIC/'representation-pairs.json')):
        view=v.read(DIAGNOSTIC/('representation-pair-'+str(index)+'.json'))
        assert item['sourcePath']==view['sourcePath'];v.check(view['sourcePath'],view['sourceSha256'])
        assert view['sourcePath'] not in STATE['views'];STATE['views'][view['sourcePath']]=(index,view)
    assert len(STATE['views'])==252
    v.validate_routes=restore_routes;v.validate_denominators=restore_denominators
    comparison,cj=comparison_adapter();v.comparison=comparison
    proof,pj=prior.adapter('validate_final_proof_routes');v.validate_final_proof_routes=proof
    run,mj=main_adapter();joins={'comparison':cj,'proofs':pj,'main':mj}
    assert joins==repair['adapterJoins'];v.save('representation-continuation-join.json',dict(joins,originalFailedValidation=repair['previousFailure']))
    run()


if __name__=='__main__':main()
