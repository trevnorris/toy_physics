#!/usr/bin/env python3
"""Finish saved certificate-route accounting and the analytic final audit."""
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
PREVIOUS=FAMILY/'analytic-acceptance-recovery-02';DIAGNOSTIC=FAMILY/'analytic-certificate-routing'
ROOT=FAMILY/'analytic-acceptance-recovery-03'
PLAN=M/'S11c_d_remaining_case_frequency_analytic_certificate_finish_plan.md'
REPAIR=M/'S11c_d_remaining_case_frequency_analytic_certificate_repair.json'
sys.dont_write_bytecode=True
spec=importlib.util.spec_from_file_location('prior_analytic_representation_finish',M/'S11c_d_remaining_case_frequency_analytic_representation_finish.py')
prior=importlib.util.module_from_spec(spec);spec.loader.exec_module(prior)
v=prior.v;owner_reader=prior.prior
STATE={'decisions':0,'ownerCalls':0}


def restore_completed_validation(checks):
    roots=v.packet(v.P/'accepted/accepted-chart/root-chart.pickle')
    baseline=v.packet(v.P/'accepted/accepted-chart/analytic-sources.pickle')
    den=v.packet(v.P/'accepted/accepted-chart/denominator-chart.pickle')
    validated={};proof_count=0
    for label,expected in checks['cases'].items():
        summary=v.read(PREVIOUS/('validated-case-'+label+'.json'))
        v.same({k:summary[k] for k in expected},expected);validated[label]=dict(expected)
        proof_count=v.read(PREVIOUS/('validated-new-records-'+label+'.json'))['proofScalarsSoFar']
        # Restore only the original consumer bookkeeping locals from already
        # validated complete records. No case/representation proof is repeated.
        folder=v.P/'new-analytic-records'/label
        if not (folder/'record-inventory.json').exists():continue
        for item in v.read(folder/'record-inventory.json').values():
            record_path=folder/item['path'];v.check(record_path,item['sha256']);record=v.packet(record_path)
            active=str(folder/'operand-checkpoints'/record_path.stem)
            for pair in record['bindingComparisons'].values():
                if pair['certificate'] is not None:
                    v.certificate_uses.extend(((active,pair['certificate']),(active,pair['mutation'])))
                for entry in pair['rootIdentities']:v.root_uses.append((active,entry['rootIdentity']))
    directory=v.P/'analytic-proof-reuse'
    for i,_ in enumerate(v.read(directory/'lift-call-routes.json')):
        op=v.packet(directory/'lift-calls'/str(i)/'input.pickle')
        entry=v.packet(directory/'lift-calls'/str(i)/'value.pickle')
        for proof in entry['proofs']:
            if proof['certificate'] is not None:v.certificate_uses.append((op['active'],proof['certificate']))
    summary=v.read(PREVIOUS/'validated-derivative-lift-routes.json')
    counts=collections.Counter({k:summary[k] for k in ('newDerivatives','reusedDerivatives','newLifts','reusedLifts') if k in summary})
    new_count=sum(x['newAnalyticImages'] for x in validated.values())
    v.save('completed-all-case-validation-reuse.json',{'cases':validated,'newRecords':new_count,'proofScalars':proof_count,
        'restoredCertificateConsumers':len(v.certificate_uses),'restoredRootConsumers':len(v.root_uses),
        'summariesSource':str(PREVIOUS),'completedCaseOrRepresentationValidationRepeated':False})
    return roots,baseline,den,validated,new_count,proof_count,counts


def restore_proof_prefix(roots,baseline,den,labels):
    trace=(PREVIOUS/'validate.stderr').read_text()
    assert 'validate_final_proof_routes' in trace and 'line 36' in trace
    # The full seed loop and initial atlas comparison precede the failed first
    # cache decision in the exact joined source. No independent summary for
    # this old, uncheckpointed prefix is invented.
    directory=v.P/'analytic-proof-reuse'
    atlas=v.packet(directory/'saved-seed-root-atlas.pickle')
    routes=v.packet(directory/'seed-root-call-routes.pickle')
    expected=v.packet(directory/'saved-certificate-atlas.pickle')
    v.save('completed-seed-and-initial-certificate-prefix-reuse.json',{
        'sourceSha256':v.check(Path(v.__file__)),'failedStderrSha256':v.check(PREVIOUS/'validate.stderr'),
        'priorCompletedOwnerViewSha256':v.check(PREVIOUS/'owner-serialization-191.json'),
        'seedCalls':len(routes),'seedIdentities':len(atlas),'certificateAtlasKeys':len(expected),
        'completionEvidence':'Exact original proof-function prefix precedes the first route membership assertion in the saved traceback.',
        'independentOldPrefixSummaryAvailable':False,'seedOrAtlasValidationRepeated':False})
    return directory,atlas,routes,expected,collections.Counter(),[]


def owner_same(owner,route_owner):
    if STATE['ownerCalls']==0:
        path=v.P/'analytic-proof-reuse/certificate-calls/0/value.pickle'
        original=v.packet(path)['owner'];route=v.read(v.P/'analytic-proof-reuse/certificate-call-routes.json')[0]['owner']
        v.same(owner,original);v.same(route_owner,route)
        v.save('completed-first-certificate-owner-reuse.json',{'valueSha256':v.check(path),
            'ownerViewSha256':v.check(PREVIOUS/'owner-serialization-191.json'),'ownerValidationRepeated':False})
    else:owner_reader.json_owner_same(owner,route_owner)
    STATE['ownerCalls']+=1


def route_matches(op,entry,route,saved_op,saved_entry,evidence):
    return (v.n.same(op,saved_op) and v.n.same(entry,saved_entry) and route==evidence['route'] and
        route['reused'] is False and entry['owner']=={'kind':'new-native-certificate','path':route['path']})


def certificate_decision(index,path,op,entry,route,expected):
    evidence=v.read(DIAGNOSTIC/('certificate-route-'+str(index)+'.json'))
    v.check(path/'input.pickle',evidence['inputSha256']);v.check(path/'value.pickle',evidence['valueSha256'])
    saved_op=v.packet(path/'input.pickle');saved_entry=v.packet(path/'value.pickle')
    v.save('certificate-live-route-'+str(index)+'.json',{
        'input':str(path/'input.pickle'),'inputSha256':evidence['inputSha256'],
        'valueSha256':evidence['valueSha256'],'diagnosticSha256':v.check(DIAGNOSTIC/('certificate-route-'+str(index)+'.json')),
        'recordedLiveReuse':route['reused'],'restoredCacheHit':evidence['restoredCacheHit'],
        'actualConsumers':[{k:x[k] for k in ('path','sha256','operandPath','operandSha256','member','point','literalNativeCall')} for x in evidence['actualConsumers']],
        'liveLookupIndependentlyReconstructed':False,'unavailableLiveMutationStringInvented':False})
    assert route_matches(op,entry,route,saved_op,saved_entry,evidence)
    assert evidence['recordedLiveReuse'] is False
    assert evidence['actualConsumers'] and all(not x['sameFullCertificate'] for x in evidence['restoredCacheCandidates'])
    # New native receipt + input-before/value-after + exact own consumer and
    # unchanged branch/caller bodies establish what ran. A restored cache key
    # cannot replay a live lookup whose original expression tree was different.
    for consumer in evidence['actualConsumers']:
        v.check(consumer['path'],consumer['sha256']);v.check(consumer['operandPath'],consumer['operandSha256'])
        v.same(v.packet(consumer['path']),entry['value'])
        operand=v.packet(consumer['operandPath'])
        assert operand['exactEqualLive']==consumer['savedExactEqualLive']
        v.same(operand['representationStrings'],tuple(consumer['liveStrings']))
        v.same(operand['representationSha256'],tuple(consumer['liveStringHashes']))
    if index==0:
        changed={
            'route':(op,entry,dict(route,reused=True)),
            'owner':(op,dict(entry,owner={'kind':'new-native-certificate','path':'wrong-owner'}),route),
            'input':(dict(op,left=op['left']+1),entry,route),
            'consumer':(op,dict(entry,value=dict(entry['value'],RESIDUAL=entry['value']['RESIDUAL']+1)),route)}
        v.save('certificate-route-mutation-operands.json',{
            'route':changed['route'][2],'owner':changed['owner'][1]['owner'],
            'inputLeft':v.sp.srepr(changed['input'][0]['left']),
            'consumerResidual':v.sp.srepr(changed['consumer'][1]['value']['RESIDUAL'])})
        controls={name:not route_matches(a,b,c,saved_op,saved_entry,evidence) for name,(a,b,c) in changed.items()}
        v.save('certificate-route-mutation-controls.json',controls);assert all(controls.values())
    STATE['decisions']+=1
    return route['reused']


def final_provenance():
    assert STATE['decisions']==STATE['ownerCalls']==20
    return {'previousValidation':str(PREVIOUS),'previousFailure':STATE['repair']['previousFailure'],
        'priorContinuationJoin':v.read(PREVIOUS/'representation-continuation-join.json'),
        'diagnosticChecksSha256':STATE['repair']['diagnosticChecksSha256'],
        'helperSha256':STATE['repair']['helperSha256'],'adapterJoins':STATE['repair']['adapterJoins'],
        'completedCasesReused':4,'completedLiveRestoredViewsReused':252,'completedRepresentationControlsReused':12,
        'actualCertificateRoutes':20,'newRouteMutationControls':4,'liveCacheLookupReconstructed':False,
        'scientificConstructionRepeated':False}


def proof_adapter():
    original=ast.parse(inspect.getsource(v.validate_final_proof_routes));changed=copy.deepcopy(original);fn=changed.body[0]
    stop=next(i for i,x in enumerate(fn.body) if isinstance(x,ast.For) and 'certificate-call-routes.json' in ast.unparse(x.iter))
    prefix=copy.deepcopy(fn.body[:stop]);fn.body[:stop]=ast.parse('directory,atlas,routes,expected,counts,covered=restore_proof_prefix(roots,baseline,den,labels)').body
    loop=fn.body[1];assert isinstance(loop,ast.For)
    at=next(i for i,x in enumerate(loop.body) if isinstance(x,ast.Assign) and ast.unparse(x.targets[0])=='reused')
    decision=copy.deepcopy(loop.body[at:at+2]);assert ast.unparse(decision[1])=="assert reused == route['reused']"
    loop.body[at:at+2]=ast.parse('reused=certificate_decision(i,path,op,entry,route,expected)').body
    class Owner(ast.NodeTransformer):
        def __init__(self):self.count=0
        def visit_Call(self,node):
            if ast.unparse(node)=="same(entry['owner'], route['owner'])":self.count+=1;node.func.id='owner_same'
            return self.generic_visit(node)
    transform=Owner();changed=transform.visit(changed);assert transform.count==1
    class ReverseOwner(ast.NodeTransformer):
        def visit_Call(self,node):
            if isinstance(node.func,ast.Name) and node.func.id=='owner_same':node.func.id='same'
            return self.generic_visit(node)
    reverse=ReverseOwner().visit(copy.deepcopy(changed));reverse.body[0].body[1].body[at:at+1]=decision
    reverse.body[0].body[:1]=prefix;assert ast.dump(reverse)==ast.dump(original)
    ns=dict(vars(v),restore_proof_prefix=restore_proof_prefix,owner_same=owner_same,certificate_decision=certificate_decision)
    exec(compile(ast.fix_missing_locations(changed),'<saved-certificate-route-continuation>','exec'),ns)
    return ns['validate_final_proof_routes'],{'wholeFunctionReverseAST':True,'ownerBoundaryCalls':1,'liveCacheDecisionReaders':1,
        'originalAST':hashlib.sha256(ast.dump(original).encode()).hexdigest(),
        'removedPrefixAST':hashlib.sha256(ast.dump(ast.Module(body=prefix,type_ignores=[])).encode()).hexdigest()}


def main_adapter():
    original=ast.parse(inspect.getsource(v.main));changed=copy.deepcopy(original);fn=changed.body[0]
    start=next(i for i,x in enumerate(fn.body) if isinstance(x,ast.Assign) and ast.unparse(x.targets[0])=='roots')
    stop=next(i for i,x in enumerate(fn.body) if isinstance(x,ast.Assign) and ast.unparse(x.targets[0])=='proof_routes')
    completed=copy.deepcopy(fn.body[start:stop])
    fn.body[start:stop]=ast.parse('roots,baseline,den,validated,new_count,proof_count,counts=restore_completed_validation(c)').body
    at=next(i for i,x in enumerate(fn.body) if isinstance(x,ast.Expr) and isinstance(x.value,ast.Call) and ast.unparse(x.value.func)=='save'
        and isinstance(x.value.args[0],ast.Constant) and x.value.args[0].value=='checks.json')
    fn.body[at:at]=ast.parse("result['validationContinuation']=final_provenance()").body
    reverse=copy.deepcopy(changed);rf=reverse.body[0]
    rf.body=[x for x in rf.body if not (isinstance(x,ast.Assign) and ast.unparse(x.targets[0])=="result['validationContinuation']")]
    rf.body[start:start+1]=completed;assert ast.dump(reverse)==ast.dump(original)
    ns=dict(vars(v),restore_completed_validation=restore_completed_validation,final_provenance=final_provenance)
    exec(compile(ast.fix_missing_locations(changed),'<saved-final-analytic-audit>','exec'),ns)
    return ns['main'],{'wholeFunctionReverseAST':True,'completedCasesAndRoutesRestored':True,
        'originalAST':hashlib.sha256(ast.dump(original).encode()).hexdigest(),
        'removedCompletedAST':hashlib.sha256(ast.dump(ast.Module(body=completed,type_ignores=[])).encode()).hexdigest()}


def main():
    repair=json.loads(REPAIR.read_text());STATE['repair']=repair;v.B=ROOT
    for p,key in ((Path(__file__),'helperSha256'),(PLAN,'planSha256'),(Path(prior.__file__),'previousHelperSha256'),
                  (Path(v.__file__),'originalValidatorSha256')):v.check(p,repair[key])
    v.check(REPAIR);owner_reader.STATE.update(repair=json.loads(owner_reader.REPAIR.read_text()),ownerCalls=192,controls=True)
    v.check(prior.REPAIR,repair['previousRepairSha256']);previous_repair=v.read(prior.REPAIR)
    v.check(prior.PLAN,previous_repair['planSha256'])
    v.check(Path(owner_reader.__file__),previous_repair['previousHelperSha256'])
    v.check(owner_reader.REPAIR,previous_repair['previousRepairSha256'])
    v.check(owner_reader.PLAN,owner_reader.STATE['repair']['newPlanSha256'])
    refs={}
    for key,origin,label in (('previousFiles',PREVIOUS,'previous-validation'),('diagnosticFiles',DIAGNOSTIC,'accepted-certificate-diagnostic')):
        for name,item in repair[key].items():
            prior.reference(origin/name,ROOT/label/name,item);refs[label+'/'+name]=dict(item,original=str(origin/name))
    for path in PREVIOUS.glob('validated-*.json'):
        prior.reference(path,ROOT/path.name,repair['previousFiles'][path.name])
    v.save('continuation-reference-inventory.json',refs)
    v.h.source.receipts.inspect_guard(DIAGNOSTIC,'diagnose')
    v.check(DIAGNOSTIC/'checks.json',repair['diagnosticChecksSha256'])
    assert (DIAGNOSTIC/'checks.json').read_bytes()==(DIAGNOSTIC/'diagnose.stdout').read_bytes()
    joins=v.read(DIAGNOSTIC/'native-certificate-route-joins.json')
    for item in joins.values():v.check(item['sourceFile'],item['sourceSha256'])
    proof,pj=proof_adapter();v.validate_final_proof_routes=proof
    run,mj=main_adapter();assert {'proofs':pj,'main':mj}==repair['adapterJoins']
    v.save('certificate-continuation-join.json',dict(repair['adapterJoins'],nativeSourceJoins=joins))
    run()


if __name__=='__main__':main()
