#!/usr/bin/env python3
"""Finish saved analytic validation at the native JSON owner boundary."""
import ast
import copy
import hashlib
import importlib.util
import inspect
import json
import sys
from pathlib import Path

REPO=Path('/var/projects/toy_physics');M=REPO/'research/pde_ledger_v3/_measurements'
PREVIOUS=REPO/'_scratch/s11c/s11c-remaining-case-frequency-20260921/analytic-acceptance'
ROOT=PREVIOUS.parent/'analytic-acceptance-recovery-01'
PLAN=M/'S11c_d_remaining_case_frequency_analytic_validation_finish_plan.md'
REPAIR=M/'S11c_d_remaining_case_frequency_analytic_validation_repair.json'
sys.dont_write_bytecode=True
spec=importlib.util.spec_from_file_location('saved_analytic_validator_original',PREVIOUS/'validate.py')
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
STATE={'ownerCalls':0,'controls':False}


def typed_metadata(value):
    if isinstance(value,dict):
        assert all(type(k) is str for k in value)
        return {'type':'dict','entries':{k:typed_metadata(x) for k,x in value.items()}}
    if isinstance(value,(tuple,list)):
        return {'type':type(value).__name__,'items':[typed_metadata(x) for x in value]}
    assert type(value) in (str,int,bool,type(None)),('owner metadata primitive',type(value).__name__)
    return {'type':type(value).__name__,'value':value}


def json_owner_same(owner,route_owner):
    assert isinstance(owner,dict) and isinstance(route_owner,dict)
    assert set(owner)<={'kind','case','key','field','nativeCall','address','sourceOwner','ownerAddress','path','index','branch','point','member'}
    # This is exactly the producer's JSON serializer boundary, exclusively for
    # owner metadata. Physical arrays, expressions, units and inputs still use
    # the unchanged strict typed comparator everywhere else.
    typed=typed_metadata(owner);encoded=json.loads(json.dumps(owner))
    index=STATE['ownerCalls'];STATE['ownerCalls']+=1
    v.save('owner-serialization-'+str(index)+'.json',{'rawTypedOwner':typed,'rawJsonOwner':route_owner,
        'nativeJsonDecodedOwner':encoded,'originalSerializerSha256':STATE['repair']['serializerSourceSha256'],
        'physicalOperandComparisonChanged':False})
    assert encoded==route_owner,'literal native owner JSON serialization'
    if not STATE['controls'] and 'ownerAddress' in owner:
        changed={}
        for name in ('case','ownerAddress','sourceOwner'):
            item=copy.deepcopy(route_owner)
            if name=='case':item[name]='wrong-case'
            elif name=='ownerAddress':item[name]=['wrong-address']+item[name]
            else:item[name]['address']=['wrong-source-owner']+item[name]['address']
            changed[name]=item
        v.save('owner-serialization-mutation-operands.json',{'original':route_owner,'changed':changed})
        controls={name:value!=encoded for name,value in changed.items()}
        v.save('owner-serialization-mutation-controls.json',controls);assert all(controls.values())
        STATE['controls']=True


def prefix_inputs():
    for name,expected in STATE['repair']['completedPrefixInputs'].items():v.check(name,expected)
    assert v.check(PREVIOUS/'validate.py')==STATE['repair']['originalValidatorSha256']
    old=STATE['repair']['failedOutcome']
    assert old['supervisor']['exitCode']==old['guard']['exitCode']==1 and old['guard']['childOutcome']['guardReason'] is None
    trace=(PREVIOUS/'validate.stderr').read_text()
    assert 'line 126' in trace and 'validate_routes' in trace and 'AssertionError: literal full saved operand identity' in trace


def restore_initial_atlas(roots,baseline,den,labels):
    prefix_inputs()
    directory=v.P/'analytic-proof-reuse';atlas=v.packet(directory/'saved-derivative-atlas.pickle')
    # Its whole initial source/owner/unit/value comparison completed before the
    # first route's JSON owner assertion. Restore the saved lookup result only.
    initial={key:entry['value'] for key,entry in atlas.items()}
    owners={key:entry['owner'] for key,entry in atlas.items()}
    v.save('completed-initial-derivative-atlas-reuse.json',{
        'atlasSha256':v.check(directory/'saved-derivative-atlas.pickle'),'keys':len(atlas),
        'originalValidatorSha256':STATE['repair']['originalValidatorSha256'],
        'originalPrefixAST':STATE['repair']['adapterJoins']['routes']['removedPrefixAST'],
        'completionEvidence':'Original traceback reached first derivative-route comparison after the complete initial atlas loop.',
        'sourceOrDerivativeConstructionRepeated':False,'initialAtlasComparisonRepeated':False})
    return roots['frequency'],directory,initial,owners,atlas


def restore_main_prefix(roots,den,chart,baseline,checks):
    prefix_inputs()
    assert v.check(v.P/'checks.json')==STATE['repair']['producerChecksSha256']
    v.save('completed-initial-unit-controls-reuse.json',{
        'originalHashPreflight':json.loads((PREVIOUS/'hash-preflight.json').read_text()),
        'originalValidatorSha256':STATE['repair']['originalValidatorSha256'],
        'completedBeforeFailedRoute':['root/chart source identity','frequency/carrier units','seed prefix receipt','three routing controls'],
        'physicalComparisonsRepeated':False,'producerBytesUnchanged':True})
    return roots['frequency']


class OwnerCalls(ast.NodeTransformer):
    def __init__(self):self.count=0
    def visit_Call(self,node):
        if ast.unparse(node)=="same(entry['owner'], route['owner'])":
            self.count+=1;node.func=ast.Name(id='json_owner_same',ctx=ast.Load())
        return self.generic_visit(node)


class ReverseOwnerCalls(ast.NodeTransformer):
    def visit_Call(self,node):
        if isinstance(node.func,ast.Name) and node.func.id=='json_owner_same':node.func.id='same'
        return self.generic_visit(node)


def adapter(name):
    original=ast.parse(inspect.getsource(getattr(v,name)));changed=copy.deepcopy(original);fn=changed.body[0]
    removed=None;position=None;added_count=0
    if name=='validate_routes':
        stop=next(i for i,x in enumerate(fn.body) if isinstance(x,ast.Assign) and ast.unparse(x.targets[0])=='routes')
        removed=copy.deepcopy(fn.body[:stop]);position=0
        fn.body[:stop]=ast.parse('w,directory,initial,owners,atlas=restore_initial_atlas(roots,baseline,den,labels)').body
        added_count=1
    elif name=='main':
        start=next(i for i,x in enumerate(fn.body) if isinstance(x,ast.Expr) and ast.unparse(x.value)=="same(chart['chart'], roots)")
        stop=next(i for i,x in enumerate(fn.body) if isinstance(x,ast.Assign) and ast.unparse(x.targets[0])=='counts')
        removed=copy.deepcopy(fn.body[start:stop]);position=start
        fn.body[start:stop]=ast.parse('w=restore_main_prefix(roots,den,chart,baseline,c)').body;added_count=1
        at=next(i for i,x in enumerate(fn.body) if isinstance(x,ast.Expr) and isinstance(x.value,ast.Call)
                and ast.unparse(x.value.func)=='save' and isinstance(x.value.args[0],ast.Constant) and x.value.args[0].value=='checks.json')
        fn.body[at:at]=ast.parse("result['validationContinuation']=CONTINUATION_JOIN").body
    transform=OwnerCalls();changed=transform.visit(changed)
    reverse=ReverseOwnerCalls().visit(copy.deepcopy(changed));rf=reverse.body[0]
    if name=='main':
        rf.body=[x for x in rf.body if not (isinstance(x,ast.Assign) and ast.unparse(x.targets[0])=="result['validationContinuation']")]
    if removed is not None:rf.body[position:position+added_count]=removed
    assert ast.dump(reverse)==ast.dump(original),('whole validator function reverse AST',name)
    expected={'validate_routes':2,'validate_final_proof_routes':1,'main':0};assert transform.count==expected[name]
    ns=dict(vars(v),json_owner_same=json_owner_same,restore_initial_atlas=restore_initial_atlas,restore_main_prefix=restore_main_prefix)
    ns['CONTINUATION_JOIN']=STATE.get('join')
    exec(compile(ast.fix_missing_locations(changed),'<saved-analytic-validation-continuation:'+name+'>','exec'),ns)
    return ns[name],{'wholeFunctionReverseAST':True,'originalAST':hashlib.sha256(ast.dump(original).encode()).hexdigest(),
                     'ownerMetadataCalls':transform.count,'completedPrefixRestorations':int(removed is not None),
                     'removedPrefixAST':hashlib.sha256(ast.dump(ast.Module(body=removed or [],type_ignores=[])).encode()).hexdigest()}


def main():
    repair=json.loads(REPAIR.read_text());STATE['repair']=repair
    assert v.f.digest(Path(__file__))==repair['newHelperSha256'] and v.f.digest(PLAN)==repair['newPlanSha256']
    assert v.f.digest(PREVIOUS/'validate.py')==repair['originalValidatorSha256']
    v.B=ROOT
    refs={}
    for name,item in repair['originalFiles'].items():
        original=PREVIOUS/name;v.check(original,item['sha256']);assert original.stat().st_size==item['bytes']
        path=ROOT/'original-validation'/name;path.parent.mkdir(parents=True,exist_ok=True)
        assert not path.exists() and not path.is_symlink();path.symlink_to(original)
        v.check(path,item['sha256']);refs[name]={'original':str(original),'resolved':str(path.resolve()),**item}
    for path in (Path(__file__),PLAN,REPAIR):v.check(path)
    v.save('original-validation-reference-inventory.json',refs)
    serializer=inspect.getsource(v.f.save)
    assert "json.dumps(value, indent=2)" in serializer and "write_text" in serializer
    assert v.f.digest(Path(inspect.getsourcefile(v.f.save)))==repair['serializerSourceSha256']
    STATE['join']={'originalValidator':str(PREVIOUS/'validate.py'),'originalValidatorSha256':repair['originalValidatorSha256'],
        'continuationHelperSha256':repair['newHelperSha256'],'planSha256':repair['newPlanSha256'],
        'originalFailedOutcome':repair['failedOutcome'],'ownerSerializerSourceSha256':repair['serializerSourceSha256'],
        'onlyOwnerMetadataSerialization':True,'physicalTypedComparatorUnchanged':True,'adapterJoins':repair['adapterJoins']}
    v.save('continuation-join.json',STATE['join'])
    route,rj=adapter('validate_routes');v.validate_routes=route
    proof,pj=adapter('validate_final_proof_routes');v.validate_final_proof_routes=proof
    finish,mj=adapter('main')
    assert {'routes':rj,'proofs':pj,'main':mj}==repair['adapterJoins']
    # The main namespace captures the two continuations after installation.
    finish()


if __name__=='__main__':main()
