#!/usr/bin/env python3
"""Continue saved certificate validation across its JSON path boundary."""
import ast
import copy
import hashlib
import importlib.util
import inspect
import json
from pathlib import Path
import sys

REPO=Path('/var/projects/toy_physics');M=REPO/'research/pde_ledger_v3/_measurements'
FAMILY=REPO/'_scratch/s11c/s11c-remaining-case-frequency-20260921'
PREVIOUS=FAMILY/'analytic-acceptance-recovery-03'
ROOT=FAMILY/'analytic-acceptance-recovery-04'
PLAN=M/'S11c_d_remaining_case_frequency_analytic_path_finish_plan.md'
REPAIR=M/'S11c_d_remaining_case_frequency_analytic_path_repair.json'
REUSED=('completed-all-case-validation-reuse.json',
        'completed-seed-and-initial-certificate-prefix-reuse.json',
        'completed-first-certificate-owner-reuse.json','certificate-live-route-0.json')
STATE={'writes':{}}
sys.dont_write_bytecode=True
spec=importlib.util.spec_from_file_location('prior_analytic_certificate_finish',
    M/'S11c_d_remaining_case_frequency_analytic_certificate_finish.py')
prior=importlib.util.module_from_spec(spec);spec.loader.exec_module(prior)
v=prior.v
original_save=v.save;original_provenance=prior.final_provenance


def decision_adapter():
    original=ast.parse(inspect.getsource(prior.certificate_decision))
    class Paths(ast.NodeTransformer):
        def __init__(self):self.fields=[]
        def visit_Call(self,node):
            if ast.unparse(node.func)=='v.packet' and len(node.args)==1:
                arg=node.args[0]
                if (isinstance(arg,ast.Subscript) and isinstance(arg.value,ast.Name)
                        and arg.value.id=='consumer' and isinstance(arg.slice,ast.Constant)
                        and arg.slice.value in ('path','operandPath')):
                    self.fields.append(arg.slice.value)
                    node.args[0]=ast.Call(func=ast.Name(id='Path',ctx=ast.Load()),args=[arg],keywords=[])
            return self.generic_visit(node)
    transform=Paths();changed=transform.visit(copy.deepcopy(original))
    assert transform.fields==['path','operandPath']
    class Reverse(ast.NodeTransformer):
        def __init__(self):self.count=0
        def visit_Call(self,node):
            if ast.unparse(node.func)=='v.packet' and len(node.args)==1:
                arg=node.args[0]
                if isinstance(arg,ast.Call) and ast.unparse(arg.func)=='Path':
                    assert ast.unparse(arg.args[0]) in ("consumer['path']","consumer['operandPath']")
                    node.args[0]=arg.args[0];self.count+=1
            return self.generic_visit(node)
    reverse=Reverse();restored=reverse.visit(copy.deepcopy(changed))
    assert reverse.count==2 and ast.dump(restored)==ast.dump(original)
    namespace=dict(vars(prior))
    exec(compile(ast.fix_missing_locations(changed),'<saved-certificate-path-reader>','exec'),namespace)
    return namespace['certificate_decision'],{'wholeFunctionReverseAST':True,'pathFields':transform.fields,
        'originalAST':hashlib.sha256(ast.dump(original).encode()).hexdigest(),
        'physicalOperandsOrGuardsChanged':False}


def save(name,value):
    path=ROOT/name
    assert not Path(name).is_absolute() and '..' not in Path(name).parts
    assert all(not p.is_symlink() for p in (path.parent,*path.parent.parents))
    if name in REUSED:
        assert name not in STATE['writes']
        expected=STATE['repair']['previousFiles'][name]
        v.check(PREVIOUS/name,expected['sha256']);v.check(path,expected['sha256'])
        assert path.is_symlink() and path.readlink()==PREVIOUS/name
        v.same(value,v.read(PREVIOUS/name))
        # All four saved views contain primitive JSON values. The exact native
        # serializer must reproduce their complete bytes before skipping a write.
        encoded=(json.dumps(value,indent=2)+'\n').encode()
        assert encoded==(PREVIOUS/name).read_bytes()
        STATE['writes'][name]={'source':str(PREVIOUS/name),'sha256':expected['sha256'],
            'fullTypedRequestedSavedIdentity':True,'nativeSerializedByteIdentity':True,'writeSkipped':True}
        original_save('path-boundary-reused-write-'+str(len(STATE['writes']))+'.json',STATE['writes'][name])
    else:original_save(name,value)


def final_provenance():
    assert set(STATE['writes'])==set(REUSED)
    result=original_provenance()
    result['pathBoundaryContinuation']={'directory':str(ROOT),'previousDirectory':str(PREVIOUS),
        'previousFailure':STATE['repair']['previousFailure'],'helperSha256':STATE['repair']['helperSha256'],
        'planSha256':STATE['repair']['planSha256'],'repairSha256':STATE['repairSha256'],
        'readerJoin':STATE['repair']['readerJoin'],'preservedPreviousFiles':len(STATE['repair']['previousFiles']),
        'completedMetadataWritesReused':STATE['writes'],'scientificConstructionRepeated':False}
    return result


def main():
    repair=json.loads(REPAIR.read_text());STATE['repair']=repair
    prior.ROOT=ROOT;v.B=ROOT
    for path,key in ((Path(__file__),'helperSha256'),(PLAN,'planSha256'),
                     (Path(prior.__file__),'previousHelperSha256'),(prior.PLAN,'previousPlanSha256'),
                     (prior.REPAIR,'previousRepairSha256')):v.check(path,repair[key])
    STATE['repairSha256']=v.check(REPAIR)
    refs={}
    for name,item in repair['previousFiles'].items():
        prior.prior.reference(PREVIOUS/name,ROOT/'previous-path-validation'/name,item)
        refs[name]=dict(item,original=str(PREVIOUS/name))
    for name in REUSED:prior.prior.reference(PREVIOUS/name,ROOT/name,repair['previousFiles'][name])
    original_save('path-boundary-reference-inventory.json',refs)
    trace=(PREVIOUS/'validate.stderr').read_text()
    assert "AttributeError: 'str' object has no attribute 'open'" in trace
    assert "v.same(v.packet(consumer['path']),entry['value'])" in trace
    assert not (PREVIOUS/'checks.json').exists()
    for name in ('active.json','validate.invocation.json','resource-guard/outcome.json','resource-guard/child-outcome.json'):
        assert v.read(PREVIOUS/name)['exitCode']==1
    assert v.read(PREVIOUS/'resource-guard/child-outcome.json')['guardReason'] is None
    decision,join=decision_adapter();assert join==repair['readerJoin']
    original_save('path-boundary-continuation-join.json',dict(join,
        previousHelperSha256=repair['previousHelperSha256'],previousFailure=repair['previousFailure'],
        originalMainUnchanged=True,originalProofAdapterUnchanged=True))
    v.save=save;prior.certificate_decision=decision;prior.final_provenance=final_provenance
    prior.main()


if __name__=='__main__':main()
