#!/usr/bin/env python3
"""Supply the saved FORM manifest to the unchanged native matrix tail."""
import ast,builtins,copy,dis,json,shutil,sys,types
from pathlib import Path
import S11c_d_remaining_case_profile_response_recover as prior
h,f,m=prior.h,prior.f,prior.m
OLD=f.STORE/'s11c-remaining-case-profile-20260920/production/complete'
PLAN=f.M/'S11c_d_remaining_case_profile_response_production_recovery_plan.md'
PROOF=f.M/'S11c_d_remaining_case_profile_response_manifest_repair.json'


def matrix_tail(manifest):
    saved=json.loads((OLD/'inputs.json').read_text())
    f.require(manifest['input']==saved['input'] and manifest['settings']==h.binding.restore_settings(saved['settings'],manifest['settings']),
              'exact saved FORM input and settings passed into native matrix tail')
    native,join=h.matrix_tail()
    missing={i.argval for i in dis.get_instructions(native) if i.opname=='LOAD_GLOBAL' and i.argval not in native.__globals__ and not hasattr(builtins,i.argval)}
    f.require(missing=={'manifest'},'only missing native matrix-tail global is actual manifest')
    scope=dict(native.__globals__,manifest=manifest)
    result=types.FunctionType(native.__code__,scope,native.__name__,native.__defaults__,native.__closure__)
    f.require(result.__code__ is native.__code__ and result.__globals__['manifest'] is manifest,'unchanged whole matrix-tail bytecode and exact manifest object')
    return result,join


def load(base,resume=None):
    f.require(resume is not None,'production-only reuse of accepted focused inputs')
    evidence=json.loads(PROOF.read_text());old=json.loads((OLD/'inputs.json').read_text())
    f.require(evidence['status']=='PASSED_ACTUAL_MATRIX_TAIL_MANIFEST_REPAIR' and f.digest(OLD/'inputs.json')==evidence['originalManifestSha256'],'accepted saved-manifest routing repair')
    for n,v in old['sourceFiles'].items():f.require(f.digest(f.ROOT/n)==f.digest(OLD/'source'/n)==v,'unchanged failed production current/frozen sources')
    for n,v in old['inputPackets'].items():f.require(f.digest(Path(n))==v,'unchanged failed production original inputs')
    for n,v in evidence['originalArtifacts'].items():f.require(f.digest(OLD/n)==v['sha256'],'unchanged original completed production artifact')
    manifest,labels=prior.load(base,resume)
    for n,v in evidence['originalArtifacts'].items():f.require(f.digest(base/n)==v['sha256'],'every original production operand reproduced byte-for-byte')
    m.retain(OLD/'inputs.json',base/'failed-production-inputs.json',manifest,evidence['originalManifestSha256'])
    for path in (Path(__file__),PLAN,PROOF):
        name=str(path.resolve().relative_to(f.ROOT));value=f.digest(path);manifest['sourceFiles'][name]=value
        dest=base/'source'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,dest)
        f.require(f.digest(dest)==value,'new frozen routing repair source')
    manifest['completedProductionInputReuse']={'directory':str(OLD),'manifestSha256':evidence['originalManifestSha256'],
        'artifacts':len(evidence['originalArtifacts']),'newOriginalMatrices':0,'newOriginalSolves':0,'scope':'Original production stopped before its first context call; no integration or assembly is repeated.'}
    f.save(base/'inputs.json',manifest);return manifest,labels


def coordinator():
    tree=ast.parse(Path(h.__file__).read_text());original=h.response.function(tree,'main');node=copy.deepcopy(original)
    target=ast.parse('matrix_tail()',mode='eval').body;replacement=ast.parse('matrix_tail(manifest)',mode='eval').body
    class Replace(ast.NodeTransformer):
        def __init__(self,a,c):self.a,self.c,self.count=a,c,0
        def visit_Call(self,n):
            if ast.dump(n)==ast.dump(self.a):self.count+=1;return copy.deepcopy(self.c)
            return self.generic_visit(n)
    change=Replace(target,replacement);node=change.visit(node);f.require(change.count==1,'one production factory manifest argument')
    undo=Replace(replacement,target);restored=undo.visit(copy.deepcopy(node))
    f.require(undo.count==1 and ast.dump(restored)==ast.dump(original),'whole original FORM coordinator reverse AST join')
    return h.response.compile_function(node,dict(vars(h),load=load,matrix_tail=matrix_tail))


if __name__=='__main__':
    f.require('--resume-focused' in sys.argv and '--focused' not in sys.argv,'production recovery resumes accepted focus only')
    coordinator()()
