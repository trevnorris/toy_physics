#!/usr/bin/env python3
"""Resume FORM input checks with exact unordered profile-cache joins."""
import ast,copy,json,shutil,types
from pathlib import Path
from types import SimpleNamespace
import S11c_d_remaining_case_profile_response as h
f,m=h.f,h.m
OLD=f.STORE/'s11c-remaining-case-profile-20260920/focused/complete'
PLAN=f.M/'S11c_d_remaining_case_profile_response_recovery_plan.md'
REPAIR=f.M/'S11c_d_remaining_case_profile_response_cache_repair.json'


def cache_evidence(left,right):
    """Compare complete cache entries; never reorder limits within an integral."""
    a,c=left['profiles'],right['profiles']
    f.require(type(a) is list and type(c) is list and len(a)==len(c),'complete actual profile cache lists')
    for entries,units in ((a,left['profileUnits']),(c,right['profileUnits'])):
        f.require(len(entries)==len(units),'one unit entry per actual cached profile')
        for i,entry in enumerate(entries):
            f.require(set(entry)=={'bound','unit'},'actual complete profile cache entry schema')
            f.require(not any(h.binding.native.same(entry['bound'],earlier['bound']) for earlier in entries[:i]),'unique actual integral cache keys')
            matches=[v for key,v in units.items() if h.binding.native.same(key,entry['bound'])]
            f.require(len(matches)==1 and h.binding.native.same(matches[0],entry['unit']),'actual cache entry joins complete unit map')
    mapping=[]
    for entry in a:
        matches=[i for i,other in enumerate(c) if h.binding.native.same(entry,other)]
        f.require(len(matches)==1,'exact original integral and unit match');mapping.append(matches[0])
    f.require(set(mapping)==set(range(len(c))) and h.binding.native.same(left['profileUnits'],right['profileUnits']),'bijective full profile cache and unit equality')
    return {'left':a,'right':c,'leftUnits':left['profileUnits'],'rightUnits':right['profileUnits'],'mapping':mapping,
            'originalListEquality':h.binding.native.same(a,c),'fullLiteralEntryMatches':len(mapping),'fullUnitMapEquality':True,
            'scope':'Cache enumeration order only; original symbolic integrands, fields, coefficients, units and every nested ordered limit remain exact.'}


def cache_join(left,right,base,address):
    directory=base/'profile-cache-joins';directory.mkdir(exist_ok=True)
    name='--'.join(address);path=directory/(name+'.pickle');raw=directory/(name+'-raw.pickle')
    actual={'left':left['profiles'],'right':right['profiles'],'leftUnits':left['profileUnits'],'rightUnits':right['profileUnits']}
    if path.exists():
        prior=f.unpickle(raw);f.require(m.same(actual,prior),'same actual profile cache operands at every reused row');return True
    f.atomic_pickle(raw,actual);result=cache_evidence(left,right);f.atomic_pickle(path,result);return True


def repaired_matrix_preflight():
    tree=ast.parse(Path(h.matrices.__file__).read_text());original=h.response.function(tree,'preflight');node=copy.deepcopy(original)
    expected=ast.parse("binding.same(bound['profiles'],old['bound']['profiles'])",mode='eval').body
    replacement=ast.parse("profile_cache_join(bound,old['bound'],base,(label,item['fromCase']))",mode='eval').body
    class Replace(ast.NodeTransformer):
        def __init__(self,a,c):self.a,self.c,self.count=a,c,0
        def visit_Call(self,n):
            if ast.dump(n)==ast.dump(self.a):self.count+=1;return copy.deepcopy(self.c)
            return self.generic_visit(n)
    change=Replace(expected,replacement);node=change.visit(node);f.require(change.count==1,'one actual profile-cache comparison adapter')
    undo=Replace(replacement,expected);restored=undo.visit(copy.deepcopy(node))
    f.require(undo.count==1 and ast.dump(restored)==ast.dump(original),'whole preflight reverse AST join')
    return h.response.compile_function(node,dict(vars(h.matrices),profile_cache_join=cache_join))


def load(base,resume=None):
    original=json.loads((OLD/'inputs.json').read_text());proof=json.loads(REPAIR.read_text())
    f.require(proof['status']=='PASSED_ACTUAL_PROFILE_CACHE_ORDER_REPAIR' and proof['wholePreflightReverseAstJoin'],'accepted actual cache comparison repair')
    for name,value in original['sourceFiles'].items():f.require(f.digest(f.ROOT/name)==f.digest(OLD/'source'/name)==value,'unchanged original/current/frozen FORM source')
    for name,value in original['inputPackets'].items():f.require(f.digest(Path(name))==value,'original input pre/post')
    for name,value in proof['originalArtifacts'].items():f.require(f.digest(OLD/name)==value['sha256'],'every original completed input artifact')
    f.require(proof['originalManifestSha256']==f.digest(OLD/'inputs.json'),'exact original completed loading')
    pins=dict(original['sourceFiles'])
    for path in (Path(__file__),PLAN,REPAIR):pins[str(path.resolve().relative_to(f.ROOT))]=f.digest(path)
    manifest=dict(original,runDirectory=str(base),sourceFiles=pins,inputPackets=dict(original['inputPackets']),copiedInputs={})
    if resume:
        focused=json.loads((resume/'checks.json').read_text());accepted=json.loads(h.FCP.read_text())
        f.require(accepted['status']=='ACCEPTED_CASE_FORM_RESPONSE_INPUTS' and accepted['checksSha256']==f.digest(resume/'checks.json'),'accepted recovered focused check')
        f.require(focused['sourceFiles']==pins and focused['mode']=='focused','exact recovered focused sources')
        for name,value in focused['artifacts'].items():m.retain(resume/name,base/name,manifest,value['sha256'])
        manifest['inputPackets'][str(resume/'checks.json')]=f.digest(resume/'checks.json');manifest['inputPackets'][str(h.FCP)]=f.digest(h.FCP)
        manifest['completedFocusedReuse']={'directory':str(resume),'checksSha256':f.digest(resume/'checks.json'),'artifacts':len(focused['artifacts'])}
    else:
        for name,value in proof['originalArtifacts'].items():m.retain(OLD/name,base/name,manifest,value['sha256'])
        m.retain(OLD/'inputs.json',base/'original-inputs.json',manifest)
        manifest['completedInputReuse']={'directory':str(OLD),'manifestSha256':f.digest(OLD/'inputs.json'),'artifacts':len(proof['originalArtifacts']),'originalCopiedInputs':len(original['copiedInputs']),'newBindings':0,'newMatrices':0,'newSolves':0}
    manifest['inputPackets'][str(OLD/'inputs.json')]=f.digest(OLD/'inputs.json')
    for name,value in pins.items():
        dest=base/'source'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/name,dest);f.require(f.digest(dest)==value,'current frozen recovery sources')
    manifest['settings']=h.binding.restore_settings(manifest['settings'],f.unpickle(base/'interiors/accepted-finite-system.pickle')['settings'])
    f.save(base/'inputs.json',manifest);return manifest,tuple(json.loads(h.BCP.read_text())['cases'])


def coordinator():
    native=repaired_matrix_preflight();proxy=SimpleNamespace(**dict(vars(h.matrices),preflight=native))
    check=types.FunctionType(h.preflight.__code__,dict(vars(h),matrices=proxy),h.preflight.__name__,h.preflight.__defaults__,h.preflight.__closure__)
    namespace=dict(vars(h),load=load,preflight=check)
    result=types.FunctionType(h.main.__code__,namespace,h.main.__name__,h.main.__defaults__,h.main.__closure__)
    f.require(result.__code__ is h.main.__code__ and check.__code__ is h.preflight.__code__,'unchanged complete FORM coordinator and input-check bytecode')
    return result


if __name__=='__main__':coordinator()()
