#!/usr/bin/env python3
"""Finish saved uniform input preparation, retaining all completed source joins."""
import ast,copy,json,shutil,types
from pathlib import Path
from types import SimpleNamespace
import S11c_d_remaining_case_uniform as u
f,m,b=u.f,u.m,u.b
OLD=f.STORE/'s11c-remaining-case-uniform-20260920/focused-retry-01/complete'
PLAN=f.M/'S11c_d_remaining_case_uniform_recovery_plan.md'
REPAIR=f.M/'S11c_d_remaining_case_uniform_evaluator_repair.json'


def pins_and_original():
    original=json.loads((OLD/'inputs.json').read_text());proof=json.loads(REPAIR.read_text())
    path=Path(u.__file__);name=str(path.relative_to(f.ROOT));f.require(f.digest(path)==proof['newHelperSha256'] and original['sourceFiles'][name]==proof['oldHelperSha256'],'exact evaluator helper transition')
    f.require(proof['wholeFileReverseAstJoin'] and proof['cutoffDependentMutationRejected'],'actual complete evaluator repair acceptance')
    for n,v in original['sourceFiles'].items():
        f.require(f.digest(OLD/'source'/n)==v,'original frozen uniform source')
        if n!=name:f.require(f.digest(f.ROOT/n)==v,'unchanged consumed uniform source')
    for n,v in original['inputPackets'].items():f.require(f.digest(Path(n))==v,'original uniform operands')
    for n,v in original['copiedInputs'].items():f.require(f.digest(OLD/n)==v,'original copies unchanged')
    pins=dict(original['sourceFiles']);pins[name]=proof['newHelperSha256']
    for p in (Path(__file__).resolve(),PLAN,REPAIR):pins[str(p.relative_to(f.ROOT))]=f.digest(p)
    return pins,original


def validate_completed(base):
    old=OLD/'source/_measurements/S11c_d_remaining_case_uniform.py'
    f.require(m.node(old,'input_routes')==m.node(Path(u.__file__),'input_routes'),'unchanged completed full input-route validator')
    routes=json.loads((base/'input-routes.json').read_text());count=0
    for address,item in routes['routes'].items():
        label,end=address.rsplit('__',1);target=base/'cases'/label/end.lower()
        f.require(f.digest(target/'mode-inputs.pickle')==item['inputSha256'] and f.digest(target/'modal.pickle')==item['modalSha256'],'actual completed route operands')
        pairs=f.unpickle(target/'uniform-input-pairs.pickle');f.require(set(pairs)==set(u.INPUT_KEYS),'full original saved input pairs')
        if item['reusedOriginalUniformResponse']:
            raw=f.unpickle(target/'uniform-original-pencil-raw.pickle');entries=f.unpickle(target/'uniform-original-pencil-certificate.pickle');final=f.unpickle(target/'uniform-original-pencil-pair.pickle')
            f.require(m.same(raw['pair'],final['pair']) and all(v==0 for v in final['residual']),'actual saved raw/proven pencil pairs')
            f.require(len(entries)==25 and all(v['expandedNumerator']==0 for v in entries),'completed exact numerator certificates')
            for index,entry in enumerate(entries):
                f.require(entry['index']==index and m.same(entry['left'],raw['pair'][0][index]) and m.same(entry['right'],raw['pair'][1][index]),'actual coefficient address in completed certificate')
            count+=len(entries)
    f.require(len(routes['routes'])==12 and count==250 and len(routes['wrongInputControls'])==4,'complete reused address and proof census')
    saved=f.unpickle(base/'new-uniform/prepared-evaluation.pickle');inputs=f.unpickle(base/'cases'/u.NEW/'right/mode-inputs.pickle')
    f.require(saved['sourceModalSha256']==f.digest(base/'cases'/u.NEW/'right/modal.pickle'),'actual completed modal preparation source')
    f.require(all(m.same(saved[k],inputs[k]) for k in ('physical','relation','binding')),'exact completed preparation inputs')
    expected=tuple(inputs['pairing']['FREQUENCY_LEGS'])+tuple(inputs['pairing']['NORMAL_LEGS'])+tuple(inputs['pairing']['BULK_LEGS'])
    f.require(tuple(saved['variables'][:-1])==expected,'actual native six-frequency/momentum coordinates')
    return {'routes':12,'zeroNumeratorCertificatesReused':count,'unchangedValidatorAst':True,'completedSymbolicPreparationReused':True}


def load(base,resume=None):
    pins,original=pins_and_original()
    manifest={'runDirectory':str(base),'sourceFiles':pins,'inputPackets':dict(original['inputPackets']),'copiedInputs':{},'input':original['input'],'scope':original['scope']}
    manifest['inputPackets'][str(OLD/'inputs.json')]=f.digest(OLD/'inputs.json')
    if resume:
        focused=json.loads((resume/'checks.json').read_text());accepted=json.loads(u.FCP.read_text())
        f.require(accepted['status']=='ACCEPTED_FOUR_CASE_UNIFORM_INPUTS' and accepted['checksSha256']==f.digest(resume/'checks.json'),'accepted completed recovered focus')
        f.require(focused['sourceFiles']==pins and focused['mode']=='focused','exact recovered focus sources')
        for n,v in focused['artifacts'].items():m.retain(resume/n,base/n,manifest,v['sha256'])
        manifest['completedFocusedReuse']={'directory':str(resume),'checksSha256':f.digest(resume/'checks.json'),'artifacts':len(focused['artifacts'])}
    else:
        for p in sorted(OLD.rglob('*')):
            if not p.is_file() or 'source' in p.relative_to(OLD).parts or p==OLD/'inputs.json':continue
            m.retain(p,base/p.relative_to(OLD),manifest)
        m.retain(OLD/'inputs.json',base/'original-inputs.json',manifest)
        m.retain(OLD/'source/_measurements/S11c_d_remaining_case_uniform.py',base/'original-uniform-helper.py',manifest)
        validation=validate_completed(base);f.save(base/'completed-input-validation-reuse.json',validation)
        manifest['completedPreparationReuse']={**validation,'directory':str(OLD),'copiedArtifacts':len(manifest['copiedInputs'])}
    for n,v in pins.items():
        dst=base/'source'/n;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/n,dst);f.require(f.digest(dst)==v,'current recovered source snapshot')
    f.save(base/'inputs.json',manifest);return manifest


def prepare(base,manifest):
    tree=ast.parse(Path(u.__file__).read_text());original=u.h.function(tree,'prepare')
    first=next(i for i,n in enumerate(original.body) if isinstance(n,ast.Assign) and isinstance(n.value,ast.Call) and isinstance(n.value.func,ast.Name) and n.value.func.id=='uniform_evaluator')
    tail=copy.deepcopy(original.body[first:])
    # Exact original per-mode replay and all guards, after the saved preparation.
    node=ast.FunctionDef(name='finish_prepared_modes',args=ast.arguments(posonlyargs=[],args=[ast.arg(n) for n in ('target','inputs','modal','saved','keys','pair')],vararg=None,kwonlyargs=[],kw_defaults=[],kwarg=None,defaults=[]),body=tail,decorator_list=[])
    f.require(ast.dump(ast.Module(body=tail,type_ignores=[]))==ast.dump(ast.Module(body=original.body[first:],type_ignores=[])),'literal full remaining preparation tail')
    target=base/'new-uniform';saved=f.unpickle(target/'prepared-evaluation.pickle');inputs=f.unpickle(base/'cases'/u.NEW/'right/mode-inputs.pickle');modal,_=f.unpickle(base/'cases'/u.NEW/'right/modal.pickle')
    symbols={str(s):s for s in inputs['physical'].free_symbols}
    pair=SimpleNamespace(modes=SimpleNamespace(k=symbols['s11cdSpectralNormalMomentum'],q=symbols['s11cdBulkRadical']))
    finish=u.h.compile_function(node,dict(vars(u)));return finish(target,inputs,modal,saved,tuple(saved['bound']),pair)


def routes(base):return json.loads((base/'input-routes.json').read_text())


if __name__=='__main__':
    namespace=dict(vars(u),load=load,input_routes=routes,prepare=prepare)
    main=types.FunctionType(u.main.__code__,namespace,u.main.__name__,u.main.__defaults__,u.main.__closure__)
    f.require(main.__code__ is u.main.__code__,'unchanged complete uniform coordinator bytecode')
    main()
