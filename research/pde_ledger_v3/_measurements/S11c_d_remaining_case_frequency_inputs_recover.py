#!/usr/bin/env python3
"""Resume saved input loading after the repository-relative publication lookup."""
import ast
import copy
import hashlib
import inspect
import json
from pathlib import Path
import shutil

import S11c_d_remaining_case_frequency_inputs as h

f = h.f
REPO = f.ROOT.parents[1]
PREVIOUS = f.STORE/'s11c-remaining-case-frequency-20260921/inputs'
PLAN = f.M/'S11c_d_remaining_case_frequency_inputs_recovery_plan.md'
REPAIR = f.M/'S11c_d_remaining_case_frequency_inputs_path_repair.json'
original_load, original_reference = h.load, h.reference


def loader():
    original = ast.parse(inspect.getsource(original_load)); changed = copy.deepcopy(original)
    fn = changed.body[0]
    assignments = [n for n in fn.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='published' for t in n.targets)]
    f.require(len(assignments)==1,'one actual publication path assignment')
    value=assignments[0].value
    f.require(isinstance(value,ast.BinOp) and isinstance(value.op,ast.Div) and ast.unparse(value.left)=='f.ROOT','original ledger-relative lookup')
    old_left=copy.deepcopy(value.left);value.left=ast.Name(id='publication_repository_root',ctx=ast.Load())
    reverse=copy.deepcopy(changed)
    assignment=next(n for n in reverse.body[0].body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='published' for t in n.targets))
    assignment.value.left=old_left
    f.require(ast.dump(reverse)==ast.dump(original),'whole loader reverse AST with one publication root change')
    namespace=dict(vars(h),reference=reuse_reference,publication_repository_root=REPO)
    exec(compile(ast.fix_missing_locations(changed),str(Path(h.__file__)) + ':publication-root', 'exec'),namespace)
    return namespace['load'], {'wholeLoaderReverseAST':True,'publicationRootEdits':1,
        'originalLoaderAstSha256':hashlib.sha256(ast.dump(original).encode()).hexdigest(),
        'wholeOriginalMainAstSha256':h.body(h.main),'originalReferenceAstSha256':h.body(original_reference)}


def reuse_reference(base,manifest,origin,relative,target,expected):
    dest=base/target;old=PREVIOUS/'complete'/target;source=origin/relative
    if old.is_symlink():
        f.require(dest.is_symlink() and old.readlink()==dest.readlink()==source and old.resolve()==dest.resolve()==source.resolve(), 'exact completed raw and resolved reference address')
        f.require(f.digest(old)==f.digest(dest)==f.digest(source)==expected,'unchanged saved and requested reference bytes')
        manifest['inputPackets'][str(source)]=expected
        manifest['referencedInputs'][target]={'original':str(source),'resolvedOriginal':str(source.resolve()),'sha256':expected,'bytes':source.stat().st_size}
        return
    f.require(not old.exists(), 'no unknown completed regular input')
    return original_reference(base,manifest,origin,relative,target,expected)


def load(base):
    repair=json.loads(REPAIR.read_text())
    f.require(f.digest(Path(h.__file__))==repair['originalHelperSha256']==f.digest(PREVIOUS/'helper-source.py'),'original current/frozen helper unchanged')
    for name,digest in repair['originalSourcePins'].items():f.require(f.digest(Path(name))==digest,'original launch source and checkpoint identity')
    actual={str(p.relative_to(PREVIOUS/'complete')) for p in (PREVIOUS/'complete').rglob('*') if p.is_file()}
    f.require(actual==set(repair['completedReferences']),'complete saved loader file census')
    reused={}
    for name,v in repair['completedReferences'].items():
        old=PREVIOUS/'complete'/name;dest=base/name
        f.require(old.is_symlink() and str(old.readlink())==v['rawTarget'] and str(old.resolve())==v['resolvedTarget'] and f.digest(old)==v['sha256'],'original completed reference')
        f.require(not dest.exists() and not dest.is_symlink(),'fresh recovery reference')
        dest.parent.mkdir(parents=True,exist_ok=True);dest.symlink_to(old.readlink())
        reused[name]={'originalReference':str(old),**v}
    f.save(base/'completed-loader-reference-reuse.json',reused)
    logs={}
    for path in PREVIOUS.rglob('*'):
        if not path.is_file() or 'complete' in path.relative_to(PREVIOUS).parts or 'completion-watcher' in path.relative_to(PREVIOUS).parts:continue
        rel=str(path.relative_to(PREVIOUS));dest=base/'original-loader-logs'/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,dest)
        f.require(f.digest(path)==f.digest(dest),'exact preserved failed loader log/source');logs[str(path)]=f.digest(path)
    fn,join=loader();f.save(base/'loader-path-join.json',join)
    # Both raw addresses are retained. The published checksum and annex object
    # establish the repository-root route without changing a physical operand.
    cp_path=f.M/'S11c_d_remaining_case_coordinate_output_checkpoint.json';cp=json.loads(cp_path.read_text())
    raw=cp['publication']['path'];right=REPO/raw;wrong=f.ROOT/raw
    pair={'publicationPath':raw,'repositoryRoot':str(REPO),'ledgerRoot':str(f.ROOT),'actualPath':str(right),'failedPath':str(wrong),'publicationSha256':cp['publication']['sha256']}
    f.save(base/'publication-address-pair.json',pair)
    f.require(raw==repair['publicationPath'] and right.is_symlink() and not wrong.exists() and f.digest(right)==cp['publication']['sha256'],'actual published file route, no operand waiver')
    manifest,labels=fn(base)
    manifest['completedLoaderReuse']={'originalRunDirectory':str(PREVIOUS),'references':len(reused),'originalHelperSha256':repair['originalHelperSha256'],'wholeOriginalMainUnchanged':True,'noPreviouslyStartedSourceCatalogue':True}
    for name,v in reused.items():manifest['inputPackets'][v['originalReference']]=v['sha256']
    manifest['inputPackets'].update(logs)
    for p in (Path(__file__).resolve(),PLAN,REPAIR):
        name=str(p.relative_to(f.ROOT));f.require(name not in manifest['sourceFiles'],'new recovery provenance source')
        manifest['sourceFiles'][name]=f.digest(p);dest=base/'source'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
    f.save(base/'inputs.json',manifest)
    return manifest,labels


if __name__=='__main__':
    h.load=load
    h.main()
