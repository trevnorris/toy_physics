#!/usr/bin/env python3
"""Resume the boundary loader with a distinct recovery bookkeeping filename."""
import ast
import copy
import json
from pathlib import Path
import shutil

import S11c_d_remaining_case_coordinate_boundary_recover as prior

f, modes, matrices, original = prior.f, prior.modes, prior.matrices, prior.original
PREVIOUS = f.STORE/'s11c-remaining-case-coordinate-20260921/boundary/production-recovery-01/complete'
PLAN = f.M/'S11c_d_remaining_case_coordinate_boundary_loader_recovery_plan.md'
REPAIR = f.M/'S11c_d_remaining_case_coordinate_boundary_loader_repair.json'


def loader_adapter():
    old = original.tree(prior.load); new = copy.deepcopy(old)
    targets = [n for n in ast.walk(new) if isinstance(n,ast.Constant) and n.value=='completed-input-reuse.json']
    f.require(len(targets)==1,'one recovery bookkeeping filename in loader')
    targets[0].value='boundary-recovery-input-reuse.json'
    reverse=copy.deepcopy(new)
    targets=[n for n in ast.walk(reverse) if isinstance(n,ast.Constant) and n.value=='boundary-recovery-input-reuse.json']
    f.require(len(targets)==1,'one distinct output filename');targets[0].value='completed-input-reuse.json'
    f.require(ast.dump(reverse)==ast.dump(old),'whole loader reverse AST: one bookkeeping filename')
    return original.compile_function(new,vars(prior)), {'originalLoaderAST':original.ast_sha(old),
        'repairedLoaderAST':original.ast_sha(new),'reverseAST':True,'filenameEdits':1,
        'priorMain':matrices.body(prior.main),'originalMain':matrices.body(original.main),
        'finiteAndCoordinatorRepair':matrices.body(prior.repaired_functions)}


def main():
    compiled, joins=loader_adapter();repair=json.loads(REPAIR.read_text())
    f.require(joins==repair['staticJoins'],'runtime loader join equals reviewed source proof')
    native_main=prior.main
    def load(base,resume):
        last=json.loads((PREVIOUS/'inputs.json').read_text())
        f.require(f.digest(PREVIOUS/'inputs.json')==repair['previousManifestSha256'], 'preserved previous loader manifest')
        for n,s in last['sourceFiles'].items():
            f.require(f.digest(f.ROOT/n)==f.digest(PREVIOUS/'source'/n)==s,'previous current/frozen helper and native joins')
        mismatches={n:{'expected':s,'actual':f.digest(PREVIOUS/n)} for n,s in last['copiedInputs'].items() if f.digest(PREVIOUS/n)!=s}
        f.require(mismatches==repair['previousCopyMismatches'],'only the recorded inherited JSON collision')
        keep=modes.retain;routes=[]
        def reuse(src,dst,manifest,expected=None):
            src,dst=Path(src),Path(dst);name=str(dst.relative_to(base));requested=f.digest(src)
            f.require(expected is None or requested==expected,'exact requested producer input')
            previous=PREVIOUS/name
            if previous.is_file() and f.digest(previous)==requested:
                keep(previous,dst,manifest,requested)
                manifest['inputPackets'][str(src)]=requested
                routes.append({'destination':name,'requestedSource':str(src),'reusedSource':str(previous),'sha256':requested})
            else:
                f.require(name=='completed-input-reuse.json' and requested==mismatches[name]['expected'],
                          'only inherited bookkeeping restored from its immutable original producer')
                keep(src,dst,manifest,requested)
                routes.append({'destination':name,'requestedSource':str(src),'sha256':requested,
                               'previousCollisionSha256':mismatches[name]['actual']})
        modes.retain=reuse
        try:manifest,labels,accepted=compiled(base,resume)
        finally:modes.retain=keep
        f.require(len(routes)==len(last['copiedInputs'])==221 and
                  sum('reusedSource' in v for v in routes)==220,'all intact previous loader copies reused')
        keep(PREVIOUS/'inputs.json',base/'previous-recovery-inputs.json',manifest,repair['previousManifestSha256'])
        keep(PREVIOUS/'completed-input-reuse.json',base/'previous-recovery-diagnostics/completed-input-reuse.json',
             manifest,mismatches['completed-input-reuse.json']['actual'])
        for n,s in repair['previousLogs'].items():
            keep(PREVIOUS.parent/n,base/'previous-recovery-logs'/n,manifest,s)
        for p in (Path(__file__).resolve(),PLAN,REPAIR):
            manifest['sourceFiles'][str(p.relative_to(f.ROOT))]=f.digest(p)
        for n,s in manifest['sourceFiles'].items():
            target=base/'source'/n
            if target.exists():f.require(f.digest(target)==s,'previous source snapshot kept byte-for-byte')
            else:
                target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f.ROOT/n,target)
            f.require(f.digest(target)==s,'explicit new loader helper source join')
        manifest['completedRecoveryLoaderReuse']={'directory':str(PREVIOUS),'previousManifestSha256':repair['previousManifestSha256'],
            'copiedRoutes':len(routes),'intactCopiesReused':220,'restoredInheritedMetadata':1,
            'newBookkeepingFile':'boundary-recovery-input-reuse.json','newScientificConstructionInPreviousRecovery':0}
        f.save(base/'inputs.json',manifest)
        f.save(base/'boundary-loader-recovery.json',{'joins':joins,'reuse':manifest['completedRecoveryLoaderReuse'],'routes':routes})
        matrices.hash_check(base,manifest)
        return manifest,labels,accepted
    prior.load=load
    f.require(prior.main is native_main,'prior recovery main callable unchanged')
    native_main()


if __name__=='__main__':main()
