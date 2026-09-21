#!/usr/bin/env python3
"""Resume saved boundary inputs after one finite-basis matrix-product repair."""
import ast
import copy
import hashlib
import json
from pathlib import Path
import shutil

import S11c_d_remaining_case_coordinate_boundary as original

f, modes, matrices = original.f, original.modes, original.matrices
OLD = f.STORE/'s11c-remaining-case-coordinate-20260921/boundary/production/complete'
PLAN = f.M/'S11c_d_remaining_case_coordinate_boundary_recovery_production_plan.md'
REPAIR = f.M/'S11c_d_remaining_case_coordinate_boundary_basis_repair.json'
WRITE_REUSE = ('material-families/LEFT/input.pickle',
               'material-families/LEFT/saved-chart-map-reuse.pickle',
               'material-families/LEFT/finite/input.pickle')


def repaired_functions():
    finite = original.tree(original.finite_construct); repaired = copy.deepcopy(finite)
    targets = [n for n in ast.walk(repaired) if isinstance(n, ast.BinOp)
               and isinstance(n.op, ast.Mult) and ast.unparse(n.left) == 'U' and ast.unparse(n.right) == 'a']
    f.require(len(targets) == 1, 'single full-basis elementwise operation repair')
    targets[0].op = ast.MatMult()
    mkdir = [n for n in ast.walk(repaired) if isinstance(n, ast.Call) and ast.unparse(n.func) == 'folder.mkdir']
    f.require(len(mkdir) == 1 and ast.unparse(mkdir[0]) == 'folder.mkdir(parents=True)', 'single saved finite directory')
    mkdir[0].keywords.append(ast.keyword(arg='exist_ok', value=ast.Constant(True)))
    reverse = copy.deepcopy(repaired)
    target = next(n for n in ast.walk(reverse) if isinstance(n,ast.DictComp)
                  and isinstance(n.value,ast.BinOp) and ast.unparse(n.value) == 'U @ a')
    target.value.op = ast.Mult()
    next(n for n in ast.walk(reverse) if isinstance(n,ast.Call) and ast.unparse(n.func) == 'folder.mkdir').keywords.pop()
    f.require(ast.dump(reverse) == ast.dump(finite), 'whole finite adapter reverse AST: one product and saved directory')
    construct = original.tree(original.construct); resumed = copy.deepcopy(construct)
    mkdir2 = [n for n in ast.walk(resumed) if isinstance(n,ast.Call) and ast.unparse(n.func) == 'directory.mkdir']
    f.require(len(mkdir2) == 1 and ast.unparse(mkdir2[0]) == 'directory.mkdir(parents=True)', 'single saved family directory')
    mkdir2[0].keywords.append(ast.keyword(arg='exist_ok',value=ast.Constant(True)))
    reverse2 = copy.deepcopy(resumed)
    next(n for n in ast.walk(reverse2) if isinstance(n,ast.Call) and ast.unparse(n.func) == 'directory.mkdir').keywords.pop()
    f.require(ast.dump(reverse2) == ast.dump(construct), 'whole coordinator reverse AST: saved family directory only')
    finite_fn = original.compile_function(repaired, vars(original))
    construct_fn = original.compile_function(resumed, dict(vars(original), finite_construct=finite_fn))
    return finite_fn, construct_fn, {'finiteOriginalAST': original.ast_sha(finite),
        'finiteRepairedAST': original.ast_sha(repaired), 'finiteReverseAST': True,
        'constructOriginalAST': original.ast_sha(construct), 'constructResumedAST': original.ast_sha(resumed),
        'constructReverseAST': True, 'productEdits': 1, 'savedDirectoryEdits': 2,
        'nativeAdaptersUnchanged': matrices.body(original.native_adapters),
        'originalMainUnchanged': matrices.body(original.main)}


def load(base, resume):
    old = json.loads((OLD/'inputs.json').read_text())
    repair = json.loads(REPAIR.read_text())
    f.require(f.digest(Path(original.__file__)) == repair['originalHelperSha256'] and
              f.digest(OLD/'inputs.json') == repair['originalManifestSha256'], 'exact failed constructor and completed input manifest')
    origin, focus, checkpoint = matrices.accepted(original.CP, 'ACCEPTED_CASE_MATERIAL_BOUNDARY_INPUTS')
    f.require(origin == resume and checkpoint['checksSha256'] == old['completedFocusReuse']['checksSha256'],
              'same accepted focused inputs')
    manifest = dict(old, runDirectory=str(base), sourceFiles=dict(old['sourceFiles']),
                    inputPackets=dict(old['inputPackets']), copiedInputs={})
    for name, sha in old['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name) == f.digest(OLD/'source'/name) == sha, 'unchanged original/current/frozen consumed source')
    for name, item in repair['completedFiles'].items():
        modes.retain(OLD/name, base/name, manifest, item['sha256'])
    modes.retain(OLD/'inputs.json', base/'original-production-inputs.json', manifest, repair['originalManifestSha256'])
    for name, sha in repair['originalLogs'].items():
        modes.retain(OLD.parent/name, base/'original-production-logs'/name, manifest, sha)
    modes.retain(Path(original.__file__), base/'original-boundary-source.py', manifest, repair['originalHelperSha256'])
    for name, item in focus['artifacts'].items():
        f.require(f.digest(base/name) == item['sha256'], 'all 191 accepted focused artifacts retained exactly')
    for path in (Path(__file__).resolve(), PLAN, REPAIR):
        manifest['sourceFiles'][str(path.relative_to(f.ROOT))] = f.digest(path)
    for name, sha in manifest['sourceFiles'].items():
        target = base/'source'/name; target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(f.ROOT/name,target);f.require(f.digest(target)==sha,'recovery source snapshot')
    manifest['completedInputReuse'] = {'directory':str(OLD),'originalManifestSha256':repair['originalManifestSha256'],
        'files':len(repair['completedFiles']),'focusArtifacts':len(focus['artifacts']),
        'savedInputCheckpoints':list(WRITE_REUSE),
        'originalScienceDisposition':'Stopped before any finite current contraction or boundary coordinate solve; no complete candidate transform packet existed.'}
    f.save(base/'inputs.json',manifest);f.save(base/'completed-input-reuse.json',manifest['completedInputReuse'])
    matrices.hash_check(base,manifest)
    return manifest,tuple(focus['result']['cases']),focus['result']


def main():
    native_main = original.main
    finite, construct, joins = repaired_functions()
    f.require(joins == json.loads(REPAIR.read_text())['staticReverseAST'],
              'actual runtime adapter joins match the reviewed static source proof')
    original.finite_construct = finite; original.construct = construct
    original.load = load
    # Only three already saved input writes may be repeated by the coordinator.
    # They must join their entire requested payload and remain byte-identical.
    old_write, old_save = f.atomic_pickle, f.save
    def retain_input(path,value):
        target = Path(path)
        if target.exists():
            matches = [n for n in WRITE_REUSE if str(target).endswith('/'+n)]
            f.require(len(matches)==1,'no other existing packet write permitted')
            n = matches[0]; root = Path(str(target)[:-len(n)])
            expected = json.loads(REPAIR.read_text())['completedFiles'][n]['sha256']
            saved = f.unpickle(target);before = f.digest(target);comparable = value
            source_join = None
            if n == 'material-families/LEFT/saved-chart-map-reuse.pickle':
                original_source = OLD/'baseline-material/left-material-boundary.pickle'
                current_source = root/'baseline-material/left-material-boundary.pickle'
                f.require(saved['source'] == str(original_source) and value['source'] == str(current_source),
                          'exact original/copied chart-map source addresses')
                source_sha = json.loads(REPAIR.read_text())['completedFiles']['baseline-material/left-material-boundary.pickle']['sha256']
                f.require(f.digest(original_source) == f.digest(current_source) == source_sha,
                          'saved chart-map source copied byte-for-byte')
                comparable = dict(value, source=saved['source'])
                source_join = {'original':str(original_source),'copied':str(current_source),'sha256':source_sha}
            folder = root/'checkpoint-write-joins';folder.mkdir(exist_ok=True)
            old_write(folder/(n.replace('/','__')+'.pickle'),{'original':saved,'requested':value,'sourceJoinedRequest':comparable,
                      'sourceAddressJoin':source_join,'source':str(OLD/n),'sha256':before})
            f.require(before == expected and modes.same(saved,comparable) and f.digest(target)==before,
                      'entire original/requested saved input equality before byte reuse')
            return
        return old_write(target,value)
    def retain_native_joins(path,value):
        target = Path(path)
        if target.name == 'material-constructor-joins.json' and target.exists():
            f.require(json.loads(target.read_text()) == value,'unchanged complete native constructor joins')
            old_save(target.parent/'basis-recovery-joins.json',joins)
            return
        return old_save(target,value)
    f.atomic_pickle = retain_input;f.save = retain_native_joins
    f.require(native_main is original.main,'whole original main callable unchanged')
    native_main()


if __name__ == '__main__':main()
