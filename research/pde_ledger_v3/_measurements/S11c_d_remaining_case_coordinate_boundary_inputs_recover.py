#!/usr/bin/env python3
"""Resume input copying after two checked saved-packet schema corrections."""
import ast
import copy
import hashlib
import json
from pathlib import Path
import shutil
from types import SimpleNamespace

import S11c_d_remaining_case_coordinate_boundary_inputs as h

f = h.f
PLAN = f.M/'S11c_d_remaining_case_coordinate_boundary_recovery_plan.md'
REPAIR = f.M/'S11c_d_remaining_case_coordinate_boundary_schema_repair.json'
ORIGIN = f.STORE/'s11c-remaining-case-coordinate-20260921/boundary/focused'
ORIGINAL_RETAIN = h.modes.retain
ORIGINAL_MAIN = h.main


def adapt_tree(original):
    """Exactly one copy-enumeration expression and one census address change."""
    changed = copy.deepcopy(original)
    load = next(n for n in changed.body if isinstance(n, ast.FunctionDef) and n.name == 'load')
    loops = [n for n in ast.walk(load) if isinstance(n, ast.For)
             and ast.unparse(n.iter) == "('input-state', 'frame', 'profile-cache-pair')"]
    assert len(loops) == 1
    loops[0].iter = ast.parse("('input-state', 'frame') if label == BASELINE else ('input-state', 'frame', 'profile-cache-pair')", mode='eval').body
    prepare = next(n for n in changed.body if isinstance(n, ast.FunctionDef) and n.name == 'prepare')
    hits = [n for n in ast.walk(prepare) if isinstance(n, ast.Call) and ast.unparse(n) == "len(binding['binding']['rows'])"]
    assert len(hits) == 1
    hits[0].args[0] = ast.parse("binding['binding']['bound']['rows']", mode='eval').body
    reverse = copy.deepcopy(changed)
    load_back = next(n for n in reverse.body if isinstance(n, ast.FunctionDef) and n.name == 'load')
    loops_back = [n for n in ast.walk(load_back) if isinstance(n, ast.For) and isinstance(n.iter, ast.IfExp)
                  and ast.unparse(n.iter.test) == 'label == BASELINE']
    assert len(loops_back) == 1
    loops_back[0].iter = ast.parse("('input-state', 'frame', 'profile-cache-pair')", mode='eval').body
    prepare_back = next(n for n in reverse.body if isinstance(n, ast.FunctionDef) and n.name == 'prepare')
    hits_back = [n for n in ast.walk(prepare_back) if isinstance(n, ast.Call)
                 and ast.unparse(n) == "len(binding['binding']['bound']['rows'])"]
    assert len(hits_back) == 1
    hits_back[0].args[0] = ast.parse("binding['binding']['rows']", mode='eval').body
    assert ast.dump(reverse) == ast.dump(original)
    return changed, {'wholeFileReverseAst': True, 'copyEnumerationEdits': 1, 'rowCensusAddressEdits': 1,
        'originalAstSha256': hashlib.sha256(ast.dump(original).encode()).hexdigest(),
        'adaptedAstSha256': hashlib.sha256(ast.dump(changed).encode()).hexdigest()}


def recovered_load(base):
    repair = json.loads(REPAIR.read_text())
    f.require(f.digest(Path(h.__file__)) == repair['originalHelperSha256'], 'original helper remains byte-identical')
    f.require(not (ORIGIN/'complete/inputs.json').exists() and
              not (ORIGIN/'complete/preparation').exists(), 'original focus stopped during copying only')
    saved = {'runDirectory': str(base), 'inputPackets': {}, 'copiedInputs': {}}
    actual_files = {str(p.relative_to(ORIGIN/'complete')): f.digest(p)
                    for p in (ORIGIN/'complete').rglob('*') if p.is_file()}
    f.require(actual_files == repair['originalCompletedCopies'], 'exact entire saved partial load inventory')
    for name, sha in actual_files.items(): ORIGINAL_RETAIN(ORIGIN/'complete'/name, base/name, saved, sha)
    for name in repair['originalLogs']:
        ORIGINAL_RETAIN(ORIGIN/name, base/'original-focus-logs'/name, saved, repair['originalLogs'][name])
    ORIGINAL_RETAIN(Path(h.__file__), base/'original-input-source.py', saved, repair['originalHelperSha256'])
    existing_joins = {}
    def retain(src, dst, manifest, expected=None):
        relative = str(dst.relative_to(base))
        if dst.exists():
            f.require(relative in actual_files and f.digest(dst) == f.digest(src) == actual_files[relative]
                      and (expected is None or expected == actual_files[relative]), 'saved/requested complete input byte identity')
            manifest['inputPackets'][str(src)] = actual_files[relative]
            manifest['inputPackets'][str(ORIGIN/'complete'/relative)] = actual_files[relative]
            manifest['copiedInputs'][relative] = actual_files[relative]
            existing_joins[relative] = {'original': str(ORIGIN/'complete'/relative), 'requested': str(src),
                                        'sha256': actual_files[relative]}
        else: ORIGINAL_RETAIN(src, dst, manifest, expected)
    tree, join = adapt_tree(ast.parse(Path(h.__file__).read_text()))
    f.require(join == repair['astJoin'], 'reviewed exact schema adapter')
    load = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'load')
    namespace = dict(vars(h), modes=SimpleNamespace(retain=retain))
    exec(compile(ast.fix_missing_locations(ast.Module(body=[load], type_ignores=[])), __file__, 'exec'), namespace)
    manifest, labels = namespace['load'](base)
    f.require(set(existing_joins) == set(actual_files), 'every original completed copy reused in its actual request')
    for key in ('inputPackets', 'copiedInputs'):
        for name, sha in saved[key].items():
            f.require(name not in manifest[key] or manifest[key][name] == sha, 'same recovery provenance identity')
            manifest[key][name] = sha
    for path in (Path(__file__).resolve(), PLAN, REPAIR):
        name = str(path.relative_to(f.ROOT)); sha = f.digest(path); manifest['sourceFiles'][name] = sha
        target = base/'source'/name; target.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(path, target)
        f.require(f.digest(target) == sha, 'frozen recovery implementation')
    manifest['completedInputReuse'] = {'originalDirectory': str(ORIGIN/'complete'),
        'originalLoadManifestCompleted': False, 'copies': len(actual_files), 'joins': existing_joins,
        'astJoin': join, 'originalMainBytecodeUnchanged': h.main.__code__ is ORIGINAL_MAIN.__code__,
        'nativeScientificBodiesUnchanged': True, 'newScientificConstruction': False}
    f.save(base/'completed-input-reuse.json', manifest['completedInputReuse'])
    f.save(base/'inputs.json', manifest); h.matrices.hash_check(base, manifest)
    return manifest, labels


def main():
    tree, join = adapt_tree(ast.parse(Path(h.__file__).read_text()))
    repair = json.loads(REPAIR.read_text()); f.require(join == repair['astJoin'], 'exact reviewed full schema join')
    prepare = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'prepare')
    namespace = dict(vars(h))
    exec(compile(ast.fix_missing_locations(ast.Module(body=[prepare], type_ignores=[])), __file__, 'exec'), namespace)
    h.prepare = namespace['prepare']; h.load = recovered_load
    f.require(h.main.__code__ is ORIGINAL_MAIN.__code__, 'entire original main unchanged')
    h.main()


if __name__ == '__main__': main()
