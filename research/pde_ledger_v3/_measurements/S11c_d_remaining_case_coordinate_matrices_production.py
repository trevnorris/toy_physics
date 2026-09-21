#!/usr/bin/env python3
"""Preserve the accepted focus inventory while adding production case results."""
import ast
import copy
import inspect
import json
from pathlib import Path
import shutil

import S11c_d_remaining_case_coordinate_matrices as h

f = h.f
PLAN = f.M/'S11c_d_remaining_case_coordinate_matrices_production_plan.md'
ORIGINAL_LOAD = h.load
ORIGINAL_CONSTRUCT = h.construct


def production_construct():
    original = ast.parse(inspect.getsource(ORIGINAL_CONSTRUCT)).body[0]
    changed = copy.deepcopy(original)
    calls = [n for n in ast.walk(changed) if isinstance(n, ast.Call)
             and ast.unparse(n.func) == 'f.save' and n.args
             and ast.unparse(n.args[0]) == "base / 'case-inventory.json'"]
    f.require(len(calls) == 1, 'one growing production inventory write')
    calls[0].args[0].right.value = 'production-case-inventory.json'
    reverse = copy.deepcopy(changed)
    reverse_calls = [n for n in ast.walk(reverse) if isinstance(n, ast.Call)
                     and ast.unparse(n.func) == 'f.save' and n.args
                     and ast.unparse(n.args[0]) == "base / 'production-case-inventory.json'"]
    f.require(len(reverse_calls) == 1, 'one separate production inventory address')
    reverse_calls[0].args[0].right.value = 'case-inventory.json'
    f.require(ast.dump(reverse) == ast.dump(original), 'whole production coordinator reverse AST')
    namespace = dict(vars(h))
    exec(compile(ast.fix_missing_locations(ast.Module(body=[changed], type_ignores=[])), __file__, 'exec'), namespace)
    return namespace['construct'], {'wholeCoordinatorReverseAst': True,
        'changedConstants': 1, 'old': 'case-inventory.json', 'new': 'production-case-inventory.json',
        'originalCoordinatorAstSha256': h.body(ORIGINAL_CONSTRUCT),
        'originalMainAstSha256': h.body(h.main), 'nativeBodiesUnchanged': True}


def load(base, resume):
    f.require(resume is not None, 'accepted saved material matrix focus required')
    manifest, labels = ORIGINAL_LOAD(base, resume)
    _, join = production_construct()
    original_inventory = resume/'case-inventory.json'
    f.require(f.digest(base/'case-inventory.json') == f.digest(original_inventory), 'accepted focus inventory unchanged')
    for path in (Path(__file__).resolve(), PLAN):
        name = str(path.relative_to(f.ROOT)); sha = f.digest(path)
        f.require(name not in manifest['sourceFiles'], 'explicit new inventory wrapper provenance')
        manifest['sourceFiles'][name] = sha
        target = base/'source'/name; target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target); f.require(f.digest(target) == sha, 'frozen inventory wrapper')
    manifest['productionInventoryJoin'] = {**join, 'acceptedFocusInventorySha256': f.digest(original_inventory),
        'acceptedFocusInventory': str(original_inventory),
        'scope': 'Only the growing production inventory filename changes; every focused artifact remains byte-identical.'}
    f.save(base/'production-inventory-join.json', manifest['productionInventoryJoin'])
    f.save(base/'inputs.json', manifest)
    return manifest, labels


def main():
    h.construct, _ = production_construct()
    h.load = load
    h.main()


if __name__ == '__main__': main()
