#!/usr/bin/env python3
"""Finish A9 navigation after a dynamic-class metadata label failed.

Keep the original reader and its two completed summaries unchanged. Only the
type-module display accepts None, and completed summary reads are skipped.
"""
import ast
import copy
import json
from pathlib import Path

import S11c_d_A9_saved_reader as native

HERE = Path(__file__).resolve()
FAILED = native.M.parents[2] / '_scratch/s11c/s11c-d-A9-inputs-20260923'
INSPECTION = FAILED / 'failure-inspection.json'
PLAN = native.M / 'S11c_d_A9_saved_reader_finish_plan.md'
COMPLETED = {}


def reuse_completed(base, index, canonical, consumers):
    if index not in COMPLETED:
        return None
    item = COMPLETED[index]
    source = Path(item['path'])
    if native.digest(source) != item['sha256']:
        raise ValueError('changed completed metadata')
    saved = json.loads(source.read_text())
    if saved['canonical'] != canonical or saved['consumers'] != consumers:
        raise ValueError('completed metadata full-input mismatch')
    target = base / source.name
    target.symlink_to(source.resolve(strict=True))
    return {'path': str(target), 'sha256': item['sha256'],
            'reusedSavedSummary': str(source)}


def patched_functions():
    original = Path(native.__file__).read_text()
    tree = ast.parse(original)
    functions = {n.name: n for n in tree.body if isinstance(n, ast.FunctionDef)}
    # Exact text edits affect metadata labels, never a physical input or value.
    display = 'type(value).__module__ +'
    fixed = original.replace(display, 'str(type(value).__module__) +')
    if original.count(display) != 2:
        raise ValueError('unexpected metadata type-label source')
    fixed_tree = ast.parse(fixed)
    for node in fixed_tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in ('label', 'describe'):
            exec(compile(ast.Module(body=[node], type_ignores=[]), str(HERE), 'exec'), vars(native))
    main = copy.deepcopy(functions['main'])
    insertion = ast.parse('''adopted = reuse_completed(base, index, canonical, consumers)
if adopted is not None:
    summaries.append(adopted)
    continue
''').body
    loops = [n for n in main.body if isinstance(n, ast.For)
             and ast.unparse(n.target) == '(index, (canonical, consumers))']
    if len(loops) != 1:
        raise ValueError('unexpected saved-reader loop')
    loop = loops[0]
    original_loop = copy.deepcopy(loop)
    loop.body[0:0] = insertion
    # Reverse the only loop change and compare the complete original body.
    reverse = copy.deepcopy(loop)
    del reverse.body[:len(insertion)]
    if ast.dump(reverse) != ast.dump(original_loop):
        raise ValueError('metadata continuation loop changed other work')
    source_list = next(n for n in main.body if isinstance(n, ast.Assign)
                       and ast.unparse(n.targets[0]) == 'source_paths')
    source_list.value.elts.extend([ast.Name(n, ast.Load()) for n in ('FINISH_HELPER', 'FINISH_PLAN', 'FAILURE_INSPECTION')])
    checks = next(n for n in main.body if isinstance(n, ast.Assign)
                  and ast.unparse(n.targets[0]) == 'checks').value
    index = next(i for i, k in enumerate(checks.keys) if k.value == 'nativeCodecRestorations')
    checks.values[index] = ast.parse('len(groups)-len(COMPLETED)', mode='eval').body
    checks.keys.append(ast.Constant('completedMetadataSummariesReused'))
    checks.values.append(ast.parse('len(COMPLETED)', mode='eval').body)
    native.reuse_completed = reuse_completed
    native.COMPLETED = COMPLETED
    native.FINISH_HELPER, native.FINISH_PLAN, native.FAILURE_INSPECTION = HERE, PLAN, INSPECTION
    exec(compile(ast.fix_missing_locations(ast.Module(body=[main], type_ignores=[])), str(HERE), 'exec'), vars(native))


def main():
    evidence = json.loads(INSPECTION.read_text())
    if native.digest(Path(native.__file__)) != evidence['originalHelperSha256']:
        raise ValueError('original failed helper changed')
    for item in evidence['completedSummaries']:
        COMPLETED[item['index']] = item
    if set(COMPLETED) != {0, 1}:
        raise ValueError('actual failed-run completed slots differ')
    patched_functions()
    native.main()


if __name__ == '__main__':
    main()
