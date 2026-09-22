#!/usr/bin/env python3
"""Finish the saved-input reader, retaining its completed metadata atlas."""
import ast
import copy
import hashlib
import json
from pathlib import Path

import S11c_d_remaining_case_frequency_end_continuation_inputs as original

M = original.M
PREVIOUS = original.REPO / '_scratch/s11c/s11c-remaining-case-frequency-20260921/end-continuation-inputs'
PLAN = M / 'S11c_d_remaining_case_frequency_end_continuation_inputs_recovery_plan.md'
REPAIR = M / 'S11c_d_remaining_case_frequency_end_continuation_inputs_reader_repair.json'
ORIGINAL_SHA = '9192296c7d060e176c321b4c0db4af888a038027a364e655ec489208f42a238f'


def adapters():
    path = Path(original.__file__)
    original.require(original.digest(path) == ORIGINAL_SHA, 'immutable failed reader')
    module = ast.parse(path.read_text())
    prohibit = next(n for n in module.body if isinstance(n, ast.FunctionDef) and n.name == 'prohibit')
    patched = copy.deepcopy(prohibit)
    removed = [n for n in patched.body if isinstance(n, ast.For) and ast.unparse(n.iter) == "('diff', 'subs', 'xreplace')"]
    original.require(len(removed) == 1, 'single overbroad pickle-class restriction')
    slot = patched.body.index(removed[0])
    patched.body.pop(slot)
    reverse = copy.deepcopy(patched)
    reverse.body.insert(slot, copy.deepcopy(removed[0]))
    original.require(ast.dump(reverse) == ast.dump(prohibit), 'whole prohibit reversal')

    inspect = next(n for n in module.body if isinstance(n, ast.FunctionDef) and n.name == 'inspect_inputs')
    finished = copy.deepcopy(inspect)
    start = next(i for i, n in enumerate(finished.body) if isinstance(n, ast.Assign) and ast.unparse(n.targets[0]) == 'binding_routes')
    stop = next(i for i, n in enumerate(finished.body) if isinstance(n, ast.Expr) and isinstance(n.value, ast.Call)
                and ast.unparse(n.value.func) == 'save' and any(isinstance(a, ast.Constant) and a.value == 'saved-binding-and-tangent-routes.json' for a in n.value.args)) + 1
    prefix = copy.deepcopy(finished.body[start:stop])
    replacement = ast.parse('restore_binding_routes(reader, base)').body
    finished.body[start:stop] = replacement
    reverse = copy.deepcopy(finished)
    reverse.body[start:start + len(replacement)] = prefix
    original.require(ast.dump(reverse) == ast.dump(inspect), 'whole inspection reversal; completed binding reads skipped')
    compiled = ast.fix_missing_locations(ast.Module(body=[patched, finished], type_ignores=[]))
    env = dict(original.__dict__, restore_binding_routes=restore_binding_routes)
    exec(compile(compiled, '<saved-input codec and completed-metadata continuation>', 'exec'), env)
    return env['prohibit'], env['inspect_inputs'], {'wholeProhibitReverseAST': True,
            'wholeInspectionReverseAST': True, 'removedOnlyBasicMethodPatching': True,
            'originalMainAndComparisonBodiesUnchanged': True, 'completedBindingMetadataRestored': True}


def saved(reader, name):
    repair = reader.json(REPAIR)
    path = PREVIOUS / 'complete' / name
    return reader.json(path, repair['completedMetadata'][name]['sha256'])


def restore_routes(reader, value):
    if isinstance(value, dict):
        if {'logical', 'canonical', 'sha256', 'bytes'} <= value.keys():
            current = reader.retain(value['logical'], value['sha256'])
            original.require(current == value, 'exact completed logical/canonical metadata route')
        else:
            for item in value.values():
                restore_routes(reader, item)
    elif isinstance(value, list):
        for item in value:
            restore_routes(reader, item)


def restore_atlas(reader, base):
    result = {}
    for tag in ('focused', 'matrix', 'contour16', 'contour32'):
        result[tag] = saved(reader, 'baseline/' + tag + '.json')
        restore_routes(reader, result[tag])
    original.save(base, 'completed-baseline-atlas-reuse.json', {
        tag: {'metadata': reader.retain(PREVIOUS / 'complete/baseline' / (tag + '.json')),
              'mapFiles': len(value['savedMaps'])} for tag, value in result.items()})
    return result


def restore_binding_routes(reader, base):
    value = saved(reader, 'saved-binding-and-tangent-routes.json')
    restore_routes(reader, value)
    original.save(base, 'completed-binding-route-reuse.json', {
        'metadata': reader.retain(PREVIOUS / 'complete/saved-binding-and-tangent-routes.json'),
        'restoredBindingRoutes': list(value), 'newBindingOrDerivativeCalls': 0})


def main():
    repair = json.loads(REPAIR.read_text())
    prohibit, inspect, joins = adapters()
    old_join = original.source_join

    def source_join(reader, cp, origin):
        result = old_join(reader, cp, origin)
        for path in (Path(__file__).resolve(), PLAN, REPAIR):
            reader.retain(path)
        for name, item in repair['preservedFiles'].items():
            reader.retain(PREVIOUS / name, item['sha256'])
        result['readerContinuation'] = {'joins': joins, 'repair': reader.retain(REPAIR),
                                       'originalFailure': repair['failure'],
                                       'pickleCodecMethodsRemainNative': True,
                                       'scientificApplicationCallsStillDisabled': True}
        return result

    original.source_join = source_join
    original.prohibit = prohibit
    original.numerical_atlas = restore_atlas
    original.inspect_inputs = inspect
    original.main()


if __name__ == '__main__':
    main()
