#!/usr/bin/env python3
"""Finish saved row review after the scalar type boundary; no science replay."""
import ast
import copy
from pathlib import Path
import S11c_d_remaining_case_frequency_remainder_rows_review as original

NAME = 'S11c_d_remaining_case_frequency_remainder_rows_review_finish'
ORIGINAL_SHA = '750a519a144462a304e024287902679bfa2dc2f475b04c4e8770bcef79d23ae6'
PLAN_SHA = '58fd761310954af66d68023e6226354194c1b982ec735e658d8ffe8313380e45'
FAILED = original.F/'remainder-rows-review'
CURRENT = original.F/'remainder-rows-review-recovery-01'
INVENTORY_SHA = '6400956e3250c1812e9b35d4def8d00b620f68f5ea15f0c8a258f149e57a8a94'
OUTCOME_SHA = '52ed9e3ff26ee2f52486514f29d20ade210d76f58de1146cb9f1074bc83fff52'
RESTORED = ('c', 'arts', 'cp', 'owner', 'unit', 'partition', 'factors',
            'factor_slots', 'fourier_count', 'fourier_value_count',
            'fourier_slots', 'phase_owners', 'mixed_count', 'completed_pairs')


def restore(reader, journal, join):
    require, saved = original.require, original.saved
    old = reader.retain(original.__file__, ORIGINAL_SHA)
    reader.retain(original.M/'S11c_d_remaining_case_frequency_remainder_rows_review_plan.md', PLAN_SHA)
    inventory = reader.json(CURRENT/'failed-file-inventory.json', INVENTORY_SHA)
    for path, rec in inventory.items():
        got = reader.retain(path, rec['sha256'])
        require(got['canonical'] == rec['canonical'] and got['bytes'] == rec['bytes'], 'unchanged failed review file')
    outcome = reader.json(CURRENT/'failed-review-outcome.json', OUTCOME_SHA)
    require(outcome['guard']['exitCode'] == outcome['native']['exitCode'] == 1 and
            outcome['guard']['childOutcome']['guardReason'] is None and
            outcome['completedSummaries'] == {'factors': 8, 'fourier': 132, 'mixed': 82}, 'exact failed review boundary')
    trace = (FAILED/'validate.stderr').read_text()
    require('line 118, in main' in trace and "array(packet(kernel['constantProduct']),())" in trace and
            trace.endswith('AssertionError: full finite saved array/type/shape\n'), 'actual first unfinished comparison trace')
    c = reader.json(original.ROOT/'complete/checks.json', original.CHECKS_SHA)
    arts = c['artifacts']
    def packet(rec): return reader.packet(rec.get('path', rec.get('logical')), rec['sha256'])
    inputs = reader.json(arts['inputs.json']['path'], arts['inputs.json']['sha256'])
    for path, rec in inputs['consumedRoutes'].items(): reader.retain(path, rec['sha256'])
    # The completed original producer-outcome checks are source/trace joined;
    # keep the exact final producer logs in the continuation posthash inventory.
    for name in ('resource-guard/outcome.json', 'resource-guard/child-outcome.json',
                 'resource-guard/effective-limits.json', 'resource-guard/limit-validation.json',
                 'resource-guard/resource-samples.jsonl', 'resource-guard/stderr',
                 'frequency_remainder_rows.invocation.json', 'frequency_remainder_rows.stdout',
                 'frequency_remainder_rows.stderr', 'guard.stdout', 'guard.stderr'):
        reader.retain(original.ROOT/name)
    cp = reader.json(inputs['acceptedPreparation']['logical'], inputs['acceptedPreparation']['sha256'])
    owner = packet(next(v['ownInput'] for v in cp['sourceActions'] if v['rowIndex'] == 51))
    unit = owner['row']['factors'][0]['unit']
    k, q = (lim[0] for lim in owner['row']['limits'])
    partition = packet(arts['literal-factor-partition.pickle'])
    summaries = {}
    for group in ('factors', 'fourier', 'mixed'):
        records = []
        for path in sorted((FAILED/'complete'/group).glob('*.json'), key=lambda p: int(p.stem)):
            rec = reader.retain(path, inventory[str(path)]['sha256'])
            summary = reader.json(path, rec['sha256'])
            dest = journal.base/group/path.name
            dest.parent.mkdir(exist_ok=True)
            require(not dest.exists() and not dest.is_symlink(), 'fresh exact completed-summary reference')
            dest.symlink_to(Path(rec['canonical']))
            journal.artifacts[group+'/'+path.name] = {'path': str(dest), 'sha256': rec['sha256'], 'bytes': rec['bytes']}
            reader.retain(dest, rec['sha256'])
            records.append(summary)
        summaries[group] = records
    factor_slots = {}
    for summary in summaries['factors']:
        arg = packet(summary['input']); group = arg['group']; env = arg['environment']
        fi = partition['groups'][group].index(arg['expression'])
        coordinates = [None] if group == 'constant' else list(env[q] if group == 'input' else env[k][:, 0])
        for index, x in enumerate(coordinates):
            factor_slots[group, fi, None if x is None else float(x)] = {'packet': summary['value'], 'keys': [] if group == 'constant' else [index]}
    fourier_slots = {}; phase_owners = {}; fourier_value_count = 0
    for summary in summaries['fourier']:
        arg = packet(summary['input']); phase = summary['phase']
        if phase['disposition'] == 'NEW_PHASE': phase_owners[phase['input']['path']] = phase
        for index, x in enumerate(arg['frequencies']):
            fourier_slots[summary['sourceIndex'], float(x)] = {'packet': summary['value'], 'keys': [index]}
        fourier_value_count += summary['count']
    factors = len(summaries['factors']); fourier_count = len(summaries['fourier']); mixed_count = len(summaries['mixed'])
    completed_pairs = sum(s['shape'][0]*s['shape'][1] for s in summaries['mixed'])
    journal.json('completed-review-and-failure-reuse.json', {
        'originalReviewer': old, 'failedInventory': reader.retain(CURRENT/'failed-file-inventory.json', INVENTORY_SHA),
        'failedOutcome': reader.retain(CURRENT/'failed-review-outcome.json', OUTCOME_SHA),
        'trace': reader.retain(FAILED/'validate.stderr'), 'wholeMainJoin': join,
        'completedSummaryCounts': {key: len(value) for key, value in summaries.items()},
        'completedMixedPairs': completed_pairs, 'completedCoverageGuard': 'Literal original source and actual later traceback; no invented old receipt or coverage recomputation.',
        'restoredNames': list(RESTORED), 'restoration': 'Saved actual inputs and immutable completed metadata summaries only; no live cache inference.',
        'newScientificCalls': 0})
    values = locals()
    return tuple(values[name] for name in RESTORED)


def restore_coarse(meta, value, packet):
    label = '1024'; folder = 'grids/'+label
    summary = meta(folder+'/summary.json'); rule = packet(summary['rule'])
    mi = value(folder+'/mixed-input.pickle'); kernel = value(folder+'/kernel-input.pickle')
    return label, folder, summary, rule, mi, kernel


def adapted(source):
    """Stdlib-only exact reverse checks of completed-prefix/first-grid removal."""
    parsed = ast.parse(source)
    old = next(n for n in parsed.body if isinstance(n, ast.FunctionDef) and n.name == 'main')
    fn = copy.deepcopy(old)
    assert len(fn.body) == 102 and ast.unparse(fn.body[82]) == 'fine_rows = []'
    array_fn = fn.body[28]
    assert isinstance(array_fn, ast.FunctionDef) and array_fn.name == 'array'
    prior_test = copy.deepcopy(array_fn.body[0].value.args[0].values[0])
    assert ast.unparse(prior_test) == 'isinstance(v, np.ndarray)'
    array_fn.body[0].value.args[0].values[0] = ast.parse('(isinstance(v, np.ndarray) or (shape == () and isinstance(v, np.complexfloating)))', mode='eval').body
    loop = fn.body[83]
    boundary = next(i for i, n in enumerate(loop.body) if ast.unparse(n) == "array(packet(kernel['constantProduct']), ())")
    prefix = copy.deepcopy(loop.body[:boundary])
    replacement = ast.parse('if panels == 1024:\n    label, folder, summary, rule, mi, kernel = restore_coarse(meta, value, packet)\nelse:\n    pass').body[0]
    replacement.orelse = copy.deepcopy(prefix)
    loop.body[:boundary] = [replacement]
    count_guard = fn.body[92]
    count_nodes = [n for n in ast.walk(count_guard) if isinstance(n, ast.Call) and ast.unparse(n) == 'int(covered.sum())']
    assert len(count_nodes) == 1
    class Counts(ast.NodeTransformer):
        def visit_Call(self, n):
            return ast.Name(id='completed_pairs', ctx=ast.Load()) if ast.unparse(n) == 'int(covered.sum())' else self.generic_visit(n)
    fn.body[92] = Counts().visit(count_guard)
    # Reverse all suffix modifications before dropping completed statements.
    reverse = copy.deepcopy(fn)
    reverse.body[28].body[0].value.args[0].values[0] = prior_test
    reverse.body[83].body[:1] = prefix
    reverse.body[92] = copy.deepcopy(old.body[92])
    assert ast.dump(reverse) == ast.dump(old)
    assignment = ast.parse(', '.join(RESTORED)+' = restore(reader, journal, ADAPTER_JOIN)').body
    fn.body = fn.body[:10] + fn.body[25:29] + assignment + fn.body[58:64] + fn.body[82:]
    module = ast.fix_missing_locations(ast.Module(body=[fn], type_ignores=[]))
    compile(module, '<saved-row-review-finish>', 'exec')
    return module, {'wholeOriginalReverseAST': True, 'completedPrefixStatements': [10, 82],
                    'completedCoarseStatements': boundary, 'preservedSuffixBegins': "array(packet(kernel['constantProduct']), ())",
                    'arrayTypeChange': 'Allow native NumPy complex scalar only at shape(); no coercion.',
                    'completedSummaryReferences': 222, 'newScientificCalls': 0}


def main():
    original.require(original.saved.digest(original.__file__) == ORIGINAL_SHA, 'immutable original reviewer')
    module, join = adapted(Path(original.__file__).read_text())
    env = dict(vars(original), __file__=str(Path(__file__).resolve()), NAME=NAME,
               restore=restore, restore_coarse=restore_coarse, ADAPTER_JOIN=join)
    exec(compile(module, '<saved-row-review-finish>', 'exec'), env)
    env['main']()


if __name__ == '__main__': main()
