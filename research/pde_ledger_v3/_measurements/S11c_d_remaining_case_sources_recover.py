#!/usr/bin/env python3
"""Validate saved case assemblies, retaining live/restored expression censuses."""
import argparse
import ast
import copy
import hashlib
import json
from pathlib import Path
import resource
import shutil
import signal
import time
import sympy as sp
import S11c_d_remaining_case_sources as c

f, engine, grades = c.f, c.engine, c.grades
ORIGIN = f.STORE/'s11c-remaining-case-sources-20260919/production/complete'
DIAGNOSTIC = ORIGIN.parent.parent/'census-diagnostic'
PLAN = f.M/'S11c_d_remaining_case_sources_recovery_plan.md'
CATALOGUE_SHA = '23b98c40936f195664246baae6ae8fa69e645cca4bee60ac09775141bd060aa0'


def case_location(base, case):
    return (base, 'accepted-') if case == c.BASELINE else (base/'cases'/'__'.join(case), '')


def load(base):
    inputs = json.loads((ORIGIN/'inputs.json').read_text())
    outcome = json.loads((ORIGIN.parent/'cases_construct.invocation.json').read_text())
    f.require(outcome == json.loads((ORIGIN.parent/'active.json').read_text())
              and outcome['exitCode'] == 1 and outcome['status'] == 'failed', 'original final failure outcome')
    f.require((ORIGIN.parent/'cases_construct.stderr').read_text().rstrip().endswith('ValueError: accepted action census'),
              'original structural-census output failure')
    f.require(f.digest(ORIGIN/'remaining-case-sources.pickle') == CATALOGUE_SHA, 'complete saved case catalogue')
    catalogue = f.unpickle(ORIGIN/'remaining-case-sources.pickle')
    f.require(catalogue['sourceFiles'] == inputs['sourceFiles'] and catalogue['inputPackets'] == inputs['inputPackets'],
              'completed catalogue/source provenance')
    for name, sha in inputs['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name) == f.digest(ORIGIN/'source'/name) == sha, ('unchanged current/frozen source', name))
    for name, sha in inputs['inputPackets'].items():
        f.require(f.digest(Path(name)) == sha, ('unchanged accepted input', name))
    inventory = json.loads((ORIGIN/'packet-inventory-before-emission.json').read_text())
    f.require(set(inventory) == {str(p.relative_to(ORIGIN)) for p in ORIGIN.rglob('*.pickle')
                               if 'source' not in p.relative_to(ORIGIN).parts}, 'every completed original packet')
    for name, sha in inventory.items():
        original, destination = ORIGIN/name, base/name
        f.require(f.digest(original) == sha, ('original packet hash', name))
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(original, destination)
        f.require(f.digest(destination) == sha, ('byte-identical packet copy', name))
    for path in ORIGIN.glob('cases/*/checks.json'):
        destination = base/path.relative_to(ORIGIN)
        shutil.copyfile(path, destination)
    pins = dict(inputs['sourceFiles'])
    for path in (Path(__file__).resolve(), PLAN):
        pins[str(path.relative_to(f.ROOT))] = f.digest(path)
    operands = dict(inputs['inputPackets'])
    operands.update({str(ORIGIN/name): sha for name, sha in inventory.items()})
    paths = [ORIGIN/'inputs.json', ORIGIN/'full.out', ORIGIN/'packet-inventory-before-emission.json',
             ORIGIN.parent/'cases_construct.invocation.json', ORIGIN.parent/'cases_construct.stderr',
             DIAGNOSTIC/'differences.json', *DIAGNOSTIC.glob('*.pickle'), *ORIGIN.glob('cases/*/checks.json')]
    operands.update({str(path): f.digest(path) for path in paths})
    for name in pins:
        destination = base/'source'/name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(f.ROOT/name, destination)
    manifest = {'sourceFiles': pins, 'inputPackets': operands, 'originalInputs': inputs,
        'originalOutcome': outcome, 'copiedPackets': inventory,
        'scope': 'Saved case validation and output only. No reduction, probe-column construction, assembly extraction, integration or solve repeated.'}
    f.save(base/'inputs.json', manifest)
    return catalogue, manifest


def census_changes(live, restored):
    f.require(len(live) == len(restored) == 25, 'complete live/restored cell census')
    differences = []
    for before, after in zip(live, restored):
        f.require(set(before) == set(after), 'same census fields')
        changed = {key: {'live': before[key], 'restored': after[key]} for key in before if before[key] != after[key]}
        if changed:
            f.require(set(changed) <= {'DAG_NODES', 'DISTINCT_DERIVATIVES'}, 'only actual representation statistics may differ')
            differences.append({'column': before['COLUMN'], 'row': before['ROW'], 'changes': changed})
    return differences


def cell_certificate(actions, row, dimensions):
    j, i = row['COLUMN'], row['ROW']
    symbols = tuple(sp.Dummy('s11cdRestoredActionCarrier'+str(n)) for n in range(len(row['GENERATORS'])))
    mapping = dict(zip(row['GENERATORS'], symbols))
    formal = engine.memo_xreplace(actions['columns'][j][i], mapping)
    reconstructed = sum(coefficient*symbol for coefficient, symbol in zip(row['COEFFICIENTS'], symbols))
    raw = formal-reconstructed
    zero = dict.fromkeys(symbols, sp.S.Zero)
    residual = sp.expand(raw)
    derivatives = tuple(sp.expand(sp.diff(formal, symbol).xreplace(zero)-coefficient)
                        for symbol, coefficient in zip(symbols, row['COEFFICIENTS']))
    selected = next((n for n, coefficient in enumerate(row['COEFFICIENTS']) if coefficient != 0), None)
    mutation = None if selected is None else sp.expand(
        sp.diff(formal-reconstructed-row['COEFFICIENTS'][selected]*symbols[selected], symbols[selected]).xreplace(zero))
    return {'column': j, 'row': i, 'symbols': symbols, 'generators': row['GENERATORS'],
        'coefficients': row['COEFFICIENTS'], 'formal': formal, 'reconstructed': reconstructed,
        'rawResidual': raw, 'normalizedResidual': residual, 'coefficientResiduals': derivatives,
        'mutationIndex': selected, 'mutationResidual': mutation,
        'rowUnit': actions['columnUnits'][(j, i)],
        'generatorUnits': tuple(dimensions.measure(v) for v in row['GENERATORS'])}


def validate_cases(base, catalogue, focused=False):
    expected_differences = json.loads((DIAGNOSTIC/'differences.json').read_text())
    source = f.unpickle(base/'accepted-reduced-action.pickle')
    baseline = f.unpickle(base/'accepted-assembly.pickle')
    summaries, cache = {}, {}
    for case in c.CASES:
        label = '__'.join(case)
        location, prefix = case_location(base, case)
        saved = f.unpickle(location/(prefix+'reduced-action.pickle'))
        actions = f.unpickle(location/(prefix+'actions.pickle'))
        assembled = f.unpickle(location/(prefix+'assembly.pickle'))
        r, dimensions = c.source.restore_context(saved)
        dimensions.__dict__.update(assembled['dimensionState'])
        values = {key: engine.named(payload, 'VALUE') for (key, _), payload in saved['payloads'].items()}
        pencil = engine.ReducedPencil(*(values[key] for key in engine.CLOSED_KEYS), r)
        f.require(pencil.strong == actions['strong'] and pencil.kernel == actions['kernel']
                  and pencil.fields == actions['fields'] and pencil.probes == actions['probes'], 'complete source/action/field identity')
        diagnostic = f.unpickle(DIAGNOSTIC/(label+'.pickle'))
        f.require(diagnostic['columns'] == actions['columns'] and diagnostic['assembly'] == assembled['result']
                  and diagnostic['sourceStrong'] == diagnostic['actionStrong'] == actions['strong']
                  and diagnostic['savedCensus'] == actions['census'], 'actual diagnostic/source packet identities')
        restored = c.source.action_census(actions['columns'], pencil)
        f.require(restored == diagnostic['restoredCensus'], 'actual restored census replay')
        differences = census_changes(actions['census'], restored)
        f.require(differences == expected_differences[label]['differences'], 'only the recorded live/restored transitions')
        bad_census = copy.deepcopy(actions['census'])
        bad_census[0]['ROW'] += 1
        rejected = False
        try:
            census_changes(bad_census, restored)
        except ValueError:
            rejected = True
        f.require(rejected, 'changed physical census address rejected')
        c.validate_assembly(assembled['result'], actions)
        cells = [row for row in assembled['result']['ROWS'] if not focused or (row['COLUMN'], row['ROW']) == (0, 3)]
        certificates = []
        for row in cells:
            certificate = cell_certificate(actions, row, dimensions)
            destination = base/'restored-cell-proofs'/label
            destination.mkdir(parents=True, exist_ok=True)
            f.atomic_pickle(destination/f"cell-{row['COLUMN']}-{row['ROW']}.pickle", certificate)
            f.require(certificate['normalizedResidual'] == 0 and all(v == 0 for v in certificate['coefficientResiduals']),
                      ('saved-column/coefficient semantic reconstruction', label, row['COLUMN'], row['ROW']))
            f.require(certificate['mutationIndex'] is None or certificate['mutationResidual'] != 0,
                      'actual one-sided coefficient mutation responds')
            certificates.append(certificate)
        constraint_count = 0
        if case != c.BASELINE and not focused:
            for ordinal, key in enumerate(engine.CLOSED_KEYS):
                record = f.unpickle(location/f'reduction-{ordinal}.pickle')
                original = dict(source['rows'][key]['value'])[record['case']]
                payload = f.unpickle(location/f'payload-{ordinal}.pickle')
                f.require(record['key'] == key and tuple(map(str, record['case'])) == case
                          and record['sourcePayload'] == original and saved['rows'] == source['rows'], 'original case/source address')
                f.require(payload == saved['payloads'][(key, record['case'])] and
                          engine.named(payload, 'VALUE') == record['reduced'], 'full saved five-slot source identity')
                value = engine.named(original, 'VALUE')
                f.require(tuple(v[0] for v in record['records']) == tuple(sorted(value.atoms(sp.Integral), key=sp.default_sort_key)),
                          'all original integral occurrences')
                replay = engine.memo_xreplace(r.branches(r.strip_local(engine.memo_xreplace(
                    value, {v[0]: v[1] for v in record['records']}))), r.normal_map)
                f.require(replay == record['reduced'], 'literal saved reduction replay; no integral evaluated')
                for native, image, constraints in record['records']:
                    for equations, matrix, rhs, solutions, determinant in constraints:
                        solutions = list(solutions)
                        residual = matrix*sp.Matrix(solutions[0])-rhs
                        f.require(len(solutions) == 1 and determinant != 0 and matrix.det() == determinant
                                  and all(sp.expand(v) == 0 for v in residual), 'saved tangential constraint/Jacobian proof')
                        constraint_count += len(residual)+1
                f.require(saved['branchBindings'][(key, record['case'])] == tuple(
                    (eq.lhs, eq.rhs) for eq in engine.named(payload, 'COMPUTED_BRANCH_BINDINGS')), 'actual branch binding identity')
            original_summary = next(item for item in catalogue['cases'] if tuple(item['case']) == case)
            f.require(constraint_count == original_summary['constraintResidualScalars'], 'complete constraint proof census')
            comparisons = f.unpickle(location/'comparisons.pickle')
            current_matches = c.exact_integral_matches(assembled['result'], baseline['result'])
            f.require(current_matches == comparisons['integrals'], 'every original exact/new integral address replay')
            for order, difference in comparisons['localDifferencesAtCommonOrders'].items():
                f.require(difference == assembled['result']['LOCAL_MATRICES'][order]-baseline['result']['LOCAL_MATRICES'][order],
                          'actual complete local matrix differences')
        summary = {'liveCensus': actions['census'], 'restoredCensus': restored, 'differences': differences,
            'cellsCertified': len(certificates), 'normalizedResidualScalars': len(certificates),
            'coefficientResidualScalars': sum(len(v['coefficientResiduals']) for v in certificates),
            'respondingCoefficientMutations': sum(v['mutationIndex'] is not None for v in certificates),
            'rawNonzeroForms': sum(v['rawResidual'] != 0 for v in certificates),
            'constraintResidualScalars': constraint_count, 'wrongAddressRejected': rejected}
        f.save(base/('validation-'+label+'.json'), summary)
        cache[case] = {'actions': actions, 'payloads': saved['payloads'], 'summary': summary}
        summaries[label] = summary
    return summaries, cache


def context_adapter(cache):
    original = next(n for n in ast.parse(Path(c.__file__).read_text()).body if getattr(n, 'name', None) == 'baseline_context')
    changed = copy.deepcopy(original)
    guards = [n for n in ast.walk(changed) if isinstance(n, ast.Call) and len(n.args) == 2
              and isinstance(n.args[1], ast.Constant) and n.args[1].value == 'accepted action census']
    f.require(len(guards) == 1, 'one representation-sensitive guard')
    old_condition = copy.deepcopy(guards[0].args[0])
    guards[0].args[0] = ast.parse('checked_census(saved, actions, pencil)', mode='eval').body
    proof = copy.deepcopy(changed)
    next(n for n in ast.walk(proof) if isinstance(n, ast.Call) and len(n.args) == 2 and
         isinstance(n.args[1], ast.Constant) and n.args[1].value == 'accepted action census').args[0] = old_condition
    f.require(ast.dump(proof) == ast.dump(original), 'whole context reverse AST join; one guard only')
    def checked_census(saved, actions, pencil):
        cases = {tuple(map(str, case)) for key, case in saved['payloads']}
        f.require(len(cases) == 1, 'one complete case context')
        expected = cache[next(iter(cases))]
        f.require(saved['payloads'] == expected['payloads'] and actions == expected['actions'],
                  'exact already validated source/columns/census/units state')
        return True
    namespace = dict(vars(c), checked_census=checked_census)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[changed], type_ignores=[])), str(Path(__file__)), 'exec'), namespace)
    return namespace['baseline_context'], {'wholeContextReverseAstJoin': True,
        'originalAstSha256': hashlib.sha256(ast.dump(original).encode()).hexdigest(),
        'adapterAstSha256': hashlib.sha256(ast.dump(changed).encode()).hexdigest()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('focused', 'finish'), required=True)
    parser.add_argument('--run-directory', type=Path, required=True)
    args = parser.parse_args()
    base = args.run_directory.resolve()
    base.relative_to(f.STORE)
    base.mkdir(parents=True, exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    started = time.monotonic()
    def timeout(*_):
        raise TimeoutError('saved-case validation/output budget; preserve all proof packets')
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(900)
    catalogue, manifest = load(base)
    validation, cache = validate_cases(base, catalogue, args.mode == 'focused')
    adapter, join = context_adapter(cache)
    f.save(base/'operand-validation.json', {'cases': validation, 'contextJoin': join,
        'liveRepresentationStringsSaved': False,
        'scope': 'Original live representation statistics retained; restored full source/action/coefficient joins established independently. No equality of unavailable live expression trees is claimed.'})
    metadata, prefix_tags = {}, 0
    if args.mode == 'finish':
        old_context = c.baseline_context
        c.baseline_context = adapter
        try:
            metadata = c.output_and_replay(base, catalogue)
        finally:
            c.baseline_context = old_context
        original = {line.partition(': ')[0]: grades._restore(line.rstrip('\n').partition(': ')[2])
                    for line in grades.decoded_lines(ORIGIN/'full.out')}
        current = {line.partition(': ')[0]: grades._restore(line.rstrip('\n').partition(': ')[2])
                   for line in grades.decoded_lines(base/'full.out')}
        differences = [tag for tag, body in original.items() if current.get(tag) != body]
        f.save(base/'emission-differences.json', differences)
        f.require(not differences and set(original) <= set(current), 'every original decoded prefix payload identical')
        prefix_tags = len(original)
    for name, sha in manifest['copiedPackets'].items():
        f.require(f.digest(ORIGIN/name) == f.digest(base/name) == sha, 'all original/copied pre/post packets')
    for name, sha in manifest['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name) == f.digest(base/'source'/name) == sha, 'current/frozen sources unchanged')
    for name, sha in manifest['inputPackets'].items():
        f.require(f.digest(Path(name)) == sha, 'all inputs unchanged')
    checks = {'mode': args.mode, 'runDirectory': str(base), 'sourceFiles': manifest['sourceFiles'],
        'inputPackets': manifest['inputPackets'], 'contextJoin': join, 'copiedPackets': manifest['copiedPackets'],
        'cases': {case: {key: value for key, value in record.items() if key not in ('liveCensus', 'restoredCensus', 'differences')}
                  for case, record in validation.items()}, 'originalPrefixTags': prefix_tags, **metadata,
        'artifacts': {str(p.relative_to(base)): {'bytes': p.stat().st_size, 'sha256': f.digest(p)}
                      for p in sorted(base.rglob('*')) if p.suffix in ('.pickle', '.out')
                      and 'source' not in p.relative_to(base).parts},
        'wallSeconds': time.monotonic()-started, 'peakRssKiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'scope': manifest['scope']}
    f.save(base/'recovery.json', {'originalDirectory': str(ORIGIN), 'contextJoin': join,
        'copiedPackets': manifest['copiedPackets'], 'originalPrefixTags': prefix_tags, 'mode': args.mode})
    f.save(base/'checks.json', checks)
    signal.alarm(0)
    print(json.dumps(checks, indent=2))


if __name__ == '__main__':
    main()
