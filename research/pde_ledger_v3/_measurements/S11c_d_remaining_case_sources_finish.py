#!/usr/bin/env python3
"""Reuse 75 completed cell proofs and finish the last case plus source output."""
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
import S11c_d_remaining_case_sources_recover as recovery

f, engine, c, grades = recovery.f, recovery.engine, recovery.c, recovery.grades
PREVIOUS = recovery.ORIGIN.parent.parent/'recovery-01/complete'
PLAN = f.M/'S11c_d_remaining_case_sources_finish_plan.md'


def load(base):
    catalogue, manifest = recovery.load(base)
    previous_inputs = json.loads((PREVIOUS/'inputs.json').read_text())
    outcome = json.loads((PREVIOUS.parent/'cases_recover.invocation.json').read_text())
    f.require(outcome == json.loads((PREVIOUS.parent/'active.json').read_text()) and outcome['exitCode'] == 1,
              'actual previous recovery exit')
    f.require((PREVIOUS.parent/'cases_recover.stderr').read_text().rstrip().endswith('ValueError: actual restored census replay'),
              'previous session-dependent census guard failure')
    for name, sha in previous_inputs['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name) == f.digest(PREVIOUS/'source'/name) == sha, ('unchanged recovery source', name))
    for name, sha in previous_inputs['inputPackets'].items():
        f.require(f.digest(Path(name)) == sha, ('previous recovery input', name))
    for name, sha in previous_inputs['copiedPackets'].items():
        f.require(f.digest(PREVIOUS/name) == f.digest(base/name) == sha, 'same complete case source packets')
    expected = {'validation-'+'__'.join(case)+'.json' for case in c.CASES[:-1]}
    f.require({p.name for p in PREVIOUS.glob('validation-*.json')} == expected, 'three completed full case validations')
    copied = {}
    for path in sorted((PREVIOUS/'restored-cell-proofs').rglob('*.pickle')):
        relative = path.relative_to(PREVIOUS)
        destination = base/relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, destination)
        sha = f.digest(path)
        f.require(f.digest(destination) == sha, 'byte-identical completed proof copy')
        copied[str(relative)] = sha
    f.require(len(copied) == 75, 'exactly 75 completed cell proofs')
    paths = [PREVIOUS/'inputs.json', PREVIOUS.parent/'cases_recover.invocation.json',
             PREVIOUS.parent/'cases_recover.stderr', *PREVIOUS.glob('validation-*.json')]
    manifest['inputPackets'].update({str(p): f.digest(p) for p in paths})
    manifest['inputPackets'].update({str(PREVIOUS/name): sha for name, sha in copied.items()})
    for path in (Path(__file__).resolve(), PLAN):
        name = str(path.relative_to(f.ROOT))
        manifest['sourceFiles'][name] = f.digest(path)
        destination = base/'source'/name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, destination)
    manifest.update(previousRecoveryOutcome=outcome, completedProofCopies=copied,
        scope='Reuse 75 complete cell proofs and their case/source checks; validate only the last 25 cells and remaining output. Representation counts are retained diagnostics, not cross-session semantic identities.')
    f.save(base/'inputs.json', manifest)
    return catalogue, manifest


def inherit_completed(base):
    summaries, cache = {}, {}
    for case in c.CASES[:-1]:
        label = '__'.join(case)
        location, prefix = recovery.case_location(base, case)
        saved = f.unpickle(location/(prefix+'reduced-action.pickle'))
        actions = f.unpickle(location/(prefix+'actions.pickle'))
        assembled = f.unpickle(location/(prefix+'assembly.pickle'))
        summary = json.loads((PREVIOUS/('validation-'+label+'.json')).read_text())
        f.require(summary['cellsCertified'] == summary['normalizedResidualScalars'] == 25
                  and summary['respondingCoefficientMutations'] == 25 and summary['wrongAddressRejected'],
                  'completed validation guard summary')
        count = 0
        for row in assembled['result']['ROWS']:
            path = base/'restored-cell-proofs'/label/f"cell-{row['COLUMN']}-{row['ROW']}.pickle"
            proof = f.unpickle(path)
            f.require((proof['column'], proof['row']) == (row['COLUMN'], row['ROW']) and
                      proof['generators'] == row['GENERATORS'] and proof['coefficients'] == row['COEFFICIENTS'],
                      'actual saved proof/cell/generator/coefficient identity')
            actual_formal = engine.memo_xreplace(actions['columns'][row['COLUMN']][row['ROW']],
                                                dict(zip(proof['generators'], proof['symbols'])))
            f.require(sp.expand(actual_formal-proof['formal']) == 0, 'completed proof applies to actual saved column')
            f.require(proof['normalizedResidual'] == 0 and all(v == 0 for v in proof['coefficientResiduals'])
                      and proof['mutationIndex'] is not None and proof['mutationResidual'] != 0,
                      'completed exact proof and actual mutation preserved')
            f.require(proof['rowUnit'] == actions['columnUnits'][(row['COLUMN'], row['ROW'])], 'saved physical row unit')
            count += len(proof['coefficientResiduals'])
        f.require(count == summary['coefficientResidualScalars'], 'every completed coefficient proof scalar')
        f.require(summary['liveCensus'] == actions['census'], 'original live census remains immutable')
        summaries[label] = summary
        cache[case] = {'actions': actions, 'payloads': saved['payloads'], 'summary': summary}
        shutil.copyfile(PREVIOUS/('validation-'+label+'.json'), base/('validation-'+label+'.json'))
    f.save(base/'completed-validation-reuse.json', {'cases': list(summaries), 'cellProofsReused': 75,
        'coefficientProofScalars': sum(v['coefficientResidualScalars'] for v in summaries.values()),
        'derivativesRecomputed': 0, 'sourceConstraintChecksInherited': True,
        'scope': 'Actual column/generator/coefficient/unit joins to completed proofs; no repeat coefficient derivatives or prior case construction/constraint validation.'})
    return summaries, cache


def record_census_variant(base, label, actions, restored, diagnostic):
    # Persist the actual in-process forms before any remaining census guard.
    forms = []
    for j, column in enumerate(actions['columns']):
        for i, value in enumerate(column):
            text = sp.srepr(value)
            forms.append({'column': j, 'row': i, 'srepr': text,
                          'sha256': hashlib.sha256(text.encode()).hexdigest()})
    f.atomic_pickle(base/'last-case-current-census.pickle', {'case': label,
        'columns': actions['columns'], 'liveCensus': actions['census'],
        'diagnosticCensus': diagnostic['restoredCensus'], 'currentCensus': restored,
        'actualInProcessRepresentationStrings': forms})
    variation = recovery.census_changes(diagnostic['restoredCensus'], restored)
    f.save(base/'last-case-current-census.json', {'case': label, 'diagnosticToCurrentDifferences': variation,
        'liveCensus': actions['census'], 'diagnosticCensus': diagnostic['restoredCensus'], 'currentCensus': restored,
        'scope': 'Only DAG-node and distinct-derivative statistics may vary. All other census fields remain exact; complete semantic coefficient proofs are still required.'})


def remaining_validator():
    original = next(n for n in ast.parse(Path(recovery.__file__).read_text()).body
                    if getattr(n, 'name', None) == 'validate_cases')
    changed = copy.deepcopy(original)
    counts = {'cases': 0, 'snapshot': 0, 'prior': 0}
    saved_calls = {}
    class Forward(ast.NodeTransformer):
        def visit_For(self, node):
            self.generic_visit(node)
            if isinstance(node.target, ast.Name) and node.target.id == 'case' and ast.unparse(node.iter) == 'c.CASES':
                counts['cases'] += 1
                node.iter = ast.parse('(c.CASES[-1],)', mode='eval').body
            return node
        def visit_Call(self, node):
            self.generic_visit(node)
            if len(node.args) == 2 and isinstance(node.args[1], ast.Constant):
                label = node.args[1].value
                if label == 'actual restored census replay':
                    counts['snapshot'] += 1
                    saved_calls['snapshot'] = copy.deepcopy(node)
                    return ast.parse('record_census_variant(base, label, actions, restored, diagnostic)', mode='eval').body
                if label == 'only the recorded live/restored transitions':
                    counts['prior'] += 1
                    saved_calls['prior'] = copy.deepcopy(node)
                    node.args[0] = ast.parse("census_changes(actions['census'], diagnostic['restoredCensus']) == expected_differences[label]['differences']", mode='eval').body
            return node
    Forward().visit(changed)
    f.require(counts == {'cases': 1, 'snapshot': 1, 'prior': 1}, 'only last-case selection and two statistical comparison guards')
    class Reverse(ast.NodeTransformer):
        def visit_For(self, node):
            self.generic_visit(node)
            if isinstance(node.target, ast.Name) and node.target.id == 'case' and ast.unparse(node.iter) == '(c.CASES[-1],)':
                node.iter = ast.parse('c.CASES', mode='eval').body
            return node
        def visit_Call(self, node):
            self.generic_visit(node)
            if isinstance(node.func, ast.Name) and node.func.id == 'record_census_variant':
                return copy.deepcopy(saved_calls['snapshot'])
            if len(node.args) == 2 and isinstance(node.args[1], ast.Constant) and node.args[1].value == 'only the recorded live/restored transitions':
                return copy.deepcopy(saved_calls['prior'])
            return node
    reversed_body = Reverse().visit(copy.deepcopy(changed))
    f.require(ast.dump(reversed_body) == ast.dump(original), 'whole remaining-validator reverse AST join')
    namespace = dict(vars(recovery), record_census_variant=record_census_variant)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[changed], type_ignores=[])), str(Path(__file__)), 'exec'), namespace)
    return namespace['validate_cases'], {'wholeValidatorReverseAstJoin': True,
        'changes': counts, 'originalAstSha256': hashlib.sha256(ast.dump(original).encode()).hexdigest(),
        'derivedAstSha256': hashlib.sha256(ast.dump(changed).encode()).hexdigest()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('finish',), required=True)
    parser.add_argument('--run-directory', type=Path, required=True)
    args = parser.parse_args()
    base = args.run_directory.resolve()
    base.relative_to(f.STORE)
    base.mkdir(parents=True, exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    started = time.monotonic()
    def timeout(*_):
        raise TimeoutError('remaining case validation/output budget; preserve completed proofs')
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(900)
    catalogue, manifest = load(base)
    validator, validator_join = remaining_validator()
    validation, cache = inherit_completed(base)
    new_validation, new_cache = validator(base, catalogue, False)
    f.require(set(new_cache) == {c.CASES[-1]}, 'only remaining case validated')
    validation.update(new_validation)
    cache.update(new_cache)
    context, context_join = recovery.context_adapter(cache)
    f.save(base/'operand-validation.json', {'cases': validation, 'validatorJoin': validator_join,
        'contextJoin': context_join, 'cellProofsReused': 75, 'newCellProofs': 25,
        'scope': 'All actual source/coefficient semantics checked. Live and reconstructed expression statistics retained as diagnostics, not presumed session-invariant.'})
    old_context = c.baseline_context
    c.baseline_context = context
    try:
        metadata = c.output_and_replay(base, catalogue)
    finally:
        c.baseline_context = old_context
    original = {line.partition(': ')[0]: grades._restore(line.rstrip('\n').partition(': ')[2])
                for line in grades.decoded_lines(recovery.ORIGIN/'full.out')}
    current = {line.partition(': ')[0]: grades._restore(line.rstrip('\n').partition(': ')[2])
               for line in grades.decoded_lines(base/'full.out')}
    differences = [tag for tag, value in original.items() if current.get(tag) != value]
    f.save(base/'emission-differences.json', differences)
    f.require(not differences and set(original) <= set(current), 'every original decoded prefix payload identical')
    for name, sha in manifest['copiedPackets'].items():
        f.require(f.digest(recovery.ORIGIN/name) == f.digest(base/name) == sha, 'original/copied case packets unchanged')
    for name, sha in manifest['completedProofCopies'].items():
        f.require(f.digest(PREVIOUS/name) == f.digest(base/name) == sha, 'inherited proof packets unchanged')
    for name, sha in manifest['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name) == f.digest(base/'source'/name) == sha, 'current/frozen source hashes')
    for name, sha in manifest['inputPackets'].items():
        f.require(f.digest(Path(name)) == sha, 'all input hashes unchanged')
    checks = {'mode': args.mode, 'runDirectory': str(base), 'sourceFiles': manifest['sourceFiles'],
        'inputPackets': manifest['inputPackets'], 'copiedPackets': manifest['copiedPackets'],
        'completedProofCopies': manifest['completedProofCopies'], 'validatorJoin': validator_join,
        'contextJoin': context_join, 'cellProofsReused': 75, 'newCellProofs': 25,
        'cases': {label: {key: value for key, value in item.items() if key not in ('liveCensus', 'restoredCensus', 'differences')}
                  for label, item in validation.items()}, 'originalPrefixTags': len(original), **metadata,
        'artifacts': {str(p.relative_to(base)): {'bytes': p.stat().st_size, 'sha256': f.digest(p)}
                      for p in sorted(base.rglob('*')) if p.suffix in ('.pickle', '.out')
                      and 'source' not in p.relative_to(base).parts},
        'wallSeconds': time.monotonic()-started, 'peakRssKiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'scope': manifest['scope']}
    f.save(base/'recovery.json', {'sourceDirectory': str(recovery.ORIGIN), 'previousValidationDirectory': str(PREVIOUS),
        'validatorJoin': validator_join, 'contextJoin': context_join, 'originalPrefixTags': len(original),
        'cellProofsReused': 75, 'newCellProofs': 25})
    f.save(base/'checks.json', checks)
    signal.alarm(0)
    print(json.dumps(checks, indent=2))


if __name__ == '__main__':
    main()
