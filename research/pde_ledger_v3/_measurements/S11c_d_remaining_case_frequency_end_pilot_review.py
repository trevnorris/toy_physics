#!/usr/bin/env python3
"""Bounded saved-result review; no evaluator, solver or map is invoked."""
import argparse
import hashlib
import json
from pathlib import Path
import pickle
import resource
import signal
import time

import numpy as np

M = Path(__file__).resolve().parent
REPO = M.parent.parent.parent
ROOT = REPO / '_scratch/s11c/s11c-remaining-case-frequency-20260921/end-pilot'
SOURCE = ROOT / 'complete'
CHECKS_SHA = '85dc47fc086ed1a25e3aae0dac954472d3dad042e8cccd14e7a2761476e5dd8e'
HELPER_SHA = '1d66333a4d11bba2a2198b7f93ae0ba03687659ef3a09a7f4a7d994b2e155e51'
PLAN_SHA = '0521c56777f9bce8ca8a41854502529c36bc4ae0971c4d6904a3893ef4bec54c'


def require(value, message):
    if not value:
        raise AssertionError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def same(a, b):
    if type(a) is not type(b):
        return False
    if isinstance(a, np.ndarray):
        return a.dtype == b.dtype and a.shape == b.shape and np.array_equal(a, b)
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(same(a[k], b[k]) for k in a)
    if isinstance(a, (list, tuple)):
        return len(a) == len(b) and all(same(x, y) for x, y in zip(a, b))
    return bool(a == b)


def magnitude(value):
    require(np.all(np.isfinite(value)), 'finite saved numerical values')
    return float(np.max(np.abs(value), initial=0.))


class Review:
    def __init__(self, base):
        self.base, self.routes, self.hashes = base, {}, {}

    def retain(self, path, expected=None):
        path = Path(path).absolute()
        real = path.resolve(strict=True)
        if str(real) not in self.hashes:
            self.hashes[str(real)] = digest(real)
        record = dict(logical=str(path), canonical=str(real),
                      directLink=str(path.readlink()) if path.is_symlink() else None,
                      bytes=real.stat().st_size, sha256=self.hashes[str(real)])
        require(expected is None or record['sha256'] == expected, ('hash', str(path)))
        require(str(path) not in self.routes or self.routes[str(path)] == record, 'stable saved route')
        self.routes[str(path)] = record
        return record

    def json(self, path, expected=None):
        record = self.retain(path, expected)
        return json.loads(Path(record['canonical']).read_text())

    def artifact(self, name):
        record = self.artifacts[name]
        require(record['path'] == str(SOURCE / name), 'literal producer artifact address')
        route = self.retain(SOURCE / name, record['sha256'])
        require(route['bytes'] == record['bytes'], 'artifact size')
        if 'canonicalValue' in record:
            require(route['canonical'] == record['canonicalValue'], 'immediate point byte reference')
        if name.endswith('.json'):
            return json.loads(Path(route['canonical']).read_text())
        with Path(route['canonical']).open('rb') as stream:
            return pickle.load(stream)

    def save(self, name, value):
        path = self.base / name
        require('..' not in Path(name).parts and not path.exists() and not path.is_symlink(), 'fresh review metadata')
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('x') as out:
            out.write(json.dumps(value, indent=2) + '\n')

    def postcheck(self):
        for path, expected in self.hashes.items():
            require(digest(path) == expected, ('post hash', path))
        for path, record in list(self.routes.items()):
            require(self.retain(path) == record, 'post link/size identity')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-directory', type=Path, required=True)
    base = parser.parse_args().run_directory.resolve()
    base.relative_to(REPO / '_scratch/s11c')
    base.mkdir(parents=True, exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS, (2 * 1024**3, 2 * 1024**3))
    signal.alarm(900)
    start = time.monotonic()
    r = Review(base)
    checks = r.json(SOURCE / 'checks.json', CHECKS_SHA)
    r.artifacts = checks['artifacts']
    outcome = r.json(ROOT / 'resource-guard/outcome.json')
    child = r.json(ROOT / 'resource-guard/child-outcome.json')
    invocation = r.json(ROOT / 'frequency_end_pilot.invocation.json')
    limits = r.json(ROOT / 'resource-guard/effective-limits.json')
    validation = r.json(ROOT / 'resource-guard/limit-validation.json')
    require(outcome['exitCode'] == child['exitCode'] == invocation['exitCode'] == 0 and
            child['guardReason'] is None and outcome['limitsVerified'] and validation['verified'], 'final clean producer guard')
    require(limits['memory.max'] == '2147483648' and limits['memory.swap.max'] == '0' and limits['pids.max'] == '32'
            and limits['nice'] == 15 and len(limits['affinity']) == 1 and set(limits['threads'].values()) == {'1'}, 'required containment')
    for name in ('frequency_end_pilot.stderr', 'guard.stderr', 'resource-guard/stderr'):
        path = ROOT / name
        require(r.retain(path)['bytes'] == 0, 'empty producer strict stderr')
    require(r.retain(ROOT / 'frequency_end_pilot.stdout')['sha256'] == CHECKS_SHA, 'stdout/checks byte identity')
    telemetry = ROOT / 'resource-guard/resource-samples.jsonl'
    r.retain(telemetry)
    samples = [json.loads(line) for line in telemetry.read_text().splitlines()]
    require(bool(samples), 'durable producer telemetry')
    for sample in samples:
        require(int(sample['memory.swap.current']) == 0 and int(sample['memory.peak']) <= 2*1024**3, 'zero swap and bounded memory')
        events = dict(line.split() for line in sample['memory.events'].splitlines())
        require(all(int(events[key]) == 0 for key in ('max', 'oom', 'oom_kill')), 'zero cap/OOM events')
    r.retain(M / 'S11c_d_remaining_case_frequency_end_pilot.py', HELPER_SHA)
    r.retain(M / 'S11c_d_remaining_case_frequency_end_pilot_plan.md', PLAN_SHA)
    r.retain(Path(__file__).resolve())
    r.retain(M / 'S11c_d_remaining_case_frequency_end_pilot_review_plan.md')
    r.json(ROOT / 'launch.json')
    inputs = r.artifact('inputs.json')
    for path, record in inputs['consumedRoutes'].items():
        require(r.retain(path, record['sha256']) == record, 'producer consumed source/input route')
    # Hash every new artifact, but deserialize only the bounded result evidence.
    for name, record in r.artifacts.items():
        route = r.retain(SOURCE / name, record['sha256'])
        require(route['bytes'] == record['bytes'], 'complete new artifact inventory')
    native = r.artifact('native-callers.json')
    require(all(native['adapter'][k] is True for k in ('equationReverseAST', 'solveReverseAST', 'continuationReverseAST',
            'exactFinalFrequencyAddress', 'nativeJacobianPackUnpackPolynomialMapsUnchanged')), 'recorded whole native source joins')
    scope = r.artifact('pilot-scope.json')
    require(scope['clusters'] == [0, 4, 6, 16, 17] and scope['directions'] == 7, 'whole bounded cluster scope')
    r.save('validated-source-and-scope.json', {'producerChecksSha256': CHECKS_SHA, 'scope': scope,
           'sourceFiles': {p: v for p, v in inputs['consumedRoutes'].items() if p.endswith('.py') or p.endswith('.md')},
           'producerGuard': outcome, 'artifactCount': len(r.artifacts)})

    v = r.artifact('evaluator-check/value.pickle')
    e = {
        'savedTangentRelativeDifference': magnitude(v['tangentDifference']) / (1 + magnitude(v['directionalDerivative'])),
        'finiteDifferenceResidual': magnitude(v['scaledFiniteDifferenceResidual']),
        'rationalValueRelativeDifference': magnitude(v['rationalValueDifference']) / (1 + magnitude(v['pencil'])),
        'rationalTangentRelativeDifference': magnitude(v['rationalTangentDifference']) / (1 + magnitude(v['directionalDerivative'])),
        'waveFrechetResidual': magnitude(r.artifact('evaluator-check/wave-value.pickle')['scaledDifference']),
    }
    require(e['savedTangentRelativeDifference'] < 1e-8 and e['rationalValueRelativeDifference'] < 1e-9
            and e['rationalTangentRelativeDifference'] < 1e-8 and e['finiteDifferenceResidual'] < 2e-6
            and e['waveFrechetResidual'] < 2e-6, 'saved new evaluator comparison evidence')
    changed = v['actualChangedCoefficientValue']
    original = v['rationalEntryValues']
    require(changed.shape == original.shape == (5, 5) and changed.dtype == original.dtype and
            changed[0, 0] == original[0, 0] + 1 and all(changed[i, j] == original[i, j]
            for i in range(5) for j in range(5) if (i, j) != (0, 0)), 'actual numerical coefficient mutation and unaffected entries')
    e['coefficientMutationResponded'] = True
    require(e == checks['evaluatorCheck'], 'evaluator summary derives saved residuals')
    r.save('validated-evaluator.json', e)
    frequencies = []
    for number in range(checks['coefficientFrequencyEvaluations']):
        arg = r.artifact(f'coefficient-evaluations/{number}/input.pickle')
        val = r.artifact(f'coefficient-evaluations/{number}/value.pickle')
        require(arg['frequency'] not in frequencies and arg['table'] == scope['source'], 'actual unique frequency/full table cache')
        frequencies.append(arg['frequency'])
        require(len(val[0]) == 25 and {tuple(t[:2]) for t in val[0]} == {(i, j) for i in range(5) for j in range(5)}, 'all actual rational entry addresses')
        require(val[1:] == (None, None), 'explicit separately routed acoustic evaluator slots')
        for row in val[0]:
            for terms in row[2:]:
                for _, coefficient in terms:
                    magnitude(coefficient)
    require(set(frequencies) == {1-.0075j, 1+0j, 1-.01j, 1-.005j}, 'fixed pilot frequency set')
    clusters, total_guards = {}, 0
    for index in scope['clusters']:
        folder = f'clusters/{index}'
        seed_input = r.artifact(folder + '/seed-input.pickle')
        seed = seed_input['nativeSeed']
        descriptor = seed_input['inputDescriptor']
        require(seed_input['source'] == scope['source'] and descriptor['savedWholeBasisKeys'] and
                not any(v['completeSavedCallerMatch'] for v in descriptor['baselineCandidates']), 'new complete own seed route')
        require(seed['index'] == index and all(same(seed['R'], seed['originalMode'][key])
                for key in descriptor['savedWholeBasisKeys']), 'actual saved whole seed basis')
        n = seed['R'].shape[1]
        size = 5*n + 2*n*n
        require(seed['R'].shape == (5, n) and n == (2 if index in (16, 17) else 1), 'complete cluster dimension')
        prepared = r.artifact(folder + '/prepared.pickle')
        gauge = r.artifact(folder + '/gauge-value.pickle')
        require(same(gauge['value'], prepared['gauge']), 'actual gauge value route')
        states = []
        guard_count = 0
        for number, frequency in enumerate((1+0j, 1-.01j, 1-.005j, 1-.01j)):
            prefix = folder + f'/solves/{number}'
            arg = r.artifact(prefix + '/input.pickle')
            state = r.artifact(prefix + '/value.pickle')
            receipt = r.artifact(prefix + '/completed.json')
            require(receipt['value'] == r.artifacts[prefix + '/value.pickle'], 'actual immediate state receipt')
            require(arg['frequency'] == state['frequency'] == frequency and arg['table'] == scope['source'] and
                    arg['seed'] == r.artifacts[folder + '/seed-input.pickle'], 'whole numerical caller join')
            shapes = {'x': (size,), 'R': (5, n), 'K': (n, n), 'Q': (n, n), 'tangent': (size,),
                      'residual': (size,), 'commutator': (n, n), 'jacobian': (size, size)}
            for key, shape in shapes.items():
                require(state[key].shape == shape, ('state shape', key))
                magnitude(state[key])
            for field, key in [('fixedGauge', 'gauge'), ('unknownScale', 'unknownScale'), ('rowScale', 'rowScale'), ('waveScale', 'waveScale')]:
                require(same(state[field], prepared[key]), 'fixed native preparation')
            require(same(state['x'][:5*n].reshape(5, n), state['R']) and
                    same(state['x'][5*n:5*n+n*n].reshape(n, n), state['K']) and
                    same(state['x'][5*n+n*n:].reshape(n, n), state['Q']), 'stored state coordinate packing')
            require(magnitude(state['residual']) < 2e-12 and magnitude(state['commutator']) < 1e-9, 'saved completed state residuals')
            last = state['history'][-1]
            require(last['residual'] == magnitude(state['residual']) and last['scaledJacobianCondition'] < 1e10
                    and last['minimumDenominatorSingularValue'] > 1e-12, 'saved regular-state history')
            for g in range(receipt['guardReceipts']):
                guard = r.artifact(prefix + f'/guard-{g}.pickle')
                require(bool(guard['condition']), 'actual native guard response')
            guard_count += receipt['guardReceipts']
            if number == 0:
                require(same(arg['initial'], state['x']) and same(state['R'], seed['R']) and
                        magnitude(r.artifact(folder + '/seed-unchanged.pickle')) == 0, 'accepted basis not corrected')
            else:
                link = folder + ('/coarse-path/point-1.pickle' if number == 1 else f'/fine-path/point-{number-1}.pickle')
                require(r.retain(SOURCE / link)['canonical'] == str((SOURCE / prefix / 'value.pickle').resolve()), 'actual saved path state byte reference')
            states.append(state)
        comparison = r.artifact(folder + '/comparison.pickle')
        require(same(comparison['coarse'], states[1]) and same(comparison['fine'], states[3]), 'full saved comparison states')
        difference = magnitude(comparison['scaledDifference'])
        require(difference < 1e-6, 'saved selected step comparison')
        clusters[index] = (seed, states)
        total_guards += guard_count
        r.save(f'validated-cluster-{index}.json', {'index': index, 'directions': n, 'states': 4,
               'nonseedPoints': 3, 'nativeGuardReceipts': guard_count, 'stepDifference': difference,
               'endResidual': magnitude(states[3]['residual']), 'condition': states[3]['history'][-1]['scaledJacobianCondition']})
    map_values = {}
    for label, state_number in [('coarse', 1), ('fine', 3)]:
        arg = r.artifact(f'maps/{label}-input.pickle')
        result = r.artifact(f'maps/{label}-value.pickle')
        require(arg['xposition'] == 64. and arg['source'] == scope['source'] and len(arg['clusters']) == 5, 'own map caller')
        for entry, (seed, states) in zip(arg['clusters'], clusters.values()):
            require(same(entry['seed'], seed) and same(entry['state'], states[state_number]), 'map consumes actual completed full cluster')
        for name in ('right', 'trace', 'extraction', 'observationAtCommonOrigin', 'inverseResidual', 'traceResidual'):
            require(result[name].shape == (5, 5), 'full outgoing map dimensions')
            magnitude(result[name])
        for name in ('incomingRight', 'forcingAtCommonOrigin', 'directIncomingSubtraction'):
            require(result[name].shape == (5, 2), 'complete incoming doublet map dimensions')
            magnitude(result[name])
        require(magnitude(result['inverseResidual']) < 1e-10 and magnitude(result['traceResidual']) < 1e-10
                and np.isfinite(result['condition']) and result['condition'] < 1e10, 'saved end-map residuals/conditioning')
        map_values[label] = result
        r.save(f'validated-map-{label}.json', {'condition': result['condition'], 'inverseResidual': magnitude(result['inverseResidual']),
               'traceResidual': magnitude(result['traceResidual']), 'receipt': r.artifact(f'maps/{label}-completed.json')})
    comparison = r.artifact('maps/comparison.pickle')
    relative = {k: magnitude(v) / (1 + magnitude(map_values['fine'][k])) for k, v in comparison['differences'].items()}
    require(relative == comparison['relative'] == checks['mapRelativeDifferences'] and max(relative.values()) < 1e-4, 'saved map comparison summaries')
    r.postcheck()
    r.save('validated-paths.json', r.routes)
    result = {'status': 'PASSED_BOUNDED_SAVED_END_PILOT_REVIEW', 'producerChecksSha256': CHECKS_SHA,
              'clusters': len(clusters), 'directions': sum(v[0]['R'].shape[1] for v in clusters.values()),
              'states': 20, 'nonseedPoints': 15, 'nativeGuardReceipts': total_guards, 'newScientificCalls': 0,
              'coefficientFrequencies': [repr(v) for v in frequencies], 'evaluatorCheck': e, 'mapRelativeDifferences': relative,
              'consumedLogicalPaths': len(r.routes), 'consumedPhysicalFiles': len(r.hashes), 'allConsumedHashesUnchanged': True,
              'wallSeconds': time.monotonic() - start,
              'scope': 'Saved bounded toy-model continuation at 1-0.01i; no pole, global domain or full scattering acceptance.'}
    require(result['clusters'] == checks['clusters'] and result['directions'] == checks['directions'] and
            result['states'] == checks['newCorrectedStates'] and result['nonseedPoints'] == checks['newNonseedPathPoints'], 'actual pilot census')
    r.save('checks.json', result)
    signal.alarm(0)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
