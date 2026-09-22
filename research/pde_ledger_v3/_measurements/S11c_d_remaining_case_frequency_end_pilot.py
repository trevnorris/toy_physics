#!/usr/bin/env python3
"""One bounded numerical pilot for the distinct RHOBR right end.

The accepted rational expressions are evaluated numerically with forward dual
arithmetic. No symbolic derivative or lambdify is reconstructed. The existing
whole-cluster equations, fixed gauge, Newton solver and maps remain the native
algorithm; only coefficient evaluation and the acoustic expression reader are
new. Every completed point and comparison is durable before later guards.
"""
import argparse
import ast
import copy
import inspect
import json
from pathlib import Path
import pickle
import resource
import signal
import textwrap
import time

import numpy as np
import sympy as sp
import S11c_d_frequency_end as native
import S11c_d_remaining_case_frequency_end_continuation_inputs as saved

M, REPO = saved.M, saved.REPO
PLAN = M / 'S11c_d_remaining_case_frequency_end_pilot_plan.md'
READY = REPO / '_scratch/s11c/s11c-remaining-case-frequency-20260921/end-continuation-inputs-recovery-01'
READY_SHA = '7374e0187c8aefca5817e182c06bb7785ec7eb13ab069b66f815d6b246db0263'
READY_INVENTORY_SHA = 'bf7513d06f18e8271c59fecb51c0c27dfce1b24e67a51456ff06a7cbfa21f024'
OWNER = 'LAB_HELD__RHOBR_CONSTANT'
TARGET = 1. - .01j
require = saved.require
norm = native.norm


class Journal:
    def __init__(self, base):
        self.base = base
        self.artifacts = {}

    def write(self, name, value):
        path = self.base / name
        require('..' not in Path(name).parts and not path.exists() and not path.is_symlink(), 'fresh pilot output')
        parent = path.parent
        while parent != self.base:
            require(not parent.is_symlink(), 'no write through input reference')
            parent = parent.parent
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(path.suffix + '.new')
        with temporary.open('xb') as out:
            pickle.dump(value, out, protocol=pickle.HIGHEST_PROTOCOL)
        temporary.rename(path)
        record = {'path': str(path), 'sha256': saved.digest(path), 'bytes': path.stat().st_size}
        self.artifacts[name] = record
        return record

    def json(self, name, value):
        saved.save(self.base, name, value)
        path = self.base / name
        self.artifacts[name] = {'path': str(path), 'sha256': saved.digest(path), 'bytes': path.stat().st_size}


class NumericExpressions:
    """A numerical evaluator, not a producer of symbolic derivative operands."""
    def __init__(self, coordinates):
        self.coordinates = coordinates
        self.nodes = []
        self.memo = {}

    def add(self, expression):
        if expression in self.memo:
            return self.memo[expression]
        if expression in self.coordinates:
            node = ('variable', self.coordinates.index(expression))
        elif expression is sp.I:
            node = ('constant', 1j)
        elif isinstance(expression, sp.Rational):
            node = ('constant', int(expression.p) / int(expression.q))
        elif isinstance(expression, sp.Float):
            node = ('constant', float(expression))
        elif isinstance(expression, sp.Add):
            node = ('add', tuple(self.add(v) for v in expression.args))
        elif isinstance(expression, sp.Mul):
            node = ('multiply', tuple(self.add(v) for v in expression.args))
        elif isinstance(expression, sp.Pow) and isinstance(expression.exp, sp.Integer):
            node = ('power', self.add(expression.base), int(expression.exp))
        else:
            raise TypeError(('unsupported saved numerical expression', type(expression), repr(expression)))
        index = len(self.nodes)
        self.nodes.append(node)
        self.memo[expression] = index
        return index

    def scalar(self, roots, point, direction=(0., 0., 0.)):
        values = []
        derivatives = []
        for node in self.nodes:
            kind = node[0]
            if kind == 'constant':
                value, derivative = complex(node[1]), 0j
            elif kind == 'variable':
                value, derivative = point[node[1]], direction[node[1]]
            elif kind == 'add':
                value = sum(values[i] for i in node[1])
                derivative = sum(derivatives[i] for i in node[1])
            elif kind == 'multiply':
                value, derivative = 1 + 0j, 0j
                for i in node[1]:
                    derivative = derivative * values[i] + value * derivatives[i]
                    value *= values[i]
            else:
                i, power = node[1:]
                value = values[i] ** power
                derivative = 0j if power == 0 else power * values[i] ** (power - 1) * derivatives[i]
            values.append(value)
            derivatives.append(derivative)
        return [values[i] for i in roots], [derivatives[i] for i in roots]

    def matrix_wave(self, root, point, direction):
        # The saved acoustic expression is polynomial in the independent K,Q.
        # This directly evaluates its arithmetic; no old Poly/terms object is
        # recreated. Matrix products keep the literal expression order.
        w, K, Q = point
        dw, dK, dQ = direction
        I = np.eye(K.shape[0], dtype=complex)
        Z = np.zeros_like(I)
        variables = (w * I, K, Q)
        tangents = (dw * I, dK, dQ)
        values, derivatives = [], []
        for node in self.nodes:
            kind = node[0]
            if kind == 'constant':
                value, derivative = complex(node[1]) * I, Z
            elif kind == 'variable':
                value, derivative = variables[node[1]], tangents[node[1]]
            elif kind == 'add':
                value = sum((values[i] for i in node[1]), start=Z.copy())
                derivative = sum((derivatives[i] for i in node[1]), start=Z.copy())
            elif kind == 'multiply':
                value, derivative = I, Z
                for i in node[1]:
                    derivative = derivative @ values[i] + value @ derivatives[i]
                    value = value @ values[i]
            else:
                i, power = node[1:]
                require(power >= 0, 'polynomial saved acoustic wave only')
                value = np.linalg.matrix_power(values[i], power)
                derivative = sum((np.linalg.matrix_power(values[i], j) @ derivatives[i] @
                                  np.linalg.matrix_power(values[i], power - 1 - j)
                                  for j in range(power)), start=Z.copy())
            values.append(value)
            derivatives.append(derivative)
        return values[root], derivatives[root]


class Coefficients:
    def __init__(self, journal, table, source_route):
        self.journal, self.table, self.route = journal, table, source_route
        source = table['source']
        coordinates = (source['frequency'], source['momentum'], source['radical'])
        journal.write('numerical-evaluator/input.pickle', {'source': source_route, 'coordinates': coordinates,
                      'coefficientAddresses': [(r['row'], r['column']) for r in table['entries']],
                      'algorithm': 'numeric forward dual arithmetic; no symbolic diff/Poly/lambdify'})
        self.program = NumericExpressions(coordinates)
        self.rows = []
        for row in table['entries']:
            terms = row['numeratorTerms'] + row['denominatorTerms']
            require(all(not c.free_symbols - {coordinates[0]} for _, c in terms), 'own coefficients depend only on frequency')
            roots = [self.program.add(c) for _, c in terms]
            self.rows.append((row['row'], row['column'], len(row['numeratorTerms']), [p for p, _ in terms], roots))
        self.wave = NumericExpressions(coordinates)
        self.wave_root = self.wave.add(source['wave'])
        self.cache = {}
        journal.write('numerical-evaluator/program.pickle', {'nodes': self.program.nodes, 'rows': self.rows,
                      'waveNodes': self.wave.nodes, 'waveRoot': self.wave_root})

    def get(self, w):
        if w not in self.cache:
            number = len(self.cache)
            folder = 'coefficient-evaluations/' + str(number)
            self.journal.write(folder + '/input.pickle', {'frequency': w, 'table': self.route,
                               'direction': (1., 0., 0.), 'program': self.journal.artifacts['numerical-evaluator/program.pickle']})
            roots = [i for _, _, _, _, indices in self.rows for i in indices]
            values, derivatives = self.program.scalar(roots, (w, 0j, 0j), (1., 0., 0.))
            entries, offset = [], 0
            for i, j, m, powers, indices in self.rows:
                n = len(indices)
                v, d = values[offset:offset+n], derivatives[offset:offset+n]
                entries.append((i, j, list(zip(powers[:m], v[:m])), list(zip(powers[m:], v[m:])),
                                list(zip(powers[:m], d[:m])), list(zip(powers[m:], d[m:]))))
                offset += n
            self.cache[w] = (entries, None, None)
            self.journal.write(folder + '/value.pickle', self.cache[w])
        return self.cache[w]


def native_adapters():
    equation = ast.parse(textwrap.dedent(inspect.getsource(native.Pair.equation))).body[0]
    adapted = copy.deepcopy(equation)
    first = next(i for i, n in enumerate(adapted.body) if isinstance(n, ast.Assign) and
                 ast.unparse(n.value) == 'polynomial(wave, K, Q, dK, dQ)')
    original_wave = copy.deepcopy(adapted.body[first:first+2])
    require(isinstance(original_wave[1], ast.If) and ast.unparse(original_wave[1].body[0]) ==
            'dW += dw * polynomial(dwave, K, Q)[0]', 'exact native wave derivative slot')
    adapted.body[first:first+2] = ast.parse('W,dW=self.wave_equation(w,K,Q,dK,dQ,dw)').body
    reversed_equation = copy.deepcopy(adapted)
    reversed_equation.body[first:first+1] = original_wave
    require(ast.dump(reversed_equation) == ast.dump(equation), 'whole native rational invariant equation unchanged outside acoustic evaluator')
    solve = ast.parse(textwrap.dedent(inspect.getsource(native.Pair.solve))).body[0]
    class Observe(ast.NodeTransformer):
        def visit_Call(self, node):
            if ast.unparse(node.func) == 'f.require':
                return ast.copy_location(ast.Call(func=ast.Attribute(value=ast.Name(id='self', ctx=ast.Load()),
                       attr='check', ctx=ast.Load()), args=node.args + [ast.Call(func=ast.Name(id='locals', ctx=ast.Load()), args=[], keywords=[])], keywords=[]), node)
            return self.generic_visit(node)
    observed = Observe().visit(copy.deepcopy(solve))
    class Reverse(ast.NodeTransformer):
        def visit_Call(self, node):
            if ast.unparse(node.func) == 'self.check':
                return ast.copy_location(ast.Call(func=ast.Attribute(value=ast.Name(id='f', ctx=ast.Load()),
                       attr='require', ctx=ast.Load()), args=node.args[:-1], keywords=[]), node)
            return self.generic_visit(node)
    require(ast.dump(Reverse().visit(copy.deepcopy(observed))) == ast.dump(solve), 'whole native Newton/guards unchanged with persistence')
    continuation = ast.parse(textwrap.dedent(inspect.getsource(native.continue_pair))).body[0]
    continued = copy.deepcopy(continuation)
    continued.name = 'continued'
    continued.args.args.append(ast.arg(arg='label'))
    replacements = {
        'seedstate[\'frequency\'] + (target - seedstate[\'frequency\']) * j / count':
            "target if j == count else seedstate['frequency'] + (target - seedstate['frequency']) * j / count",
        'pair.solve(w, initial)': 'pair.corrected(w, initial, label)',
        'f.atomic_pickle(path, state)': 'pair.reference_point(path, state)',
    }
    counts = {key: 0 for key in replacements}
    class ContinueRoute(ast.NodeTransformer):
        def visit(self, node):
            if isinstance(node, ast.expr) and ast.unparse(node) in replacements:
                key = ast.unparse(node)
                counts[key] += 1
                return ast.copy_location(ast.parse(replacements[key], mode='eval').body, node)
            return super().visit(node)
    continued = ContinueRoute().visit(continued)
    require(set(counts.values()) == {1}, 'exact native continuation routes')
    inverse = {value: key for key, value in replacements.items()}
    class ContinueReverse(ast.NodeTransformer):
        def visit(self, node):
            if isinstance(node, ast.expr) and ast.unparse(node) in inverse:
                return ast.copy_location(ast.parse(inverse[ast.unparse(node)], mode='eval').body, node)
            return super().visit(node)
    reverse = ContinueReverse().visit(copy.deepcopy(continued))
    reverse.name = continuation.name
    reverse.args.args.pop()
    require(ast.dump(reverse) == ast.dump(continuation), 'whole native continuation reverse AST')
    module = ast.fix_missing_locations(ast.Module(body=[adapted, observed, continued], type_ignores=[]))
    env = dict(native.__dict__)
    exec(compile(module, '<native end equations with numerical coefficient reader>', 'exec'), env)
    return env['equation'], env['solve'], env['continued'], {'equationReverseAST': True, 'solveReverseAST': True,
            'continuationReverseAST': True, 'exactFinalFrequencyAddress': True,
            'adapterSource': ast.unparse(module), 'coefficientAlgorithm': 'numerical dual arithmetic of exact saved expressions',
            'nativeJacobianPackUnpackPolynomialMapsUnchanged': True}


class PilotPair(native.Pair):
    def __init__(self, journal, table, seed, seed_route, coefficients, gauge_cache):
        self.journal, self.table, self.seed = journal, table, seed
        self.seed_route, self.reader = seed_route, coefficients
        self.n = seed['R'].shape[1]
        self.I = np.eye(self.n, dtype=complex)
        folder = 'clusters/' + str(seed['index'])
        self.folder = folder
        self.solve_number = 0
        left = seed['R'].conj().T
        gram = left @ seed['R']
        journal.write(folder + '/gauge-input.pickle', {'gram': gram, 'left': left, 'seed': seed_route,
                      'scope': 'constant real-seed coordinate gauge; not a complex-frequency flux probability'})
        matches = [v for a, b, v in gauge_cache if saved.same(a, gram) and saved.same(b, left)]
        self.gauge = matches[0] if matches else np.linalg.solve(gram, left)
        gauge_cache.append((gram, left, self.gauge))
        journal.write(folder + '/gauge-value.pickle', {'value': self.gauge, 'reusedExactNewCall': bool(matches)})
        # The accepted own physical seed matrix supplies the row-scale input.
        # No old symbolic binding or seed matrix evaluation is replayed.
        P = seed['originalMode']['pencil']
        self.row_scale = 1 + np.max(abs(P), axis=1)
        self.wave_scale = 1 + abs(seed['q'])**2 + 100 * abs(seed['k'])**2
        self.unknown_scale = np.concatenate((np.full(5*self.n, max(1., norm(seed['R']))),
                np.full(self.n**2, max(1., abs(seed['k']))), np.full(self.n**2, max(1., abs(seed['q'])))))
        journal.write(folder + '/prepared.pickle', {'gauge': self.gauge, 'rowScale': self.row_scale,
                      'waveScale': self.wave_scale, 'unknownScale': self.unknown_scale, 'seed': seed_route})

    def coefficients(self, w):
        return self.reader.get(w)

    def wave_equation(self, w, K, Q, dK, dQ, dw):
        return self.reader.wave.matrix_wave(self.reader.wave_root, (w, K, Q), (dw, dK, dQ))

    def check(self, condition, message, local):
        record = {k: v for k, v in local.items() if k != 'self'}
        record.update(condition=condition, message=message)
        self.journal.write(self.current + '/guard-' + str(self.guard_number) + '.pickle', record)
        self.guard_number += 1
        require(condition, message)

    def corrected(self, w, initial, reason):
        self.current = self.folder + '/solves/' + str(self.solve_number)
        self.solve_number += 1
        self.guard_number = 0
        self.journal.write(self.current + '/input.pickle', {'frequency': w, 'initial': initial,
                           'seed': self.seed_route, 'table': self.reader.route, 'reason': reason})
        start = time.monotonic()
        state = self.solve(w, initial)
        value = self.journal.write(self.current + '/value.pickle', state)
        self.journal.json(self.current + '/completed.json', {'wallSeconds': time.monotonic() - start,
                           'value': value, 'guardReceipts': self.guard_number})
        require(state['history'][-1]['scaledJacobianCondition'] < 1e10,
                'pilot conditioning must support useful floating-point precision')
        self.last_state, self.last_value = state, value
        return state

    def reference_point(self, path, state):
        require(state is self.last_state, 'actual immediately completed native point')
        path.relative_to(self.journal.base)
        require(not path.exists() and not path.is_symlink() and not path.parent.is_symlink(), 'fresh point reference')
        target = Path(self.last_value['path']).resolve(strict=True)
        require(saved.digest(target) == self.last_value['sha256'], 'completed point byte identity')
        path.symlink_to(target)
        self.journal.artifacts[str(path.relative_to(self.journal.base))] = {
            'path': str(path), 'sha256': self.last_value['sha256'], 'bytes': target.stat().st_size,
            'canonicalValue': str(target)}


def continue_path(pair, state, target, step, label):
    folder = pair.journal.base / pair.folder / label.replace(' ', '-')
    folder.mkdir()
    return CONTINUE(folder, pair, state, target, step, label)[0]


def evaluator_check(journal, table, bound_tangent, transport, coefficients):
    source = table['source']
    coordinates = (source['frequency'], source['momentum'], source['radical'])
    point = (1. - .0075j, .125 + .02j, .25 + 2.5j)
    journal.write('evaluator-check/input.pickle', {'table': coefficients.route, 'point': point,
                  'boundTangent': bound_tangent, 'transport': transport, 'step': 1e-6})
    program = NumericExpressions(coordinates)
    pencil = [program.add(v) for v in source['livePencil']]
    tangent = [program.add(v) for v in bound_tangent]
    tr = program.add(transport)
    t = program.scalar([tr], point)[0][0]
    value, derivative = program.scalar(pencil, point, (1., 0., t))
    saved_tangent = program.scalar(tangent, point)[0]
    h = 1e-6
    plus = tuple(v + h*d for v, d in zip(point, (1., 0., t)))
    minus = tuple(v - h*d for v, d in zip(point, (1., 0., t)))
    p = np.asarray(program.scalar(pencil, plus)[0], complex)
    q = np.asarray(program.scalar(pencil, minus)[0], complex)
    derivative = np.asarray(derivative, complex)
    finite = (p - q) / (2*h)
    residual = derivative - np.asarray(saved_tangent, complex)
    fd = (finite - derivative) / (1 + norm(derivative))
    values = np.zeros((5, 5), complex)
    tangents = np.zeros_like(values)
    def terms(rows):
        v, d = 0j, 0j
        for (a, b), c in rows:
            v += c * point[1]**a * point[2]**b
            if b:
                d += c * point[1]**a * b * point[2]**(b-1) * t
        return v, d
    for i, j, nt, dt, dnt, ddt in coefficients.get(point[0])[0]:
        N, Nq = terms(nt)
        D, Dq = terms(dt)
        Nw, _ = terms(dnt)
        Dw, _ = terms(ddt)
        values[i, j] = N / D
        tangents[i, j] = (Nw + Nq) / D - N * (Dw + Dq) / D**2
    reconstruction = values.ravel() - np.asarray(value, complex)
    table_tangent = tangents.ravel() - derivative
    mutation = values.copy()
    mutation[0, 0] += 1
    journal.write('evaluator-check/value.pickle', {'pencil': value, 'directionalDerivative': derivative,
                  'acceptedBoundTangentEvaluation': saved_tangent, 'tangentDifference': residual,
                  'finiteDifference': finite, 'scaledFiniteDifferenceResidual': fd,
                  'rationalEntryValues': values, 'rationalEntryTangents': tangents,
                  'rationalValueDifference': reconstruction, 'rationalTangentDifference': table_tangent,
                  'actualChangedCoefficientValue': mutation})
    require(norm(residual) < 1e-8 * (1 + norm(derivative)), 'new numerical evaluator agrees with saved full native tangent')
    require(norm(fd) < 2e-6, 'new numerical evaluator responds to independent centered difference')
    require(norm(reconstruction) < 1e-9 * (1+norm(value)), 'actual 25 rational entry values join own live pencil')
    require(norm(table_tangent) < 1e-8 * (1+norm(derivative)), 'actual rational coefficient duals join own full tangent')
    require(norm(mutation.ravel() - np.asarray(value, complex)) > .5, 'actual coefficient mutation responds')
    K = np.array([[point[1], .01], [0., point[1]]], complex)
    Q = np.array([[point[2], 0.], [.003, point[2]]], complex)
    dK = np.array([[.3, .2], [.1, -.1]], complex)
    dQ = np.array([[.2j, .1], [-.1, .2]], complex)
    dw = .2 + .1j
    journal.write('evaluator-check/wave-input.pickle', {'waveProgram': journal.artifacts['numerical-evaluator/program.pickle'],
                  'frequency': point[0], 'K': K, 'Q': Q, 'dK': dK, 'dQ': dQ, 'dw': dw, 'step': h,
                  'scope': 'Numerical Frechet instrument at declared-coordinate test matrices, not a physical invariant pair.'})
    wave = coefficients.wave
    W, dW = wave.matrix_wave(coefficients.wave_root, (point[0], K, Q), (dw, dK, dQ))
    Z = np.zeros_like(K)
    WP = wave.matrix_wave(coefficients.wave_root, (point[0]+h*dw, K+h*dK, Q+h*dQ), (0., Z, Z))[0]
    WM = wave.matrix_wave(coefficients.wave_root, (point[0]-h*dw, K-h*dK, Q-h*dQ), (0., Z, Z))[0]
    wave_difference = ((WP-WM)/(2*h) - dW)/(1+norm(dW))
    journal.write('evaluator-check/wave-value.pickle', {'value': W, 'frechet': dW,
                  'centeredDifference': (WP-WM)/(2*h), 'scaledDifference': wave_difference})
    require(norm(wave_difference) < 2e-6, 'new acoustic matrix evaluator centered Frechet check')
    return {'savedTangentRelativeDifference': norm(residual) / (1+norm(derivative)), 'finiteDifferenceResidual': norm(fd),
            'rationalValueRelativeDifference': norm(reconstruction)/(1+norm(value)),
            'rationalTangentRelativeDifference': norm(table_tangent)/(1+norm(derivative)),
            'waveFrechetResidual': norm(wave_difference), 'coefficientMutationResponded': True}


def prohibit_old_science():
    def forbidden(*args, **kwargs):
        raise RuntimeError('pilot cannot replay accepted symbolic construction')
    for name in ('diff', 'lambdify', 'cancel', 'factor', 'expand', 'simplify', 'solve', 'resultant', 'gcd', 'gcdex', 'integrate', 'fraction'):
        setattr(sp, name, forbidden)
    for name in ('load', 'seeds', 'focused', 'main'):
        setattr(native, name, forbidden)
    native.Pair.__init__ = forbidden


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-directory', type=Path, required=True)
    base = parser.parse_args().run_directory.resolve()
    base.relative_to(REPO / '_scratch/s11c')
    base.mkdir(parents=True, exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    signal.alarm(900)
    start = time.monotonic()
    reader, journal = saved.Reader(), Journal(base)
    cp = reader.json(saved.CP, saved.CP_SHA)
    origin = Path(cp['runDirectory'])
    ready = reader.json(READY / 'complete/checks.json', READY_SHA)
    ready_inventory = reader.json(READY / 'completed-input-artifact-inventory.json', READY_INVENTORY_SHA)
    require(ready_inventory['checksSha256'] == READY_SHA, 'completed metadata inventory/checks identity')
    reader.retain(READY / 'completion-inspection.json')
    require(ready['clusterUses'] == 40 and ready['completeBaselineCallerCandidates'] == 30, 'actual completed input inspection')
    joins = saved.source_join(reader, cp, origin)
    for path in (Path(__file__).resolve(), PLAN):
        reader.retain(path)
    global CONTINUE
    equation, solve, CONTINUE, adapter = native_adapters()
    PilotPair.equation, PilotPair.solve = equation, solve
    journal.json('native-callers.json', {'sources': joins, 'adapter': adapter})
    prohibit_old_science()

    def packet(name):
        return reader.packet(origin / name, cp['artifacts'][name]['sha256'])

    table_name = 'new-rational-end/right-rational-end.pickle'
    table_route = reader.retain(origin / table_name, cp['artifacts'][table_name]['sha256'])
    table = packet(table_name)
    selection_name = 'end-input-cases/' + OWNER + '/right/whole-cluster-inputs.pickle'
    selection = packet(selection_name)
    full_name = 'end-input-cases/' + OWNER + '/right/full-inputs.pickle'
    full = packet(full_name)
    context_name = 'end-input-cases/' + OWNER + '/context-pairs.pickle'
    reader.retain(origin / context_name, cp['artifacts'][context_name]['sha256'])
    tangent = packet('binding-operations/FREQUENCY_PENCIL_PLUS/value.pickle')
    transport = packet('binding-operations/FREQUENCY_TRANSPORT/value.pickle')
    seed_evidence = {}
    for name in ('seed-matrix-join/input.pickle', 'seed-matrix-join/normalized-residual.pickle',
                 'native-prefix-checks.json', 'saved-tangent-native-join.pickle', 'saved-transport-native-join.pickle'):
        seed_evidence[name] = reader.retain(origin / name, cp['artifacts'][name]['sha256'])
    journal.json('accepted-seed-and-tangent-source-routes.json', seed_evidence)
    require(saved.same(full['address'], (OWNER, 'RIGHT')), 'own physical end source')
    journal.json('pilot-scope.json', {'owner': [OWNER, 'RIGHT'], 'source': table_route,
        'fullInput': reader.retain(origin / full_name, cp['artifacts'][full_name]['sha256']),
        'selection': reader.retain(origin / selection_name, cp['artifacts'][selection_name]['sha256']),
        'frequencyPath': ['1+0i', '1-0.005i', '1-0.01i'], 'coarseStep': .01, 'fineStep': .005,
        'clusters': list(selection['groups']), 'directions': len(selection['selectedItems']),
        'sheet': 'Saved own lifted acoustic Q and whole selected outgoing/incoming clusters, tracked continuously without root reselection.',
        'chart': 'Fixed real-seed R gauge; full commuting K,Q invariant pair, finite regulator and approximate boundary.',
        'comparison': 'one coarse/fine path comparison only; no contour extension or certified disk',
        'observablePrecision': 'about 1% resolved scattering observables later; amplitude1e-4/current1e-6 absolute target',
        'localInstrumentLimits': {'nativeResidual': 2e-12, 'stepScaledDifference': 1e-6, 'mapRelativeDifference': 1e-4},
        'stopping': 'Stop after five clusters and one map comparison; stop early if first-cluster measured cost cannot fit remaining 900-second guard. Unstable results remain unresolved.',
        'priorCostContext': 'Accepted baseline focused 10 clusters / 30 points took about 27.6 seconds; this new evaluator has a separately measured first-cluster pilot.'})
    coefficients = Coefficients(journal, table, table_route)
    evaluator = evaluator_check(journal, table, tangent, transport, coefficients)
    modes = {v['info']['INDEX']: v for v in selection['allCandidateModes']}
    incoming = {v['RECORD_INDEX'] for v in selection['channel']['incoming']}
    clusters, gauge_cache = [], []
    for position, (index, items) in enumerate(selection['groups'].items()):
        descriptor_name = 'cases/' + OWNER + '/right/cluster-' + str(index) + '.json'
        descriptor = reader.json(READY / 'complete' / descriptor_name, ready_inventory['artifacts'][descriptor_name]['sha256'])
        require(descriptor['savedWholeBasisKeys'] and not any(x['completeSavedCallerMatch'] for x in descriptor['baselineCandidates']), 'genuinely distinct own whole-cluster input')
        mode = modes[index]
        R = mode[descriptor['savedWholeBasisKeys'][0]]
        seed = {'index': index, 'R': R, 'k': complex(mode['info']['K']), 'q': complex(mode['info']['Q']),
                'kind': items[0]['kind'], 'direction': 'incoming' if index in incoming else 'outgoing',
                'items': items, 'originalMode': mode}
        seed_route = journal.write('clusters/' + str(index) + '/seed-input.pickle', {
                     'nativeSeed': seed, 'rawK': mode['info']['K'], 'rawQ': mode['info']['Q'],
                     'source': table_route, 'inputDescriptor': descriptor,
                     'fieldUnits': full['signature']['fieldUnits'], 'equationUnits': full['signature']['equationUnits']})
        tick = time.monotonic()
        pair = PilotPair(journal, table, seed, seed_route, coefficients, gauge_cache)
        initial = pair.initial()
        state = pair.corrected(1.+0j, initial, 'new own invariant-pair Jacobian/tangent at accepted unchanged mode seed')
        delta = state['x'] - initial
        journal.write('clusters/' + str(index) + '/seed-unchanged.pickle', delta)
        require(norm(delta) == 0, 'accepted own mode reused without new root correction')
        coarse = continue_path(pair, state, TARGET, .01, 'coarse path')
        fine = continue_path(pair, state, TARGET, .005, 'fine path')
        difference = (coarse['x'] - fine['x']) / pair.unknown_scale
        journal.write('clusters/' + str(index) + '/comparison.pickle', {'coarse': coarse, 'fine': fine,
                      'scaledDifference': difference, 'elapsedSeconds': time.monotonic() - tick})
        require(norm(difference) < 1e-6, 'bounded whole-cluster step comparison')
        clusters.append({'seed': seed, 'seedState': state, 'coarse': coarse, 'fine': fine,
                         'difference': norm(difference), 'correctedCalls': pair.solve_number})
        if position == 0:
            cost = time.monotonic() - tick
            remaining = 900 - (time.monotonic() - start)
            decision = {'firstClusterSeconds': cost, 'remainingClusters': len(selection['groups']) - 1,
                        'remainingBudgetSeconds': remaining, 'estimatedRemainingWithFactorThree': 3 * cost * (len(selection['groups']) - 1),
                        'stopping': 'No automatic extension or retry; preserve completed points.'}
            journal.json('measured-pilot-cost.json', decision)
            require(decision['estimatedRemainingWithFactorThree'] + 60 < remaining, 'measured pilot cost fits bounded remaining batch')
    maps = {}
    map_calls = []
    for label in ('coarse', 'fine'):
        arguments = [{'seed': v['seed'], 'state': v[label]} for v in clusters]
        journal.write('maps/' + label + '-input.pickle', {'clusters': arguments, 'xposition': 64., 'source': table_route})
        consumed = [{'direction': v['seed']['direction'], 'kind': v['seed']['kind'],
                     'R': v['state']['R'], 'K': v['state']['K']} for v in arguments]
        matches = [value for prior, value in map_calls if saved.same(consumed, prior)]
        result = matches[0] if matches else native.maps(arguments, 64.)
        journal.write('maps/' + label + '-value.pickle', result)
        journal.json('maps/' + label + '-completed.json', {'reusedExactCompletedMapInput': bool(matches)})
        map_calls.append((consumed, result))
        maps[label] = result
        require(norm(result['inverseResidual']) < 1e-10 and norm(result['traceResidual']) < 1e-10, 'own complete end-map residuals')
    differences = {name: maps['fine'][name] - maps['coarse'][name] for name in
                   ('trace', 'forcingAtCommonOrigin', 'observationAtCommonOrigin', 'directIncomingSubtraction')}
    relative = {name: norm(value) / (1 + norm(maps['fine'][name])) for name, value in differences.items()}
    journal.write('maps/comparison.pickle', {'differences': differences, 'relative': relative})
    require(max(relative.values()) < 1e-4, 'one selected complete-map path refinement')
    reader.postcheck()
    journal.json('inputs.json', {'acceptedCheckpoint': reader.retain(saved.CP, saved.CP_SHA),
                 'savedInputInspectionChecks': reader.retain(READY / 'complete/checks.json', READY_SHA),
                 'consumedRoutes': reader.routes, 'allConsumedInputsUnchanged': True})
    checks = {'status': 'COMPLETED_BOUNDED_RHOBR_RIGHT_END_PILOT', 'owner': [OWNER, 'RIGHT'],
              'clusters': len(clusters), 'directions': sum(v['seed']['R'].shape[1] for v in clusters),
              'newCorrectedStates': sum(v['correctedCalls'] for v in clusters),
              'newNonseedPathPoints': sum(v['correctedCalls'] - 1 for v in clusters),
              'frequency': repr(TARGET), 'newSymbolicDerivativesOrLambdify': 0,
              'coefficientFrequencyEvaluations': len(coefficients.cache), 'evaluatorCheck': evaluator,
              'maximumStepDifference': max(v['difference'] for v in clusters), 'mapRelativeDifferences': relative,
              'maximumEndResidual': max(norm(v['fine']['residual']) for v in clusters),
              'maximumCondition': max(v['fine']['history'][-1]['scaledJacobianCondition'] for v in clusters),
              'minimumDenominatorSingularValue': min(v['fine']['history'][-1]['minimumDenominatorSingularValue'] for v in clusters),
              'fineMapCondition': maps['fine']['condition'], 'sourceAndInputsUnchanged': True,
              'artifacts': journal.artifacts, 'wallSeconds': time.monotonic() - start,
              'peakRssKiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              'scope': 'Bounded analog toy-model end continuation and maps at one complex frequency. No scattering result, pole, certified complex domain or numerical reuse partition is accepted here.'}
    saved.save(base, 'checks.json', checks)
    signal.alarm(0)
    print(json.dumps(checks, indent=2))


if __name__ == '__main__':
    main()
