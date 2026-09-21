#!/usr/bin/env python3
"""Supply saved native Gauss rules while retaining the physical-solver blocks."""
import argparse
import ast
import copy
import inspect
import json
from pathlib import Path
import resource
import shutil
import signal
import sys
import time

import S11c_d_remaining_case_coordinate_matrices_production as prior

h = prior.h
f, np, modes = h.f, h.np, h.modes
ORIGIN = f.STORE/'s11c-remaining-case-coordinate-20260921/matrices/production/complete'
RULE_FOCUS = f.STORE/'s11c-remaining-case-coordinate-20260921/matrices/rule-repair/complete'
PLAN = f.M/'S11c_d_remaining_case_coordinate_matrices_recovery_plan.md'
REPAIR = f.M/'S11c_d_remaining_case_coordinate_matrices_rule_repair.json'
ORIGINAL_LEGGAUSS = np.polynomial.legendre.leggauss
ACTIVE_CACHE = None


def orders_from(settings):
    values = [settings[k] for k in ('outerOrder', 'panelOrder', 'sourceNodes', 'profileNodes')]
    values += list(settings['innerOrders'])
    f.require(all(type(v) is int and v >= 2 for v in values), 'actual positive integer quadrature orders')
    return sorted(set(values))


def rule_proof(order, packet):
    x, w = packet['nodes'], packet['weights']
    f.require(packet['order'] == order and x.shape == w.shape == (order,), 'actual native rule order/shape')
    f.require(x.dtype == w.dtype == np.dtype('float64') and np.isfinite(x).all() and np.isfinite(w).all(), 'finite native rule arrays')
    f.require(np.all(np.diff(x) > 0) and np.all(abs(x) < 1) and np.all(w > 0), 'ordered interior nodes and positive weights')
    f.require(np.array_equal(x, -x[::-1]) and np.array_equal(w, w[::-1]), 'literal native rule symmetry')
    # Independent moments on the standard interval, through its actual degree.
    powers = np.arange(2*order)
    actual = (x[:, None]**powers[None, :]).T @ w
    exact = np.where(powers % 2 == 0, 2/(powers+1), 0.)
    residual = actual-exact
    f.require(float(np.max(abs(residual))) < 1e-11, 'all unit-rule polynomial moments')
    return {'order': order, 'moments': actual, 'exactMoments': exact, 'residuals': residual,
            'maximumMomentResidual': float(np.max(abs(residual)))}


class SavedGauss:
    def __init__(self, root, orders):
        self.orders = tuple(orders); self.rules = {}; self.calls = {v: 0 for v in orders}
        for order in orders:
            packet = f.unpickle(root/f'order-{order}.pickle')
            f.require(packet['order'] == order, 'saved unit rule address')
            self.rules[order] = packet

    def __call__(self, order):
        f.require(type(order) is int and order in self.orders, 'only exact approved native quadrature orders')
        packet = self.rules[order]; self.calls[order] += 1
        # Native callers may cache or mutate arrays; their independent copies
        # preserve the immutable saved source arrays.
        return packet['nodes'].copy(), packet['weights'].copy()


def numpy_sources():
    modules = (np, np.polynomial.legendre, np.polynomial.polyutils,
               np.linalg, np.linalg.linalg, np.linalg._umath_linalg)
    paths = {str(Path(module.__file__).resolve()) for module in modules}
    paths.add(str(Path(inspect.getsourcefile(ORIGINAL_LEGGAUSS)).resolve()))
    return {path: f.digest(Path(path)) for path in sorted(paths)}


def native_joins():
    _, inventory = prior.production_construct()
    tree = ast.parse(inspect.getsource(h.main))
    blocks = [n for n in ast.walk(tree) if isinstance(n, ast.Call) and ast.unparse(n.func) == 'setattr'
              and n.args and ast.unparse(n.args[0]) == 'np.linalg']
    f.require(len(blocks) == 1, 'original complete physical linear-algebra block retained')
    return {'inventoryJoin': inventory, 'originalMainAstSha256': h.body(h.main),
            'originalLeggaussAstSha256': h.body(ORIGINAL_LEGGAUSS),
            'nativeSourceRuleAstSha256': h.body(h.engine.BoundedSourceFourierQuadrature.rule),
            'nativeMomentumRuleAstSha256': h.body(h.engine.BoundedSourceFourierQuadrature.FiniteMomentum.rule),
            'originalNativeJoins': h.native_joins(), 'linearAlgebraBlockUnchanged': True,
            'scope': 'The unchanged native rule callers receive exact saved leggauss results; no physical eigensolver is enabled.'}


def rule_focus(base):
    base.mkdir(parents=True, exist_ok=False); start = time.monotonic()
    old = json.loads((ORIGIN/'inputs.json').read_text()); cp = json.loads(h.FCP.read_text())
    f.require(cp['status'] == 'ACCEPTED_CASE_MATERIAL_MATRIX_INPUTS' and
              f.digest(Path(cp['runDirectory'])/'checks.json') == cp['checksSha256'], 'accepted matrix input focus')
    f.require(old['settings'] == cp['settings'], 'all approved settings unchanged')
    f.require(not list((ORIGIN/'cases').rglob('*.pickle')), 'failure preceded any completed new scientific packet')
    inv = json.loads((ORIGIN.parent/'coordinate_matrices_construct.invocation.json').read_text())
    f.require(inv['exitCode'] == 1 and 'leggauss' in (ORIGIN.parent/'coordinate_matrices_construct.stderr').read_text(), 'actual quadrature dependency failure')
    expected = {}
    for n, sha in old['sourceFiles'].items():
        expected[str(f.ROOT/n)] = sha; expected[str(ORIGIN/'source'/n)] = sha
    for n, sha in old['copiedInputs'].items(): expected[str(ORIGIN/n)] = sha
    for n, item in cp['artifacts'].items():
        f.require(f.digest(ORIGIN/n) == item['sha256'], 'all accepted focus copies preserved')
    for p in (ORIGIN/'inputs.json', h.FCP, Path(__file__).resolve(), PLAN, Path(prior.__file__).resolve()):
        expected[str(p)] = f.digest(p)
    for name in ('active.json', 'coordinate_matrices_construct.invocation.json', 'coordinate_matrices_construct.stdout',
                 'coordinate_matrices_construct.stderr', 'guard.stdout', 'guard.stderr', 'resource-guard/outcome.json',
                 'resource-guard/child-outcome.json', 'resource-guard/effective-limits.json',
                 'resource-guard/limit-validation.json', 'resource-guard/resource-samples.jsonl'):
        path = ORIGIN.parent/name; expected[str(path)] = f.digest(path)
    expected.update(numpy_sources())
    for n, sha in expected.items(): f.require(f.digest(Path(n)) == sha, ('repair original/current pre hash', n))
    source_files = {str(p.relative_to(f.ROOT)): f.digest(p) for p in (Path(__file__).resolve(), PLAN, Path(prior.__file__).resolve())}
    f.save(base/'inputs.json', {'inputPackets': expected, 'sourceFiles': source_files, 'settings': old['settings'],
                             'originalDirectory': str(ORIGIN), 'acceptedFocusChecksSha256': cp['checksSha256']})
    for path in (Path(__file__).resolve(), PLAN, Path(prior.__file__).resolve()):
        target = base/'source'/path.relative_to(f.ROOT); target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        f.require(f.digest(target) == expected[str(path)], 'frozen rule repair source')
    joins = native_joins(); f.save(base/'native-rule-joins.json', joins)
    orders = orders_from(old['settings']); rules = base/'rules'; rules.mkdir(); proofs = {}; controls = []
    for order in orders:
        # Only standard-interval node/weight generation; no source, profile,
        # momentum operator, trial basis or physical pencil is evaluated.
        x, w = ORIGINAL_LEGGAUSS(order)
        packet = {'order': order, 'nodes': x, 'weights': w, 'numpyVersion': np.__version__,
                  'producerAstSha256': joins['originalLeggaussAstSha256']}
        f.atomic_pickle(rules/f'order-{order}.pickle', packet)
        proof = rule_proof(order, packet); f.atomic_pickle(rules/f'proof-{order}.pickle', proof); proofs[order] = proof['maximumMomentResidual']
        for kind in ('weight', 'node', 'order'):
            changed = copy.deepcopy(packet)
            if kind == 'weight': changed['weights'][0] *= 1.001
            elif kind == 'node': changed['nodes'][0] += 0.001
            else: changed['order'] += 1
            rejected = False
            try: rule_proof(order, changed)
            except (ValueError, AssertionError, RuntimeError): rejected = True
            f.require(rejected, ('actual unit-rule mutation rejected', order, kind))
            controls.append({'order': order, 'kind': kind, 'changed': changed, 'rejected': rejected})
    cache = SavedGauss(rules, orders)
    np.polynomial.legendre.leggauss = cache
    def forbidden(*a, **k): raise RuntimeError('physical linear algebra and scientific constructors disabled in Gauss cache regression')
    for name in ('solve', 'inv', 'pinv', 'svd', 'eig', 'eigh', 'eigvals', 'eigvalsh'): setattr(np.linalg, name, forbidden)
    solver_blocks = {}
    for name in ('solve', 'inv', 'pinv', 'svd', 'eig', 'eigh', 'eigvals', 'eigvalsh'):
        rejected = False
        try: getattr(np.linalg, name)()
        except RuntimeError: rejected = True
        f.require(rejected, 'physical solver entry remains disabled'); solver_blocks[name] = rejected
    h.engine.NumericalReducedAction.bind = forbidden
    h.c.Chart.__init__ = h.c.Chart.bind_image = h.c.Chart.image = h.c.Chart.basis = forbidden
    h.c.source.coordinate_change = forbidden; h.interior.assemble = h.matrices.direct_cells = forbidden
    h.material_native.MaterialProduction.prepare_basis = h.material_native.MaterialProduction.matrix_group = forbidden
    for order in orders:
        x, w = cache(order); packet = cache.rules[order]
        f.require(np.array_equal(x, packet['nodes']) and np.array_equal(w, packet['weights']), 'literal saved native Gauss return arrays')
        x[0] += 1; f.require(np.array_equal(cache(order)[0], packet['nodes']), 'caller mutations cannot change saved rules')
        bound = old['settings']['sourceBound'] if order == old['settings']['sourceNodes'] else old['settings']['profileBound']
        panels = [-bound, -10., 0., 10., bound]
        nodes, weights = h.engine.BoundedSourceFourierQuadrature.rule(panels, order)
        points = sorted(set(map(float, panels)))
        expected_nodes = np.concatenate([(a+b)/2+(b-a)*packet['nodes']/2 for a, b in zip(points, points[1:])])
        expected_weights = np.concatenate([(b-a)*packet['weights']/2 for a, b in zip(points, points[1:])])
        f.atomic_pickle(base/f'panel-{order}.pickle', {'points': points, 'order': order, 'nodes': nodes, 'weights': weights,
                        'expectedNodes': expected_nodes, 'expectedWeights': expected_weights})
        f.require(np.array_equal(nodes, expected_nodes) and np.array_equal(weights, expected_weights), 'unchanged native affine panel rule with every physical solver blocked')
    momentum = h.engine.BoundedSourceFourierQuadrature.FiniteMomentum.__new__(h.engine.BoundedSourceFourierQuadrature.FiniteMomentum)
    momentum._legendre = {}
    bound = old['settings']['momentumBound']
    for order in sorted(set([old['settings']['outerOrder'], old['settings']['panelOrder']]+old['settings']['innerOrders'])):
        nodes, weights, points = momentum.rule(-bound, bound, order)
        packet = cache.rules[order]
        expected_nodes, expected_weights = bound*packet['nodes'], bound*packet['weights']
        f.atomic_pickle(base/f'momentum-panel-{order}.pickle', {'points': points, 'order': order,
                        'nodes': nodes, 'weights': weights, 'expectedNodes': expected_nodes, 'expectedWeights': expected_weights})
        f.require(np.array_equal(nodes, expected_nodes) and np.array_equal(weights, expected_weights),
                  'unchanged native momentum panel rule with every physical solver blocked')
    for value in (True, float(orders[0]), orders[0]+1):
        rejected = False
        try: cache(value)
        except (ValueError, AssertionError, RuntimeError): rejected = True
        f.require(rejected, 'wrong rule type/order rejected'); controls.append({'changedOrder': value, 'rejected': rejected})
    f.atomic_pickle(base/'mutation-controls.pickle', controls)
    for n, sha in expected.items(): f.require(f.digest(Path(n)) == sha, ('repair pre/post hash', n))
    for n, sha in source_files.items(): f.require(f.digest(base/'source'/n) == sha, 'unchanged frozen repair source')
    result = {'status': 'VALIDATED_SAVED_NATIVE_GAUSS_RULES', 'runDirectory': str(base), 'orders': orders,
              'settings': old['settings'], 'nativeJoins': joins, 'momentResiduals': proofs, 'controls': len(controls),
              'acceptedFocusCopies': len(cp['artifacts']), 'inputPackets': expected, 'sourceFiles': source_files,
              'numpyVersion': np.__version__, 'physicalSolverRejections': solver_blocks,
              'newPhysicalIntegrationsBindingsMatricesModesCurrentsSolves': 0, 'allPhysicalLinearAlgebraBlocked': True,
              'artifacts': {str(p.relative_to(base)): {'sha256': f.digest(p), 'bytes': p.stat().st_size}
                            for p in base.rglob('*') if p.is_file() and 'source' not in p.relative_to(base).parts
                            and p not in (base/'inputs.json', base/'checks.json')},
              'wallSeconds': time.monotonic()-start}
    f.save(base/'checks.json', result); return result


def recover_load(base, resume):
    global ACTIVE_CACHE
    cp = json.loads(h.FCP.read_text()); repair = json.loads(REPAIR.read_text())
    f.require(cp['status'] == 'ACCEPTED_CASE_MATERIAL_MATRIX_INPUTS' and str(resume) == cp['runDirectory'] and
              f.digest(resume/'checks.json') == cp['checksSha256'], 'original accepted focus retained')
    f.require(repair['status'] == 'ACCEPTED_NATIVE_QUADRATURE_RULE_REUSE' and
              repair['runDirectory'] == str(RULE_FOCUS) and repair['checksSha256'] == f.digest(RULE_FOCUS/'checks.json'), 'accepted exact native Gauss rule cache')
    rule_checks = json.loads((RULE_FOCUS/'checks.json').read_text())
    old = json.loads((ORIGIN/'inputs.json').read_text())
    f.require(old['settings'] == rule_checks['settings'] and not list((ORIGIN/'cases').rglob('*.pickle')), 'saved original inputs, no completed numerical construction to repeat')
    manifest = copy.deepcopy(old); manifest.update(runDirectory=str(base), copiedInputs={})
    for n, sha in old['sourceFiles'].items(): f.require(f.digest(f.ROOT/n) == f.digest(ORIGIN/'source'/n) == sha, 'original/current/frozen helper join')
    for n, sha in old['copiedInputs'].items(): f.require(f.digest(ORIGIN/n) == sha, 'completed original load copy')
    for path in sorted(ORIGIN.rglob('*')):
        if not path.is_file() or 'source' in path.relative_to(ORIGIN).parts or path == ORIGIN/'inputs.json': continue
        modes.retain(path, base/path.relative_to(ORIGIN), manifest)
    modes.retain(ORIGIN/'inputs.json', base/'original-production-inputs.json', manifest)
    for n in ('active.json', 'coordinate_matrices_construct.invocation.json', 'coordinate_matrices_construct.stdout',
              'coordinate_matrices_construct.stderr', 'guard.stdout', 'guard.stderr', 'launch.json',
              'resource-guard/outcome.json', 'resource-guard/child-outcome.json', 'resource-guard/effective-limits.json',
              'resource-guard/limit-validation.json', 'resource-guard/resource-samples.jsonl'):
        modes.retain(ORIGIN.parent/n, base/'original-production-logs'/n, manifest)
    for n, item in cp['artifacts'].items(): f.require(f.digest(base/n) == item['sha256'], 'every accepted focused artifact reused unchanged')
    for n, item in rule_checks['artifacts'].items(): modes.retain(RULE_FOCUS/n, base/'quadrature-rule-reuse'/n, manifest, item['sha256'])
    modes.retain(RULE_FOCUS/'checks.json', base/'quadrature-rule-reuse/checks.json', manifest)
    modes.retain(RULE_FOCUS/'inputs.json', base/'quadrature-rule-reuse/inputs.json', manifest)
    for n, sha in rule_checks['inputPackets'].items():
        f.require(f.digest(Path(n)) == sha, 'exact accepted Gauss source/helper/environment'); manifest['inputPackets'][n] = sha
    for path in (Path(__file__).resolve(), PLAN, REPAIR): manifest['sourceFiles'][str(path.relative_to(f.ROOT))] = f.digest(path)
    for n, sha in manifest['sourceFiles'].items():
        target = base/'source'/n; target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(f.ROOT/n, target); f.require(f.digest(target) == sha, 'recovery frozen source')
    ACTIVE_CACHE = SavedGauss(base/'quadrature-rule-reuse/rules', rule_checks['orders'])
    np.polynomial.legendre.leggauss = ACTIVE_CACHE
    manifest['completedInputReuse'] = {'directory': str(ORIGIN), 'originalInputsSha256': f.digest(ORIGIN/'inputs.json'),
        'acceptedFocusArtifacts': len(cp['artifacts']), 'originalSourceFiles': len(old['sourceFiles']),
        'originalCopiedInputs': len(old['copiedInputs']), 'newNumericalPacketsBeforeFailure': 0,
        'nativeRuleChecksSha256': f.digest(RULE_FOCUS/'checks.json'), 'nativeJoins': rule_checks['nativeJoins'],
        'scope': 'Saved completed input load and native standard-interval rules; no physical eigensolver is re-enabled.'}
    f.save(base/'completed-input-reuse.json', manifest['completedInputReuse']); f.save(base/'inputs.json', manifest)
    h.hash_check(base, manifest)
    return manifest, tuple(cp['cases'])


def main():
    if '--rule-focus' in sys.argv:
        parser = argparse.ArgumentParser(); parser.add_argument('--rule-focus', action='store_true'); parser.add_argument('--run-directory', type=Path, required=True)
        args = parser.parse_args(); base = args.run_directory.resolve(); base.relative_to(f.STORE)
        resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3)); signal.alarm(900)
        result = rule_focus(base); signal.alarm(0); print(json.dumps(result, indent=2)); return
    coordinator, _ = prior.production_construct()
    def construct(base, manifest, labels, joins):
        try: return coordinator(base, manifest, labels, joins)
        finally:
            if ACTIVE_CACHE is not None: f.save(base/'quadrature-rule-use.json', {'calls': ACTIVE_CACHE.calls, 'orders': ACTIVE_CACHE.orders,
                                'scope': 'Exact saved native Gauss rules returned; no eigensolver executed by production.'})
    h.load = recover_load; h.construct = construct
    h.main()


if __name__ == '__main__': main()
