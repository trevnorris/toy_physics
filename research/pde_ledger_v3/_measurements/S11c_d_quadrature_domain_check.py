#!/usr/bin/env python3
"""Source-derived denominator, Abel-resolution and profile-tail operands."""
import argparse
import ast
import contextlib
import json
from pathlib import Path
import pickle
import resource
import shutil
import time

import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.special import roots_legendre

import S11c_d_numerical_action_check as native
from S11c_d_numerical_action_check import ROOT, STORE, engine, digest, save, atomic_pickle, number
from S11c_d_output_codec import decoded_lines, restore_emission_index
from ledger_fold import _restore

CHECKPOINT = ROOT/'_measurements/S11c_d_numerical_action_checkpoint.json'
PLAN = ROOT/'_measurements/S11c_d_quadrature_limits_plan.md'
PREFIX = 'QUADRATURE_DOMAIN_LAB_HELD_RHO4_CONSTANT'
SOURCES = (*native.SOURCES, Path(__file__).resolve(), Path(native.__file__).resolve(), CHECKPOINT, PLAN)


def load():
    checkpoint = json.loads(CHECKPOINT.read_text())
    base = Path(checkpoint['runDirectory'])
    for name, item in checkpoint['artifacts'].items():
        if digest(base/name) != item['sha256']:
            raise ValueError(('accepted numerical operand changed', name))
    for name, sha in checkpoint['sourceFiles'].items():
        if digest(base/'source'/name) != sha or digest(ROOT/name) != sha:
            raise ValueError(('accepted numerical constructor source changed', name))
    if digest(ROOT/checkpoint['publication']['path']) != checkpoint['publication']['sha256']:
        raise ValueError('numerical action annex payload changed')
    _, actions, pencil, assembly, _ = native.load()
    adapter = engine.NumericalReducedAction(pencil, assembly, json.loads(native.INPUT.read_text()))
    with (base/'actions.pickle').open('rb') as stream:
        previous = pickle.load(stream)
    return checkpoint, actions, pencil, assembly, adapter, previous


def denominator_records(assembly, adapter):
    r = adapter.r
    momenta = {r.normal_map[g[2]] for g in r.momentum_groups}
    coordinates = momenta | {r.regulator, r.z, r.zp, r.xi}
    seen, bases = set(), set()
    def visit(node):
        if node in seen:
            return
        seen.add(node)
        if node.is_Pow and node.exp.is_negative:
            bound = adapter.bind(node.base)
            if engine.dag_free_symbols(bound) & coordinates:
                bases.add(node.base)
        for child in node.args:
            visit(child)
    for expression in (*assembly['LOCAL_MATRICES'].values(), *assembly['NONLOCAL_INTEGRALS']):
        visit(expression)
    records, pairs = [], set()
    for base in sorted(bases, key=sp.default_sort_key):
        bound = adapter.bind(base)
        real, imaginary = sp.simplify(sp.re(bound)), sp.simplify(sp.im(bound))
        norm = sp.factor(real**2+imaginary**2)
        pair = tuple(sorted(engine.dag_free_symbols(bound) & momenta, key=sp.default_sort_key))
        if bound.has(r.regulator) and pair:
            pairs.add(pair)
        records.append({'source': base, 'bound': bound, 'real': real, 'imaginary': imaginary,
            'normSquared': norm, 'unit': engine.PHYSICAL_METADATA.dimensions.measure(base),
            'realPositive': int(real.is_positive is True), 'realNegative': int(real.is_negative is True),
            'normPositive': int(norm.is_positive is True), 'momenta': pair})
    return records, tuple(sorted(pairs, key=lambda p: tuple(map(str, p))))


def split_rule(lower, upper, center, width, order):
    """Finite quadrature panels centered on a source-derived concentration width."""
    points = {lower, upper}
    if lower < center < upper:
        points.add(center)
    distance = width
    while distance < 2*(upper-lower+abs(center)):
        for p in (center-distance, center+distance):
            if lower < p < upper:
                points.add(p)
        distance *= 2
    nodes, weights = roots_legendre(order)
    positions = sorted(points)
    return (np.concatenate([(b-a)*nodes/2+(a+b)/2 for a, b in zip(positions, positions[1:])]),
            np.concatenate([(b-a)*weights/2 for a, b in zip(positions, positions[1:])]), positions)


def abel_records(adapter, previous):
    r = adapter.r
    k = r.normal_map[r.momentum_groups[0][2]]
    center = sp.Symbol('s11cdQuadratureCenterMomentum', real=True)
    delta = sp.Symbol('s11cdQuadraturePositiveWidth', positive=True)
    dimensions = engine.PHYSICAL_METADATA.dimensions
    unit = dimensions.measure(k)
    dimensions.known[center] = unit; dimensions.known[delta] = unit
    source_integral = sp.Integral(r.abel_even, (r.transfer, -sp.oo, sp.oo))
    transformed = source_integral.transform(r.transfer, (r.ell*(k-center), k))
    density = adapter.bind(transformed.function)
    primitive = sp.integrate(density, k)
    primitive_residual = sp.simplify(sp.diff(primitive, k)-density)
    solutions = sp.solve(sp.Eq(density.subs(k, center+delta), density.subs(k, center)/2), delta)
    if len(solutions) != 1 or solutions[0].is_positive is not True:
        raise ValueError('unresolved source-derived Abel width')
    width = solutions[0]
    width_residual = sp.simplify(density.subs(k, center+width)-density.subs(k, center)/2)
    density_fn = sp.lambdify((k, center, r.regulator), density, 'numpy')
    primitive_fn = sp.lambdify((k, center, r.regulator), primitive, 'numpy')
    width_fn = sp.lambdify(r.regulator, width, 'numpy')
    setting = previous['result']['results'][0]['settings']
    bound = setting['momentumBound']
    regulators = [setting['regulator']/2**j for j in range(3)]
    ordinary, split = [], []
    for count in (12, 16, 32, 64, 128, 256, 512, 1024):
        nodes, weights = roots_legendre(count)
        nodes, weights = bound*nodes, bound*weights
        centers = (0.0, float(nodes[np.argmin(np.abs(nodes))]), 0.7*bound, 0.95*bound)
        for regulator in regulators:
            for location in centers:
                exact = float(primitive_fn(bound, location, regulator)-primitive_fn(-bound, location, regulator))
                value = complex(np.dot(weights, density_fn(nodes, location, regulator)))
                w = float(width_fn(regulator))
                ordinary.append({'count': count, 'regulator': regulator, 'center': location,
                    'width': w, 'exact': exact, 'value': value, 'residual': value-exact})
                # Use both a low and a refined panel order on exactly the
                # same source-derived panel boundaries.
                for order in (8, 16):
                    x, weights_split, panels = split_rule(-bound, bound, location, w, order)
                    value_split = complex(np.dot(weights_split, density_fn(x, location, regulator)))
                    split.append({'parentCount': count, 'order': order, 'count': len(x),
                        'regulator': regulator, 'center': location, 'width': w, 'panels': panels,
                        'exact': exact, 'value': value_split, 'residual': value_split-exact})
    return {'sourceIntegral': source_integral, 'transformed': transformed, 'density': density,
        'primitive': primitive, 'primitiveResidual': primitive_residual, 'width': width,
        'widthResidual': width_residual, 'center': center, 'momentum': k, 'momentumUnit': unit,
        'ordinary': ordinary, 'split': split, 'bound': bound}


def adaptive_complex(function, lower, upper):
    real, er = quad(lambda x: float(np.real(function(x))), lower, upper, epsabs=1e-12, epsrel=1e-12, limit=300)
    imag, ei = quad(lambda x: float(np.imag(function(x))), lower, upper, epsabs=1e-12, epsrel=1e-12, limit=300)
    return complex(real, imag), er+ei


def profile_records(assembly, adapter):
    r = adapter.r
    original = set()
    for integral in assembly['NONLOCAL_INTEGRALS']:
        original.update(g for g in integral.atoms(sp.Integral) if g.limits == ((r.xi, -sp.oo, sp.oo),))
    momenta = {r.normal_map[g[2]] for g in r.momentum_groups}
    records = []
    for integral in sorted(original, key=sp.default_sort_key):
        bound = adapter.bind(integral)
        phases = list(bound.function.atoms(sp.exp))
        if len(phases) != 1:
            raise ValueError('unresolved native profile phase')
        transfer = sp.simplify(sp.I*sp.diff(phases[0].args[0], r.xi))
        legs = sorted(engine.dag_free_symbols(transfer) & momenta, key=sp.default_sort_key)
        solved = sp.solve(transfer-r.transfer, legs[0])
        if len(solved) != 1:
            raise ValueError('profile transfer coordinate is not isolated')
        normalized = bound.function.subs(legs[0], solved[0])
        # Simplify the linear phase while retaining the bounded tanh profile
        # representation for adaptive evaluation arbitrarily far in the tails.
        normalized = normalized.replace(lambda n: n.func == sp.exp,
                                        lambda n: sp.exp(sp.expand(n.args[0])))
        if engine.dag_free_symbols(normalized)-{r.xi, r.transfer}:
            raise ValueError('profile transform retains an unbound leg after transfer isolation')
        function = sp.lambdify((r.xi, r.transfer), normalized, 'numpy', cse=True)
        amplitude = sp.Abs(normalized.subs(r.transfer, 0))
        amplitude_fn = sp.lambdify(r.xi, amplitude, 'numpy')
        samples, tails = [], []
        for cutoff in (6.0, 10.0, 14.0):
            left, el = quad(amplitude_fn, -np.inf, -cutoff, epsabs=1e-12, epsrel=1e-12, limit=300)
            right, er = quad(amplitude_fn, cutoff, np.inf, epsabs=1e-12, epsrel=1e-12, limit=300)
            tails.append({'cutoff': cutoff, 'left': float(left), 'right': float(right), 'errorEstimate': float(el+er)})
        for cutoff in (6.0, 10.0):
            for transfer_value in (-40.0, -20.0, -5.0, -1.0, -0.2, 0.0, 0.2, 1.0, 5.0, 20.0, 40.0):
                f = lambda x: function(x, transfer_value)
                left, el = adaptive_complex(f, -cutoff, 0)
                right, er = adaptive_complex(f, 0, cutoff)
                value = left+right
                for order in (16, 32, 64, 128):
                    x, w = roots_legendre(order)
                    left_nodes, right_nodes = cutoff*(x-1)/2, cutoff*(x+1)/2
                    gauss = complex(np.dot(cutoff*w/2, f(left_nodes)+f(right_nodes)))
                    samples.append({'cutoff': cutoff, 'transfer': transfer_value, 'order': order,
                        'adaptive': value, 'adaptiveErrorEstimate': el+er,
                        'gauss': gauss, 'residual': gauss-value})
        records.append({'source': integral, 'bound': bound, 'transfer': transfer,
            'legIsolation': solved[0], 'normalizedIntegrand': normalized, 'amplitude': amplitude,
            'unit': engine.PHYSICAL_METADATA.dimensions.measure(integral), 'samples': samples, 'tails': tails})
    return records


def emit_result(result, pencil, provenance):
    dimensions = engine.PHYSICAL_METADATA.dimensions
    zero = dimensions.zero
    metadata = engine.FullPencilModes.__new__(engine.FullPencilModes); metadata.r = pencil.r
    def numeric(suffix, value, unit):
        body = number(value)
        if engine.dag_size(body) < 60:
            emitted = body
        elif all(expression.is_number for _, expression in engine.leaves(body)):
            emitted = metadata.compact_fingerprint(body)
        else:
            # These are bound symbolic operands, not evaluated integrals.
            emitted = engine.carrier_fingerprint(body)
        engine.emit(PREFIX+'_'+suffix, emitted)
        engine.emit('METADATA_'+PREFIX+'_'+suffix, metadata.numeric_metadata(body, unit))
    engine.physical(PREFIX+'_PROVENANCE', provenance)
    for i, record in enumerate(result['denominators']):
        unit = record['unit']
        engine.fingerprinted(PREFIX+'_SOURCE_DENOMINATOR_'+str(i), record['source'])
        numeric('BOUND_DENOMINATOR_PARTS_'+str(i), tuple(record[k] for k in ('bound', 'real', 'imaginary')),
                lambda p: unit)
        numeric('DENOMINATOR_NORM_SQUARED_'+str(i), record['normSquared'], lambda p: tuple(2*n for n in unit))
        numeric('EXACT_SIGN_FLAGS_'+str(i), tuple(record[k] for k in ('realPositive', 'realNegative', 'normPositive')),
                lambda p: zero)
    engine.physical(PREFIX+'_ABEL_TRANSFER_PAIRS', result['pairs'])
    abel = result['abel']; k_unit = abel['momentumUnit']
    engine.fingerprinted(PREFIX+'_SOURCE_ABEL_TRANSFORM', abel['sourceIntegral'])
    engine.fingerprinted(PREFIX+'_CHANGED_VARIABLE_ABEL_TRANSFORM', abel['transformed'])
    numeric('ABEL_DENSITY', abel['density'], lambda p: tuple(-u for u in k_unit))
    numeric('ABEL_PRIMITIVE', abel['primitive'], lambda p: zero)
    numeric('ABEL_PRIMITIVE_RESIDUAL', abel['primitiveResidual'], lambda p: tuple(-u for u in k_unit))
    numeric('ABEL_WIDTH', abel['width'], lambda p: k_unit)
    numeric('ABEL_HALF_HEIGHT_RESIDUAL', abel['widthResidual'], lambda p: tuple(-u for u in k_unit))
    for name in ('ordinary', 'split'):
        rows = abel[name]
        numeric('ABEL_'+name.upper()+'_SETTINGS', [(r['count'], r['regulator'], r['center'], r['width']) for r in rows],
                lambda p: zero if p[1] < 2 else k_unit)
        numeric('ABEL_'+name.upper()+'_MASS_OPERANDS', [(r['exact'], r['value'], r['residual']) for r in rows],
                lambda p: zero)
    for i, record in enumerate(result['profiles']):
        unit = record['unit']
        engine.fingerprinted(PREFIX+'_SOURCE_PROFILE_TRANSFORM_'+str(i), record['source'])
        numeric('BOUND_PROFILE_TRANSFORM_'+str(i), record['bound'], lambda p: unit)
        numeric('PROFILE_TRANSFER_'+str(i), record['transfer'], lambda p: zero)
        numeric('PROFILE_TRANSFER_INTEGRAND_'+str(i), record['normalizedIntegrand'], lambda p: unit)
        numeric('PROFILE_SAMPLE_SETTINGS_'+str(i), [(r['cutoff'], r['transfer'], r['order']) for r in record['samples']],
                lambda p: zero)
        numeric('PROFILE_SAMPLE_OPERANDS_'+str(i), [(r['adaptive'], r['gauss'], r['residual'], r['adaptiveErrorEstimate'])
                for r in record['samples']], lambda p: unit)
        numeric('PROFILE_TAIL_SETTINGS_'+str(i), [r['cutoff'] for r in record['tails']], lambda p: zero)
        numeric('PROFILE_TAIL_OPERANDS_'+str(i), [(r['left'], r['right'], r['errorEstimate']) for r in record['tails']],
                lambda p: unit)


def recover_saved(origin, base, pencil, provenance):
    """Retain all symbolic/profile work and restore the full complex mass sums."""
    import copy
    origin = origin.resolve(); origin.relative_to(STORE)
    previous = json.loads((origin/'checks.json').read_text())
    instrument = str(Path(__file__).resolve().relative_to(ROOT))
    for name, sha in previous['sourceFiles'].items():
        if digest(origin/'source'/name) != sha:
            raise ValueError(('recovery source snapshot changed', name))
        if name != instrument and digest(ROOT/name) != sha:
            raise ValueError(('recovery consumed source changed', name))
    for name, item in previous['artifacts'].items():
        if digest(origin/name) != item['sha256']:
            raise ValueError(('recovery saved operand changed', name))
    frozen = ast.parse((origin/'source'/instrument).read_text())
    current = ast.parse(Path(__file__).read_text())
    joins = {}
    for name in ('load', 'denominator_records', 'split_rule', 'adaptive_complex', 'profile_records', 'emit_result'):
        a = next(n for n in frozen.body if getattr(n, 'name', None) == name)
        b = next(n for n in current.body if getattr(n, 'name', None) == name)
        joins[name] = ast.dump(a) == ast.dump(b)
    old_abel = next(n for n in frozen.body if getattr(n, 'name', None) == 'abel_records')
    new_abel = next(n for n in current.body if getattr(n, 'name', None) == 'abel_records')
    class RestoreComplexCast(ast.NodeTransformer):
        def visit_Call(self, node):
            self.generic_visit(node)
            if (isinstance(node.func, ast.Name) and node.func.id == 'float' and len(node.args) == 1
                and isinstance(node.args[0], ast.Call) and isinstance(node.args[0].func, ast.Attribute)
                and isinstance(node.args[0].func.value, ast.Name) and node.args[0].func.value.id == 'np'
                and node.args[0].func.attr == 'dot'):
                node.func.id = 'complex'
            return node
    joins['abelRecordsComplexCastsOnly'] = ast.dump(RestoreComplexCast().visit(old_abel)) == ast.dump(new_abel)
    if not all(joins.values()):
        raise ValueError(('recovery helper AST join', joins))
    with (origin/'domains.pickle').open('rb') as stream:
        packet = pickle.load(stream)
    if packet['provenance'] != provenance:
        raise ValueError('recovery packet source identity')
    engine.PHYSICAL_METADATA.dimensions.__dict__.update(packet['dimensionState'])
    result = dict(packet['result'])
    abel = copy.deepcopy(result['abel'])
    density = sp.lambdify((abel['momentum'], abel['center'], pencil.r.regulator), abel['density'], 'numpy')
    real_differences, imaginary_parts = [], []
    for old_rows, new_rows, method in ((result['abel']['ordinary'], abel['ordinary'], 'ordinary'),
                                      (result['abel']['split'], abel['split'], 'split')):
        for old, row in zip(old_rows, new_rows):
            if method == 'ordinary':
                nodes, weights = roots_legendre(row['count'])
                nodes, weights = abel['bound']*nodes, abel['bound']*weights
            else:
                nodes, weights, panels = split_rule(-abel['bound'], abel['bound'], row['center'], row['width'], row['order'])
                if panels != row['panels'] or len(nodes) != row['count']:
                    raise ValueError('recovery quadrature panel join')
            row['value'] = complex(np.dot(weights, density(nodes, row['center'], row['regulator'])))
            row['residual'] = row['value']-row['exact']
            real_differences.append(row['value'].real-old['value'])
            imaginary_parts.append(row['value'].imag)
    result['abel'] = abel
    shutil.copyfile(origin/'denominators.pickle', base/'denominators.pickle')
    atomic_pickle(base/'abel.pickle', abel)
    record = {'runDirectory': str(origin), 'sourceFiles': previous['sourceFiles'],
        'artifacts': previous['artifacts'], 'consumedHelperAstJoins': joins,
        'realProjectionResiduals': real_differences, 'restoredImaginaryParts': imaginary_parts,
        'denominatorPacketSha256Unchanged': digest(base/'denominators.pickle') == digest(origin/'denominators.pickle'),
        'profileOperandIdentity': result['profiles'] is packet['result']['profiles']}
    save(base/'recovery.json', record)
    return result, record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-directory', type=Path, required=True)
    parser.add_argument('--resume-from', type=Path,
                        help='reuse all saved symbolic/profile results and restore full complex Abel sums')
    args = parser.parse_args()
    base = args.run_directory.resolve(); base.relative_to(STORE); base.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    paths = tuple(dict.fromkeys(SOURCES))
    pins = {str(p.relative_to(ROOT)): digest(p) for p in paths}
    for path in paths:
        target = base/'source'/path.relative_to(ROOT); target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
    checkpoint, actions, pencil, assembly, adapter, previous = load()
    provenance = {'NUMERICAL_ACTION_CHECKPOINT_SHA256': digest(CHECKPOINT),
                  'APPROVED_INPUT_SHA256': digest(native.INPUT),
                  'NUMERICAL_ACTION_PACKET_SHA256': checkpoint['artifacts']['actions.pickle']['sha256']}
    save(base/'preflight.json', {'sourceFiles': pins, 'provenance': provenance})
    recovery = None
    if args.resume_from:
        result, recovery = recover_saved(args.resume_from, base, pencil, provenance)
        denominators, pairs, abel, profiles = (result[k] for k in ('denominators', 'pairs', 'abel', 'profiles'))
    else:
        denominators, pairs = denominator_records(assembly, adapter)
        atomic_pickle(base/'denominators.pickle', {'denominators': denominators, 'pairs': pairs})
        abel = abel_records(adapter, previous)
        atomic_pickle(base/'abel.pickle', abel)
        profiles = profile_records(assembly, adapter)
        result = {'denominators': denominators, 'pairs': pairs, 'abel': abel, 'profiles': profiles}
    atomic_pickle(base/'domains.pickle', {'result': result, 'provenance': provenance,
        'dimensionState': dict(vars(engine.PHYSICAL_METADATA.dimensions))})
    before = digest(base/'domains.pickle')
    with (base/'full.out').open('x') as stream, contextlib.redirect_stdout(stream):
        emit_result(result, pencil, provenance)
        keys = {tag: 's11cd'+''.join(w.title() for w in tag.removeprefix('PY_S11CD_').split('_'))
                for tag in engine.EMISSION_LINES if not tag.startswith('PY_S11CD_METADATA_')}
        engine.physical(PREFIX+'_WRITE_KEYS', keys)
        index = engine.emission_index(engine.EMISSION_LINES)
        zero_units = {p: (0, 0, 0) for p, _ in engine.leaves(engine.cas(index))}
        engine.physical(PREFIX+'_EMISSION_LINES', index, zero_dimensions=zero_units)
    entries = {}
    for line in decoded_lines(base/'full.out'):
        tag, _, body = line.rstrip('\n').partition(': ')
        if tag in entries:
            raise ValueError('duplicate quadrature domain emission')
        entries[tag] = _restore(body)
    seen = set(); original = engine.emit
    def compare(name, value):
        tag = 'PY_S11CD_'+name
        if tag in seen or entries.get(tag) != engine.cas(value):
            raise ValueError(('quadrature domain emission mismatch', tag))
        seen.add(tag)
    engine.emit = compare
    try:
        emit_result(result, pencil, provenance)
        engine.physical(PREFIX+'_WRITE_KEYS', keys)
        engine.physical(PREFIX+'_EMISSION_LINES', index, zero_dimensions=zero_units)
    finally:
        engine.emit = original
    if seen != set(entries) or len(keys) != len(set(keys.values())) or set(keys.values()) & set(engine.IMPORT_KEYS):
        raise ValueError('quadrature domain write-key/replay census')
    final = 'PY_S11CD_'+PREFIX+'_EMISSION_LINES'
    restore_emission_index({str(k): v for k, v in entries[final]}, list(entries)[:list(entries).index(final)])
    metadata_paths = 0
    for tag, body in entries.items():
        if not tag.startswith('PY_S11CD_METADATA_'):
            continue
        for record in body:
            if len(record) == 2 and isinstance(record[0], sp.Tuple):
                fields = {str(k): v for k, v in record[1]}; count = 1
            else:
                fields = {str(k): v for k, v in record}; count = len(fields['PATHS'])
            unit = fields['DIMENSION_L_T_M']
            if not isinstance(unit, sp.Tuple) or len(unit) != 3 or any(v.free_symbols for v in unit):
                raise ValueError(('unresolved quadrature domain unit', tag))
            if any(k not in fields for k in ('MULTIGRADE', 'EPSILON_LAMBDA_SUPPORT')):
                raise ValueError(('missing quadrature domain grades', tag))
            metadata_paths += count
    summary = {'runDirectory': str(base), 'sourceFiles': pins, 'provenance': provenance,
        'denominatorCount': len(denominators), 'abelTransferPairs': [[str(s) for s in p] for p in pairs],
        'denominatorExactSignFlags': [{k: r[k] for k in ('realPositive', 'realNegative', 'normPositive')} for r in denominators],
        'primitiveResidual': str(abel['primitiveResidual']), 'halfHeightResidual': str(abel['widthResidual']),
        'widthExpression': str(abel['width']),
        'ordinaryMassErrors': {str(n): max(abs(r['residual']) for r in abel['ordinary'] if r['count'] == n)
                               for n in sorted({r['count'] for r in abel['ordinary']})},
        'splitMassErrors': {str(n): max(abs(r['residual']) for r in abel['split'] if r['order'] == n) for n in (8, 16)},
        'profileCount': len(profiles),
        'profileComparisons': [{'maxErrorsByOrder': {str(n): max(abs(r['residual']) for r in p['samples'] if r['order'] == n)
                                for n in (16, 32, 64, 128)}, 'tails': p['tails']} for p in profiles],
        'tagCount': len(entries), 'writeKeyCount': len(keys), 'metadataPaths': metadata_paths,
        'packetSha256BeforeEmission': before, 'packetSha256AfterEmission': digest(base/'domains.pickle'),
        'wallSeconds': time.monotonic()-started, 'peakRssKiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'artifacts': {p.name: {'bytes': p.stat().st_size, 'sha256': digest(p)} for p in sorted(base.iterdir())
                      if p.suffix in ('.pickle', '.out')},
        'scope': 'Source denominator operands, positive-regulator Abel resolution and finite profile-transform/tail quadrature. '
                 'Adaptive error estimates are numerical estimates. Full action/domain convergence and the Abel weak limit remain unresolved.'}
    if recovery is not None:
        summary['recovery'] = recovery
    save(base/'checks.json', summary)
    if pins != {str(p.relative_to(ROOT)): digest(p) for p in paths} or before != digest(base/'domains.pickle'):
        raise ValueError('quadrature domain source/packet changed')
    if abel['primitiveResidual'] != 0 or abel['widthResidual'] != 0:
        raise ValueError('source Abel primitive/width residual needs inspection')
    if max(abs(r['residual']) for r in abel['split'] if r['order'] == 16) > 1e-10:
        raise ValueError('split Abel mass resolution needs inspection')
    if recovery is not None and (any(v != 0 for v in recovery['realProjectionResiduals'])
            or not recovery['denominatorPacketSha256Unchanged'] or not recovery['profileOperandIdentity']):
        raise ValueError('recovered mass/source operand join needs inspection')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
