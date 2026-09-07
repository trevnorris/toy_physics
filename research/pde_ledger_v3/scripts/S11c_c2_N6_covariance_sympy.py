#!/usr/bin/env python3
"""S11c-c2 N6 source-naturality commuting-square instrument.

1. PRINT computed objects, never state conclusions.
2. PRINT the residual, never assert it (no residual-zero exit).
3. Interpretation belongs to the record.

The prediction substitutes the supplied, differentially prolonged field map
into imported inputs.mu/eps. The actual operand uses the diagnostic material
variation. Each invocation handles one case with one joint PIT, including the
source-level control deltas. A nonzero modular numerator is one-sided evidence;
all-zero samples mean only no nonzero found, conditional on the emitted delta.

Build-leg controls: in /tmp copies change ACTUAL_A_RHO to 2 for RHOBR_CONSTANT,
or ACTUAL_JUNK to 1 for MATERIAL_ADVECTED.RHO4_CONSTANT. Prediction parameters
remain 1. The shipped actual parameters are 1 and 0 respectively. No stored
outputs, cross-run PIT joins, or full-symbolic residual-zero tests are used.
"""
from __future__ import annotations

import argparse
import hashlib
import inspect
from pathlib import Path
import sys

sys.dont_write_bytecode = True
import sympy as sp
import S11c_c2_N6_reconcile_sympy as r

n = r.n
c = n.c
BASE_EMIT = n.emit
ACTUAL_A_RHO = sp.Integer(1)
ACTUAL_JUNK = sp.Integer(0)
PREDICTION_A_RHO = sp.Integer(1)
PREDICTION_H_ALPHA = sp.Integer(1)
SOURCE_OBJECTS = frozenset((
    'SOURCE_ACTUAL', 'SOURCE_PREDICTED', 'R_COV', 'SOURCE_BASELINE',
    'R_COV_BASELINE', 'SOURCE_CONTROL_DELTA', 'R_COV_CONTROL_DELTA',
))


def emit(name, payload, **labels):
    for prefix in ('N6COV_', 'N6RC_', 'N6_'):
        if name.startswith(prefix):
            name = name[len(prefix):]
            break
    if name in SOURCE_OBJECTS and 'dimension' in labels:
        labels['dimension'] = [dict(d, target=(1, -1, 0),
            consistent=not d['computed'] or d['computed'] == [(1, -1, 0)])
            for d in labels['dimension']]
    if name == 'PIT_PROVENANCE':
        per_prime = [min(sp.S.One, sp.Rational(payload['D'], p - 1 - payload['E']))
                     ** payload['draws_per_prime_cell'] for p in payload['primes']]
        payload = dict(payload, per_prime=per_prime,
            family_times_max_per_prime=payload['family_size'] * max(per_prime),
            delta_formula='min(1, family_size * max(per_prime))',
            joint_pit_calls=1)
    BASE_EMIT('N6COV_' + name, payload, **labels)


def prolonged_phi(b, mu_e, alpha, rho, kappa_a, kappa_h):
    """Discover jet paths from DERIVATIVE_MAP; differentiate only live paths.

    The independent name census also catches imported jets absent from the
    derivative graph. Such a domain error stops construction, never classifies
    a partially substituted amplitude as a covariance residual.
    """
    rho4 = b.density_pair(rho)[0]
    gradient = tuple(b.total_derivative(rho4, i, background_depth=3)
                     for i in b.DIRECTIONS)
    adv = b.dot(b.u, gradient) / rho4
    h = b.dot(b.u, b.grad_W) / b.W_bg if alpha == 'LAB_HELD' else sp.S.Zero
    base_images = {b.theta: b.theta + kappa_a * adv,
                   b.e_W: b.e_W + kappa_h * h}
    # Each entry records (base field, multi-index, parent jet, direction).
    paths = {field: (field, (), None, None) for field in base_images}
    queue = list(paths)
    for parent in queue:
        base, index, _, _ = paths[parent]
        for direction in b.DIRECTIONS:
            child = b.DERIVATIVE_MAP[direction].get(parent)
            if isinstance(child, sp.Symbol) and child not in paths:
                paths[child] = (base, (*index, direction), parent, direction)
                queue.append(child)
    candidates = sorted((atom for atom in mu_e.free_symbols
        if atom in paths or atom.name in ('theta', 'e_W')
        or atom.name.startswith(('theta_', 'grad_theta_', 'e_W_'))), key=str)
    images = dict(base_images)

    def image_of(atom):
        if atom not in images:
            _, _, parent, direction = paths[atom]
            images[atom] = b.total_derivative(image_of(parent), direction,
                                             background_depth=3)
        return images[atom]

    for atom in candidates:
        if atom in paths:
            image_of(atom)
    phi = dict(sorted(images.items(), key=lambda item: str(item[0])))
    uncovered = [atom for atom in candidates if atom not in phi]
    max_rank = max((len(paths[atom][1]) for atom in candidates if atom in paths), default=0)
    emit('PHI_DOMAIN_CENSUS', {
        'imported_mu_sha256': n.sha(mu_e),
        'imported_free_atoms': sorted(mu_e.free_symbols, key=str),
        'imported_theta_e_jets': candidates,
        'coverage': [(atom, atom in phi, paths.get(atom)) for atom in candidates],
        'uncovered': uncovered, 'domain_count': len(candidates),
        'covered_count': sum(atom in phi for atom in candidates),
        'max_present_jet_rank': max_rank, 'map_entry_count': len(phi),
        'before_amplitude_substitution': True,
        'derivative_map_sha256': n.sha(b.DERIVATIVE_MAP),
    })
    emit('FROZEN_PHI', {
        'fixed_axes': (alpha, rho), 'rho4': rho4, 'density_gradient': gradient,
        'a_rho': adv, 'h_alpha': h, 'base_images': base_images,
        'prediction_parameters': (kappa_a, kappa_h),
        'substitution_map': tuple(phi.items()),
        'jet_paths': [(atom, paths[atom]) for atom in phi],
        'background_depth': 3, 'simultaneous': True,
        'retained_grades': n.GRADES,
    })
    if uncovered:
        raise ValueError(('uncovered imported amplitude jets', uncovered))
    if max_rank + 1 > 3:
        raise ValueError(('prolongation requires greater background depth', max_rank + 1))
    return phi


def predicted_amplitude(mu_e, phi):
    """Only the imported amplitude and supplied map enter this route."""
    return [n.bounded(mu_e.subs(phi, simultaneous=True), 'predicted_mu')]


def actual_amplitudes(a, b, inputs, alpha, rho, kappa_a, kappa_j):
    """The diagnostic tag changes only actual material theta advection."""
    _, material, tag, _ = n.constitutive(a, b, inputs, alpha, rho)
    baseline = [term.subs(tag, 1) for term in material]
    actual = [term.subs(tag, kappa_a) for term in material]
    junk = sp.Symbol('n6cov_J_mu', real=True, nonzero=True)
    c.DIMENSION_SCHEMA[junk.name] = (-1, -2, 1)
    c.dimension.cache_clear()
    actual.append(kappa_j * junk * b.e_W)
    emit('ACTUAL_CONTROL_PARAMETERS', {
        'kappa_a': kappa_a, 'kappa_j': kappa_j, 'baseline_parameters': (1, 0),
        'material_tag': tag, 'junk_symbol': junk,
        'junk_dimension': c.DIMENSION_SCHEMA[junk.name],
        'junk_wave': b.e_W, 'inserted_amplitude': kappa_j * junk * b.e_W,
        'sampler_entry': (junk.name, 'global'),
        'sampler_support': '1,...,p-1',
    })
    return actual, baseline


def source_difference(left, right):
    return {face: {wave: left[face].get(wave, sp.S.Zero)
                        - right[face].get(wave, sp.S.Zero)
                   for wave in sorted(set(left[face]) | set(right[face]), key=str)}
            for face in c.FACES}


def source_columns(comp, sources, waves):
    """Retain each wave contribution separately before any kernel/extraction."""
    output = {}
    for face in c.FACES:
        for wave in waves:
            grades = n.source_value(comp, {wave: sources[face].get(wave, sp.S.Zero)})
            for grade in n.GRADES:
                output[face, str(wave), grade] = grades.get(grade, n.number(0))
    return output


def run(args):
    import S11c_a_interface_geometry_sympy_audit as a
    import S11c_b_brane_operator_sympy_audit as b
    alpha, rho = args.anchoring, args.density
    n.CASE = {'anchoring': alpha, 'density': rho}
    n.progress('imports')
    fold, _ = c.load_model(str(n.ROOT / 'scripts/S11c_b_exports.py'),
                           str(n.ROOT / 'scripts/S11c_c1_exports.py'))
    inputs = c.bind_inputs(fold)
    mu_e = inputs.mu[alpha, rho][1] / inputs.eps
    phi = prolonged_phi(b, mu_e, alpha, rho, PREDICTION_A_RHO, PREDICTION_H_ALPHA)
    mu_pred = predicted_amplitude(mu_e, phi)
    mu_actual, mu_baseline = actual_amplitudes(a, b, inputs, alpha, rho,
                                               ACTUAL_A_RHO, ACTUAL_JUNK)
    ev = {s: inputs.geometry['face_velocity'][alpha, s, 'DELTA_W'] for s in c.FACES}
    mv = r.build_material_velocity(a, alpha)
    predicted = {s: n.source_terms(inputs, alpha, rho, s, mu_pred, ev[s]) for s in c.FACES}
    actual = {s: n.source_terms(inputs, alpha, rho, s, mu_actual, mv[s]) for s in c.FACES}
    baseline = {s: n.source_terms(inputs, alpha, rho, s, mu_baseline, mv[s]) for s in c.FACES}
    rcov = source_difference(actual, predicted)
    baseline_rcov = source_difference(baseline, predicted)
    source_delta = source_difference(actual, baseline)
    waves = sorted(set().union(*(set(src[s]) for src in (actual, predicted, baseline)
                                for s in c.FACES)), key=str)
    comp = n.Compiler(inputs)
    tables = {name: source_columns(comp, src, waves) for name, src in (
        ('SOURCE_ACTUAL', actual), ('SOURCE_PREDICTED', predicted), ('R_COV', rcov),
        ('SOURCE_BASELINE', baseline), ('R_COV_BASELINE', baseline_rcov),
        ('SOURCE_CONTROL_DELTA', source_delta))}
    tables['R_COV_CONTROL_DELTA'] = n.residual(tables['R_COV'], tables['R_COV_BASELINE'])
    slots = tuple(inputs.a(prefix + label) for label in ('plus', 'minus')
                  for prefix in ('delta_p_', 'd_w_delta_p_'))
    mu_slot = sp.Symbol('n6MaterialMu', real=True)
    c.DIMENSION_SCHEMA[mu_slot.name] = (-1, -2, 1)
    m_coeff, m_rows, m_prov = r.build_material_carrier(a, b, inputs, alpha, rho, mu_slot, slots)
    kernels = {s: n.kernel_coefficients(inputs, alpha, rho, s) for s in c.FACES}
    increment = r.closed_response(comp, inputs, m_coeff, rcov, kernels)
    tables['R_COV_INCREMENT'] = {key: value for key, value in increment.items()
                                  if key[2] in (6, 9, 12)}
    emit('PROVENANCE', {
        'source_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in
            (Path(__file__), Path(r.__file__), Path(n.__file__), Path(a.__file__),
             Path(b.__file__), Path(c.__file__), n.ROOT / 'scripts/S11c_b_exports.py',
             n.ROOT / 'scripts/S11c_c1_exports.py')},
        'builder_sha256': {fn.__name__: hashlib.sha256(inspect.getsource(fn).encode()).hexdigest()
            for fn in (prolonged_phi, predicted_amplitude, actual_amplitudes, n.constitutive,
                       b.material_pullback, b.total_derivative, n.source_terms,
                       source_columns, r.build_material_velocity, r.build_material_carrier,
                       r.closed_response, n.pit, n.Sample.atom)},
        'stage_inputs': {
            'mu_E': ('inputs.mu[alpha,rho][1]', 'inputs.eps'),
            'Phi': ('supplied_N4_field_map', 'DERIVATIVE_MAP', 'total_derivative'),
            'mu_pred': ('mu_E', 'Phi', 'subs(simultaneous=True)'),
            'mu_actual': ('diagnostic.constitutive', 'actual_kappa_a', 'actual_kappa_j'),
            'SOURCE_PREDICTED': ('mu_pred', 'V_E', 'source_terms'),
            'SOURCE_ACTUAL': ('mu_actual', 'V_M', 'source_terms'),
            'R_COV': ('SOURCE_ACTUAL', 'SOURCE_PREDICTED', 'subtract'),
            'R_COV_INCREMENT': ('C_M', 'R_COV', 'reconcile.closed_response'),
            'SOURCE_CONTROL_DELTA': ('SOURCE_ACTUAL', 'SOURCE_BASELINE', 'subtract'),
            'R_COV_CONTROL_DELTA': ('R_COV', 'R_COV_BASELINE', 'subtract'),
        },
        'prediction_function_parameters': tuple(inspect.signature(predicted_amplitude).parameters),
        'prediction_function_source': inspect.getsource(predicted_amplitude),
        'operand_sha256': {name: n.sha(value) for name, value in (
            ('mu_E', mu_e), ('Phi', phi), ('mu_pred', mu_pred), ('mu_actual', mu_actual),
            ('mu_baseline', mu_baseline), ('V_E', ev), ('V_M', mv), ('C_M', m_coeff),
            ('material_rows', m_rows), ('material_faces', m_prov),
            ('SOURCE_PREDICTED', predicted), ('SOURCE_ACTUAL', actual))},
        'velocity_sha_equal': n.sha(ev) == n.sha(mv),
        'source_column_schema': ('face', 'wave', 'grade'),
        'source_column_value': 'source coefficient times recognized wave jet at Y',
        'increment_column_schema': ('block', 'grade', 'signature', 'face'),
        'increment_signatures': sorted({key[2] for key in tables['R_COV_INCREMENT']}),
        'retained_grades': n.GRADES,
    })
    objects = {('N6COV_' + name, None): table for name, table in tables.items()}
    r.emit_circuits(objects)
    n.pit(objects, inputs, comp, args.draws, args.seed)
    n.progress('case_finished')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--anchoring', choices=c.ANCHORINGS, required=True)
    parser.add_argument('--density', choices=c.DENSITIES, required=True)
    parser.add_argument('--draws', type=int, default=8)
    parser.add_argument('--seed', type=int, default=110603)
    args = parser.parse_args()
    if args.draws < 8:
        parser.error('at least eight valid draws per prime and cell')
    previous_n, previous_r = n.emit, r.emit
    n.emit = r.emit = emit
    try:
        run(args)
    except n.BlockSize:
        return 2
    finally:
        n.emit, r.emit = previous_n, previous_r
    return 0


if __name__ == '__main__':
    sys.exit(main())
