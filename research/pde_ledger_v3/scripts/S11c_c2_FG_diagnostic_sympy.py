#!/usr/bin/env python3
"""Computed F/G operands only; interpretation is external to this instrument.

The construction, closure maps, restrictions, and retained algebra belong to
S11c_c2_selfenergy_fold_sympy_audit.  This module adds source controls and a
bilinear weak-kernel reader; it never evaluates a Fourier or DtN integral.
"""
from __future__ import annotations

import copy
from functools import lru_cache
import gc
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import resource
import sys
import time

import sympy as sp
from sympy.core.function import AppliedUndef

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'scripts'))
import S11c_c2_selfenergy_fold_sympy_audit as audit
from S11c_c2_selfenergy_fold_sympy_audit import (
    bind_inputs, load_model, build_case, extract, retained_shape,
    shape_coefficients, difference, wave_jet, regression_coordinates,
)

FORWARD = 'TRANSVERSE_TO_THICKNESS'
REVERSE = 'THICKNESS_TO_TRANSVERSE'
TARGETS = ('THETA', 'E_W', 'DIV_U')
GRADES = tuple(itertools.product((0, 1), repeat=2))
PAIRING = {
    'source': ['S11c_c2_SHARED_PHYSICS/3b', 'S11c_b_SHARED_PHYSICS/1c', 'S11c_b_SHARED_PHYSICS/3c'],
    'form': 'bilinear', 'conjugation': False, 'omega_transform': 'identity',
    'i_transform': 'identity', 'Lambda_I_transform': 'identity',
    'slot_swap': ['extract.Trial', 'extract.Test'],
    'leg_swap': [['X', 'k_output'], ['Y', 'k_input']], 'face_transform': 'identity',
    'implicit_pairing_measure': 'dX', 'transposed_operator_measure': 'dY',
    'ibp': 'in_plane_interior', 'ibp_stage': 'per_grade_weak_coefficient_reader',
    'boundary_term': 0, 'temporal_ibp': False, 'Z_Fourier_integral_evaluation': False,
}


def encode(value):
    if isinstance(value, dict):
        return {str(k): encode(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, sp.Tuple)):
        return [encode(v) for v in value]
    if isinstance(value, sp.MatrixBase):
        return {'shape': list(value.shape), 'entries': encode(list(value))}
    if isinstance(value, sp.Basic):
        return sp.sstr(value)
    return value


def emit(name, value, case=None, face=None, **labels):
    print(json.dumps({'object': name, 'case': case, 'face': face,
                      'representation': 'EULERIAN', **labels,
                      'data': encode(value)}, sort_keys=True), flush=True)


def structural(name, condition, operands):
    emit('STRUCTURAL', {'name': name, 'condition': bool(condition),
                        'operands': operands})
    if not condition:
        raise ValueError(name)


def transform(value, fn):
    return audit.tree(value, fn)


def add(a, b):
    return difference(a, transform(b, lambda e: -e))


def total(values):
    values = list(values)
    result = transform(values[0], lambda e: sp.S.Zero)
    for value in values:
        result = add(result, value)
    return result


def stage(rows, inputs):
    return transform(rows, lambda e: inputs.physical_fields(retained_shape(e, inputs)))


def weak_extract(rows, inputs):
    return transform(extract(rows, inputs), lambda e: retained_shape(e, inputs))


def slots(inputs, face=None):
    return tuple(inputs.a(prefix + label)
                 for s, label in ((1, 'plus'), (-1, 'minus'))
                 if face is None or s == face
                 for prefix in ('delta_p_', 'd_w_delta_p_'))


def carrier(rows, pressure_slots):
    return difference(rows, transform(rows, lambda e: e.subs(
        {p: sp.S.Zero for p in pressure_slots}, simultaneous=True)))


def pinned(inputs, case, model, overrides=None, substitutions=None):
    # Pin to build_case's own per-face closed carriers in the baseline route.
    imported = audit.expanded_rows(inputs.slab[case])
    source = transform(imported, lambda e: e.xreplace(overrides or {}))
    chi = model['substitutions'] if substitutions is None else substitutions
    faces = {}
    for s in audit.FACES:
        raw = carrier(source, slots(inputs, s))
        sub = {p: chi[p] for p in slots(inputs, s)}
        closed_rows = (model['faces'][s] if substitutions is None else
                       stage(transform(raw, lambda e: e.subs(sub, simultaneous=True)), inputs))
        opened_rows = stage(raw, inputs)
        ic = weak_extract(closed_rows, inputs)
        ib = transform(weak_extract(opened_rows, inputs), lambda e: -e)
        faces[s] = {'S_P': raw, 'I_closed': ic, 'I_bare': ib, 'I': add(ic, ib)}
    complete = carrier(source, slots(inputs))
    ic = weak_extract(stage(transform(complete, lambda e: e.subs(chi, simultaneous=True)), inputs), inputs)
    ib = transform(weak_extract(stage(complete, inputs), inputs), lambda e: -e)
    return faces, {'S_P': complete,
                   'I_closed': total(v['I_closed'] for v in faces.values()),
                   'I_bare': total(v['I_bare'] for v in faces.values()),
                   'I': total(v['I'] for v in faces.values()),
                   'I_closed_complete': ic, 'I_bare_complete': ib,
                   'I_complete': add(ic, ib)}


def uniform_inputs(inputs, case):
    # Enumerate the imported jet ansatz, including the c1 Fourier hats.  The
    # profiles themselves are deviations from their uniform representatives.
    mapping = {inputs.eta: sp.S.Zero, inputs.sigma: sp.S.Zero,
               inputs.a('W_bg'): inputs.values['W_0'],
               inputs.a('mu_R_bg'): inputs.values['mu_R'],
               inputs.a('e_W_bg'): inputs.a('e_W')}
    for name, atom in inputs.atoms.items():
        if (name.startswith(('W_bg_d', 'mu_R_bg_d'))
                or re.fullmatch(r'[wm]1_profile(?:_d[123](?:d[123])*)?', name)
                or name == 's11cc1_w1_profile_hat_transfer'
                or name.startswith('s11cc1_w1_profile_jet_hat_')):
            mapping[atom] = sp.S.Zero
    # Obtain both density bindings from background_density_map, not a newly
    # posited density law.  Two passes resolve its W_bg representative.
    for density in inputs.density.values():
        for eq in density.atoms(sp.Equality):
            if isinstance(eq.lhs, sp.Symbol) and eq.lhs.name.startswith(('rho_4D_bg_', 'rho_br_bg_')):
                mapping[eq.lhs] = eq.rhs.xreplace(mapping)
    specialized = copy.copy(inputs)
    specialized.profiles = {k: v.xreplace(mapping)
                            for k, v in inputs.profiles.items()}
    # Preserve symbol lookup / eps,eta,sigma identities. Specialize each source
    # collection that build_face and kernel_bridge read. These are simultaneous
    # symbol substitutions (xreplace); all keys are imported coefficient symbols.
    # The specialized Inputs carries the complete map through build_case.
    for attr in ('mu', 'response', 'kernel', 'density', 'geometry'):
        setattr(specialized, attr, transform(getattr(inputs, attr),
                lambda e: e.xreplace(mapping)))
    rows = transform(audit.expanded_rows(inputs.slab[case]),
                     lambda e: e.xreplace(mapping))
    specialized = replace_rows(specialized, case, rows)
    return specialized, mapping


def inventory(mapping):
    return [{'key': k, 'value': v} for k, v in sorted(mapping.items(), key=lambda kv: str(kv[0]))]


def fingerprint(value):
    return hashlib.sha256(sp.srepr(audit.cas(value)).encode()).hexdigest()


def replace_rows(inputs, case, rows):
    result = copy.copy(inputs)
    result.slab = dict(inputs.slab)
    names = {'U_BODY_BALANCE': 'U', 'THETA_BALANCE': 'THETA', 'E_W_BALANCE': 'E_W'}
    result.slab[case] = sp.Tuple(*(sp.Tuple(name, sp.Tuple(*(
        sp.Tuple(key, rows[names[str(name)]] if str(key) == 'EXPANDED' else value)
        for key, value in entry))) if str(name) in names else sp.Tuple(name, entry)
        for name, entry in inputs.slab[case]))
    return result


def relocation(inputs, case):
    """Move the named E_W pressure coefficient to U[0], reversing its sign.

    The unit adapter is solved from the imported row dimensions and W_0's
    declared dimension.  It is not a hand-written physical coefficient.
    """
    rows = audit.expanded_rows(inputs.slab[case])
    p = inputs.a('delta_p_plus')
    coefficient = sp.expand(rows['E_W']).coeff(p)
    source_term = coefficient * p
    source_dim = audit.dimension(source_term)
    destination_dim = audit.dimension(rows['U'][0])
    length = inputs.values['W_0']
    exponent = sp.Symbol('fgUnitExponent', real=True)
    equations = [sp.Eq(exponent * d, b - a) for d, a, b in
                 zip(audit.dimension(length), source_dim, destination_dim)]
    solutions = sp.solve(equations, (exponent,), dict=True)
    structural('unit_adapter_solution', len(solutions) == 1, {'equations': equations, 'solutions': solutions})
    adapter = length ** solutions[0][exponent]
    mutated = dict(rows)
    mutated['E_W'] = rows['E_W'] - source_term
    mutated['U'] = sp.Tuple(rows['U'][0] - adapter * source_term, *rows['U'][1:])
    data = {'face': 1, 'slot': p, 'FROM': 'slab_operator/E_W_BALANCE/EXPANDED',
            'TO': 'slab_operator/U_BODY_BALANCE/EXPANDED/0',
            'coefficient': coefficient, 'source_term': source_term,
            'unit_equations': equations, 'unit_solution': solutions,
            'unit_adapter': adapter, 'adapter_domain': sp.Ne(sp.denom(adapter), 0),
            'source_dimension': source_dim, 'destination_dimension': destination_dim,
            'relocated_dimension': audit.dimension(adapter * source_term),
            'source_grades': shape_coefficients(source_term, inputs.eta, inputs.sigma),
            'destination_grades': shape_coefficients(-adapter * source_term, inputs.eta, inputs.sigma),
            'baseline_operands': rows, 'mutated_operands': mutated,
            'literal_difference': difference(mutated, rows)}
    return replace_rows(inputs, case, mutated), data


def channels(expression, inputs):
    """Linear integral kernels plus the entire non-Integral remainder.

    Bound variables are alpha-aligned to the audit's 6/9/12-variable integral
    signatures. No integrations, symmetry reductions, or momentum restrictions
    occur here. Outside factors enter the kernel before coefficient reading.
    """
    ko = tuple(inputs.a(f's11cc1_k_output_{j}') for j in range(1, 4))
    ki = tuple(inputs.a(f's11cc1_k_input_{j}') for j in range(1, 4))
    signatures = {6: (*ko, *audit.Y), 9: (*ko, *ki, *audit.Y),
                  12: (*ko, *ki, *audit.Y, *audit.MIDDLE)}

    def read(e):
        if not e.has(sp.Integral):
            return {(): e}
        if isinstance(e, sp.Integral):
            expected = signatures.get(len(e.limits))
            if expected is None:
                structural('integral_signature', False, e.limits)
            rename = dict(zip((lim[0] for lim in e.limits), expected))
            limits = tuple((v, *tuple(x.xreplace(rename) for x in lim[1:]))
                           for v, lim in zip(expected, e.limits))
            return {limits: e.function.xreplace(rename)}
        if isinstance(e, sp.Add):
            out = {}
            for term in e.args:
                for limits, kernel in read(term).items():
                    out[limits] = out.get(limits, 0) + kernel
            return out
        if isinstance(e, sp.Mul):
            out = {(): sp.S.One}
            for factor in e.args:
                nxt = {}
                for (a, v), (b, w) in itertools.product(out.items(), read(factor).items()):
                    if a and b:
                        structural('linear_integral_carrier', False, e)
                    limits = a or b
                    nxt[limits] = nxt.get(limits, 0) + v*w
                out = nxt
            return out
        structural('linear_integral_syntax', False, e)
    out = read(expression)
    out.setdefault((), sp.S.Zero)
    return out


def potential(e, kind=None):
    base = e.expr if isinstance(e, sp.Derivative) else e
    return (isinstance(base, AppliedUndef) and
            base.func.__name__.startswith('s11cc2' + (kind or '')) and
            any(base.func.__name__.startswith('s11cc2' + k) for k in
                ((kind,) if kind else ('Test', 'Trial'))))


def polynomial(expression, kind=None):
    # Sparse polynomial arithmetic only in the potential jets. In particular,
    # do not expand the physical rational coefficients to find field monomials.
    functions = tuple({a.func for a in expression.atoms(AppliedUndef) if potential(a, kind)})

    def multiply(left, right):
        out = {}
        for (a, b), (c, d) in itertools.product(left.items(), right.items()):
            basis = a*c
            out[basis] = out.get(basis, 0) + b*d
        return out

    @lru_cache(maxsize=None)
    def read(e):
        if potential(e, kind):
            return {e: sp.S.One}
        if not functions or not e.has(*functions):
            return {sp.S.One: e}
        if isinstance(e, sp.Add):
            out = {}
            for arg in e.args:
                for basis, coefficient in read(arg).items():
                    out[basis] = out.get(basis, 0) + coefficient
            return out
        if isinstance(e, sp.Mul):
            out = {sp.S.One: sp.S.One}
            for arg in e.args:
                out = multiply(out, read(arg))
            return out
        if isinstance(e, sp.Pow) and e.exp.is_Integer and e.exp >= 0:
            out = {sp.S.One: sp.S.One}
            for _ in range(int(e.exp)):
                out = multiply(out, read(e.base))
            return out
        structural('polynomial_potential_jet_syntax', False, e)
    return sorted(read(expression).items(), key=lambda kv: sp.default_sort_key(kv[0]))


@lru_cache(maxsize=4096)
def rational_coefficient(coefficient):
    # Treat Fourier phases / outgoing symbols / coefficient-field jets as
    # opaque algebraic atoms during rational cancellation. This prevents CAS
    # expansion of a Fourier exponential, and is undone before any emission.
    opaque = {a: sp.Dummy('fgRationalAtom')
              for a in coefficient.atoms(sp.Function, sp.Derivative)}
    encoded = coefficient.xreplace(opaque)
    numerator, raw_den = sp.fraction(sp.together(encoded))
    # A cleared polynomial numerator tests the rational identity without a
    # multivariate gcd or any assumptions beyond the displayed denominator.
    numerator = sp.factor_terms(sp.expand(numerator))
    denominator = sp.factor_terms(raw_den)
    reduced = numerator / denominator
    restore = {v: k for k, v in opaque.items()}
    return tuple(e.xreplace(restore) for e in (reduced, numerator, denominator, raw_den))


def interior_ibp(expression, inputs):
    """Move only spatial derivatives off tests, with boundary term fixed to 0.

    No time IBP is performed. The remaining time jets are distinct weak-basis
    entries. inputs.dx supplies the audit's local background-jet derivation.
    """
    out = []
    for basis, coefficient in polynomial(expression, 'Test'):
        if isinstance(basis, sp.Derivative) and potential(basis, 'Test'):
            spatial = [v for v in basis.variables if v in (*audit.X, *audit.Y)]
            other = [v for v in basis.variables if v not in (*audit.X, *audit.Y)]
            for v in spatial:
                coefficient = (-inputs.dx(coefficient, audit.X.index(v)) if v in audit.X
                               else -sp.diff(coefficient, v))
            basis = sp.diff(basis.expr, *other) if other else basis.expr
        out.append(basis * coefficient)
    return sp.Add(*out)


@lru_cache(maxsize=256)
def witness_table(expression, inputs):
    rows = []
    kernel_domains = []
    remainder = sp.S.Zero
    for limits, kernel in sorted(channels(expression, inputs).items(), key=lambda kv: str(kv[0])):
        bases = sorted({p.base for p in kernel.atoms(sp.Pow) if p.exp.is_negative}, key=sp.default_sort_key)
        kernel_domains.append({'integral_limits': limits, 'denominator_bases': bases,
                               'domain': sp.And(*(sp.Ne(b, 0) for b in bases))})
        normalized = interior_ibp(kernel, inputs)
        if not limits:
            remainder = normalized
        for basis, coefficient in polynomial(normalized):
            reduced, numerator, denominator, raw_den = rational_coefficient(coefficient)
            # Zero coefficients are absent entries. The explicit remainder and
            # empty rows remain data; the predicate also reads the remainder.
            if reduced == 0:
                continue
            rows.append({'integral_limits': limits, 'weak_basis': basis,
                         'coefficient': reduced, 'cleared_numerator': numerator,
                         'denominator': denominator, 'pre_cancel_denominator': raw_den,
                         'domain': sp.And(sp.Ne(raw_den, 0), sp.Ne(denominator, 0)),
                         'test_degree': sum(n for a, n in basis.as_powers_dict().items() if potential(a, 'Test')),
                         'trial_degree': sum(n for a, n in basis.as_powers_dict().items() if potential(a, 'Trial'))})
    equations = [sp.Eq(row['cleared_numerator'], 0) for row in rows]
    return {'rows': [row for row in rows if row['test_degree'] == row['trial_degree'] == 1],
            'nonbilinear_rows': [row for row in rows if not row['test_degree'] == row['trial_degree'] == 1],
            'non_integral_remainder': remainder,
            'weak_zero': sp.sympify(all(row['cleared_numerator'] == 0 for row in rows)),
            'coefficient_zero_equations': equations,
            'kernel_domains': kernel_domains,
            'denominator_domain': sp.And(*(row['domain'] for row in rows),
                                        *(row['domain'] for row in kernel_domains))}


def multigrades(expression, inputs):
    coefficients = shape_coefficients(expression, inputs.eta, inputs.sigma)
    unused = sp.Symbol('fgUnusedGrade')
    for a, b in GRADES:
        epsilon_coefficients = shape_coefficients(coefficients.get((a, b), sp.S.Zero), inputs.eps, unused)
        for epsilon in (0, 1):
            yield {'epsilon': epsilon, 'eta': a, 'sigma_W': b}, epsilon_coefficients.get((epsilon, 0), sp.S.Zero)


def emit_blocks(name, blocks, inputs, case, face):
    checks = []
    for direction, target in itertools.product((FORWARD, REVERSE), TARGETS):
        expression = blocks[direction][target]
        emit(name + '_OPERAND', expression, case, face, block=direction + '/' + target)
        for grade, coefficient in multigrades(expression, inputs):
            table = witness_table(coefficient, inputs)
            emit(name + '_WEAK_TABLE', table, case, face,
                 block=direction + '/' + target, grade=grade)
            checks.append({'block': direction + '/' + target, 'grade': grade,
                           'result': table['weak_zero']})
    emit(name + '_WEAK_ZERO', {'result': sp.And(*(row['result'] for row in checks)),
                             'by_block_grade': checks}, case, face)


def formal_adjoint(expression, inputs):
    """Bilinear slot/leg transpose of the independently extracted weak form.

    X has an implicit in-plane pairing measure. Transposition exchanges that
    measure with explicit Y; the returned operator is again represented with
    external X and integrated Y. Face, frequency and complex scalars persist.
    """
    swap = dict(zip(audit.X, audit.Y)) | dict(zip(audit.Y, audit.X))
    for suffix in ('1', '2', '3'):
        ko, ki = (inputs.a('s11cc1_k_' + leg + '_' + suffix) for leg in ('output', 'input'))
        swap.update({ko: ki, ki: ko})
    qo, qi = (inputs.a('s11cc1_q_out_' + leg) for leg in ('output', 'input'))
    swap.update({qo: qi, qi: qo})
    result = []
    for limits, kernel in channels(expression, inputs).items():
        # Make the location of every coefficient explicit before swapping legs.
        local = {a: inputs.at_source(a).xreplace(dict(zip(audit.Y, audit.X)))
                 for a in kernel.atoms(sp.Symbol)
                 if re.match(r'(?:W_bg|mu_R_bg|[wm]1_profile)(?:$|_d)', a.name)}
        kernel = kernel.xreplace(local)
        kernel = kernel.replace(lambda a: potential(a) and isinstance(a, AppliedUndef),
            lambda a: sp.Function(a.func.__name__.replace('s11cc2Trial', 'fgSwapTest')
                .replace('s11cc2Test', 's11cc2Trial').replace('fgSwapTest', 's11cc2Test'))(*a.args))
        kernel = kernel.xreplace(swap)
        if not limits:
            # The local pairing has only one spatial integration variable;
            # rename its swapped dummy back to the displayed external X.
            kernel = kernel.xreplace(dict(zip(audit.Y, audit.X)))
        if len(limits) == 6:
            # The sole diagonal momentum is a bound dummy, not k_in=k_out.
            rename = {inputs.a(f's11cc1_k_input_{j}'): inputs.a(f's11cc1_k_output_{j}') for j in range(1, 4)}
            kernel = kernel.xreplace(rename)
        # Keep the transposed weak form here. The common coefficient reader
        # performs its in-plane IBP after grade extraction, avoiding expansion
        # of all grades together. This is the same bilinear weak operator.
        result.append(sp.Integral(kernel, *limits) if limits and kernel != 0 else kernel)
    return sp.Add(*result)


def slot_inventory(name, obj, inputs, case, face):
    symbols = audit.cas(obj).free_symbols
    emit(name, {'free_symbols': sorted(symbols, key=str),
                'P': [{'slot': p, 'present': p in symbols} for p in slots(inputs)]}, case, face)


def ordering(inputs, case, model, mapping):
    imported = audit.named(inputs.open_kernel[case], 'COMPLETE_OPERATOR_BLOCKS')
    closed = transform(imported, lambda e: inputs.physical_fields(
        e.xreplace(mapping).subs(model['substitutions'], simultaneous=True)))
    aligned = regression_coordinates(closed, inputs)
    emit('ORDERING_EXTRACT_FIRST_UNIFORM_OPERAND', aligned, case, 'assembled')
    emit('CLOSED_COUPLING_KERNEL_UNIFORM', model['closed_kernel'], case, 'assembled')
    emit('ORDERING_COMMUTATOR_UNIFORM', difference(model['closed_kernel'], aligned), case, 'assembled')


def f_diagnostic(inputs, case):
    uniform, mapping = uniform_inputs(inputs, case)
    emit('F_UNIFORM_MAP', inventory(mapping), case)
    emit('F_UNIFORM_MAP_TRANSPORT', {'operation': 'simultaneous_symbol_xreplace',
         'route': 'specialized_Inputs', 'collections': ['slab.EXPANDED', 'mu', 'response',
             'kernel', 'density', 'geometry', 'profiles'],
         'consumers': ['build_case', 'build_face', 'kernel_bridge', 'regression_coordinates']}, case)
    emit('F_UNIFORM_DENSITY_OPERANDS', {'baseline': inputs.density, 'specialized': uniform.density}, case)
    emit('F_UNIFORM_FACE_NORMAL_OPERANDS', {'baseline': inputs.geometry['face_normal'],
                                         'specialized': uniform.geometry['face_normal']}, case)
    model = build_case(uniform, *case)
    emit('F_UNIFORM_CHI', inventory(model['substitutions']), case)
    faces, assembled = pinned(uniform, case, model, mapping)
    for face, values in [*faces.items(), ('assembled', assembled)]:
        emit_blocks('F_I_CLOSED_UNIFORM', values['I_closed'], uniform, case, face)
    ordering(uniform, case, model, mapping)

    # Source-level polynomial ansatz in the pressure closure map. extract itself
    # supplies the independent test; wave_jet supplies the unrestricted trial.
    p = slots(uniform, 1)[0]
    amplitude = sp.Symbol('fgInjectedPressurePerDisplacement')
    displacement = wave_jet(inputs.a('u_1'))
    audit.NEW_DIMENSIONS[amplitude] = tuple(a-b for a, b in zip(audit.dimension(p), audit.dimension(displacement)))
    injected = dict(model['substitutions'])
    injected[p] += amplitude * displacement
    emit('F_INJECTED_CONTROL_SOURCE', {'slot': p, 'amplitude_dimension': audit.dimension(amplitude),
         'baseline_operand': model['substitutions'][p], 'mutated_operand': injected[p],
         'literal_difference': difference(injected[p], model['substitutions'][p])}, case, 1,
         control='injected')
    injected_faces, injected_sum = pinned(uniform, case, model, mapping, injected)
    for face, values in [*injected_faces.items(), ('assembled', injected_sum)]:
        emit_blocks('F_INJECTED_CONTROL_I_CLOSED_UNIFORM', values['I_closed'], uniform, case, face)

    changed_inputs, operands = relocation(inputs, case)
    emit('F_FORM_SOURCE', operands, case, 1, control='injected')
    changed_uniform, changed_mapping = uniform_inputs(changed_inputs, case)
    changed_model = build_case(changed_uniform, *case)
    changed_faces, changed_sum = pinned(changed_uniform, case, changed_model, changed_mapping)
    for face, values in [*changed_faces.items(), ('assembled', changed_sum)]:
        baseline = assembled if face == 'assembled' else faces[face]
        emit('F_FORM_I_CLOSED_UNIFORM_LITERAL_DIFFERENCE', difference(values['I_closed'], baseline['I_closed']), case, face)
        emit_blocks('F_FORM_I_CLOSED_UNIFORM', values['I_closed'], changed_uniform, case, face)


def g_diagnostic(inputs, case, baseline_model, faces, assembled):
    reverse_source_before = fingerprint(inputs.slab[case])
    for face, values in [*faces.items(), ('assembled', assembled)]:
        emit_blocks('G_INCREMENT', values['I'], inputs, case, face)
    emit('G_DIRECT_INCREMENT', difference(baseline_model['closed_kernel'], baseline_model['open_kernel']), case, 'assembled')

    # Independent source routes: direct reverse from close/extract minus open
    # extract; forward from a separately constructed per-face pressure carrier.
    # Neither block is manufactured from the other; the mutation has one owner.
    forward_inputs = copy.copy(inputs)
    forward_inputs.slab = dict(inputs.slab)
    forward_source_before = fingerprint(forward_inputs.slab[case])
    forward_model = build_case(forward_inputs, *case)
    forward_faces, forward_sum = pinned(forward_inputs, case, forward_model)
    changed_inputs, operands = relocation(inputs, case)
    emit('G_ONE_SIDED_SOURCE', operands, case, 1, control='injected')
    changed_model = build_case(changed_inputs, *case)
    changed_faces, changed_sum = pinned(changed_inputs, case, changed_model)
    # Recompute the unmodified route from its own source copy, including extract.
    forward_again = build_case(forward_inputs, *case)
    repeated_faces, repeated_sum = pinned(forward_inputs, case, forward_again)
    emit('G_ROUTE_PROVENANCE_AVAILABILITY', {
        'direct_reverse': {'source': 'slab_operator/U_BODY_BALANCE/EXPANDED',
            'route': ['build_case.closed_kernel', 'build_case.open_kernel', 'difference', REVERSE],
            'baseline_hash': reverse_source_before, 'corrupted_hash': fingerprint(changed_inputs.slab[case]),
            'available': REVERSE in baseline_model['closed_kernel']},
        'forward': {'source': 'independent_copy/slab_operator',
            'route': ['build_case.substitutions', 'S_P_face', 'extract_closed_plus_bare', FORWARD, 'bilinear_slot_leg_swap'],
            'baseline_hash': forward_source_before, 'recomputed_hash': fingerprint(forward_inputs.slab[case]),
            'available': FORWARD in forward_sum['I']},
        'separate_input_objects': inputs is not forward_inputs and changed_inputs is not forward_inputs,
        'separate_slab_maps': inputs.slab is not forward_inputs.slab and changed_inputs.slab is not forward_inputs.slab,
        'mutation_owner': 'direct_reverse', 'face_preserved': 1}, case)
    for face in (*audit.FACES, 'assembled'):
        reverse_base = (difference(baseline_model['closed_kernel'], baseline_model['open_kernel'])
                        if face == 'assembled' else faces[face]['I'])
        reverse_changed = (difference(changed_model['closed_kernel'], changed_model['open_kernel'])
                           if face == 'assembled' else changed_faces[face]['I'])
        fwd = forward_sum if face == 'assembled' else forward_faces[face]
        fwd_repeat = repeated_sum if face == 'assembled' else repeated_faces[face]
        for target in TARGETS:
            adj = formal_adjoint(fwd['I'][FORWARD][target], inputs)
            adj_repeat = formal_adjoint(fwd_repeat['I'][FORWARD][target], inputs)
            direct = reverse_base[REVERSE][target]
            changed = reverse_changed[REVERSE][target]
            objects = {'convention_table': PAIRING,
                       'direct_reverse': direct, 'independent_forward': fwd['I'][FORWARD][target],
                       'formal_adjoint': adj, 'pairing_residual': difference(direct, adj),
                       'corrupted_direct_reverse': changed,
                       'recomputed_independent_forward': fwd_repeat['I'][FORWARD][target],
                       'recomputed_forward_adjoint': adj_repeat,
                       'corrupted_pairing_residual': difference(changed, adj_repeat),
                       'direct_route_literal_difference': difference(changed, direct),
                       'forward_route_literal_difference': difference(adj_repeat, adj)}
            emit('G_PAIRING_OPERANDS', objects, case, face, target=target)
            for name in ('pairing_residual', 'corrupted_pairing_residual', 'direct_route_literal_difference', 'forward_route_literal_difference'):
                for grade, coefficient in multigrades(objects[name], inputs):
                    emit('G_PAIRING_WEAK_TABLE', witness_table(coefficient, inputs),
                         case, face, target=target, operand=name, grade=grade)


def diagnostic():
    started = time.monotonic()
    emit('METADATA', {
        'mode': os.environ.get('S11CC2_FG', 'FULL').upper(),
        'construction': str(ROOT / 'scripts/S11c_c2_selfenergy_fold_sympy_audit.py'),
        'construction_sha256': hashlib.sha256(Path(audit.__file__).read_bytes()).hexdigest(),
        'scope': {'F': 'secondary_smoke_test', 'validates_nonuniform_coefficient': False,
                  'validates_nonuniform_sign': False, 'validates_nonuniform_parity': False,
                  'G_representation': 'EULERIAN', 'cross_representation_comparison': 'N6_pending'},
        'retained_eta_sigma_grades': GRADES,
        'pairing': PAIRING,
        'weak_basis': {'trial_test': 'generic_independent_compact_support_potentials',
                       'affine_slot_terms': 'trial_degree_0', 'zero_rows': 'empty_table'},
    })
    fold, _ = load_model(str(ROOT / 'scripts/S11c_b_exports.py'), str(ROOT / 'scripts/S11c_c1_exports.py'))
    inputs = bind_inputs(fold)
    del fold
    symbols = set(inputs.atoms.values()) | set((*audit.X, *audit.Y, audit.TIME, *audit.MIDDLE, audit.MIDDLE_Q))
    emit('SYMBOL_ASSUMPTIONS', [{'symbol': a, 'assumptions': a.assumptions0}
                                for a in sorted(symbols, key=str)])
    selected = list(itertools.product(audit.ANCHORINGS, audit.DENSITIES))
    if os.environ.get('S11CC2_FG', '').upper() == 'TRIAGE':
        selected = selected[:1]
    for case in selected:
        emit('CASE_BEGIN', {'elapsed_seconds': time.monotonic() - started}, case)
        model = build_case(inputs, *case)
        faces, assembled = pinned(inputs, case, model)
        direct = difference(model['closed_kernel'], model['open_kernel'])
        emit('F_I_DIRECT', direct, case, 'assembled')
        emit('F_SPLIT_OPERANDS', assembled, case, 'assembled')
        emit('F_R_SPLIT', difference(direct, assembled['I_complete']), case, 'assembled', role='transcription_linearity')
        emit('F_R_ADD', difference(direct, total(v['I'] for v in faces.values())), case, 'assembled')
        for face, values in [*faces.items(), ('assembled', assembled)]:
            emit('F_I_BARE', values['I_bare'], case, face)
            slot_inventory('F_I_CLOSED_SLOT_INVENTORY', values['I_closed'], inputs, case, face)
            slot_inventory('F_I_BARE_SLOT_INVENTORY', values['I_bare'], inputs, case, face)
        f_diagnostic(inputs, case)
        g_diagnostic(inputs, case, model, faces, assembled)
        emit('CASE_END', {'elapsed_seconds': time.monotonic() - started,
                          'max_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}, case)
        del model, faces, assembled, direct
        # The audit memoizes large CAS operands; process cases serially and
        # release those caches only after all routes for this case are emitted.
        for fn in (retained_shape, shape_coefficients, audit.restricted,
                   audit.outgoing_spectral, audit.Inputs.dx, rational_coefficient, witness_table):
            fn.cache_clear()
        gc.collect()
    emit('COMPUTED_BINDINGS', inventory(audit.COMPUTED_BINDINGS))
    emit('END', {'cases': selected, 'elapsed_seconds': time.monotonic() - started,
                 'max_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})


if __name__ == '__main__':
    # build_case owns a legacy progress file. Restore it after this diagnostic;
    # its writes are not a new diagnostic artifact or an adjudication record.
    progress = ROOT / '_measurements/S11c_c2_sympy_progress.json'
    before = progress.read_bytes() if progress.exists() else None
    try:
        diagnostic()
    finally:
        if before is None:
            progress.unlink(missing_ok=True)
        else:
            progress.write_bytes(before)
