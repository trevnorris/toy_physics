#!/usr/bin/env python3
"""S11c-c2 N6 fixed-anchoring reconcile instrument.

1. The script may PRINT computed objects. It may NOT state conclusions.
2. PRINT the residual; do NOT assert it.
3. Interpretation belongs to the reconcile record.

Native builders and the imported Eulerian substrate supply the two operands.
The source is the diagnostic's imported-slot source_terms circuit. The affine
split is I(dC, ms) + B(m_coeff, ds) + B(dC, ds); B visits only 6/9/12.
One invocation handles one case, with one joint PIT. No review/control runs or
stored PIT inputs are used. Builder entry points remain separate for external
one-sided ablations in temporary copies.
"""
from __future__ import annotations

import argparse
import hashlib
import inspect
import itertools
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
import sympy as sp
import S11c_c2_N6_diagnostic_sympy as n

c = n.c
ROOT = n.ROOT
BASE_EMIT = n.emit
SOURCE_NAMES = frozenset('N6RC_SOURCE_' + suffix for suffix in
                         ('EULERIAN', 'MATERIAL', 'BRIDGE_RESIDUAL'))


def emit(name, payload, **labels):
    """Namespace imported metadata; supply the combined-source target units."""
    if name in SOURCE_NAMES and 'dimension' in labels:
        labels['dimension'] = [dict(d, target=(1, -1, 0),
            consistent=not d['computed'] or d['computed'] == [(1, -1, 0)])
            for d in labels['dimension']]
    BASE_EMIT(name if name.startswith('N6RC_') else 'N6RC_' + name,
              payload, **labels)


def frozen_relations(a, b, inputs, alpha, rho, slots):
    """Predeclare the permitted maps, before either residual is constructed."""
    rho4, rhobr, _ = b.density_pair(rho)
    gradient = tuple(b.total_derivative(rho4, i, background_depth=3) for i in range(3))
    g = tuple(x / rho4 for x in gradient)
    adv = b.dot(b.u, g)
    h = b.dot(b.u, b.grad_W) / b.W_bg if alpha == 'LAB_HELD' else sp.S.Zero
    maps = {b.theta: b.theta + adv, b.e_W: b.e_W + h}
    maps.update({b.grad_theta[i]: b.total_derivative(b.theta + adv, i, background_depth=3)
                 for i in range(3)})
    maps.update({b.grad_e[i]: b.total_derivative(b.e_W + h, i, background_depth=3)
                 for i in range(3)})
    parameter = sp.Symbol('n6rc_map_parameter', real=True)
    jet_bridge = tuple((a.grad_theta[i], b.grad_theta[i]) for i in range(3))
    emit('N6RC_FROZEN_RELATIONS', {
        'fixed_axes': (alpha, rho), 'rho4': rho4, 'rho_br_bg': rhobr,
        'g': g, 'a_rho': adv, 'h_alpha': h,
        'field_and_derivative_maps': tuple(maps.items()),
        'quadratic_density_jacobian': 1 + b.trace(b.grad_u),
        'density_wave_projection_degree': 2,
        'material_inverse_transpose': a.material_inverse_transpose(parameter, a.grad_u),
        'map_parameter': parameter,
        'jet_vocabulary_bridge': jet_bridge,
        'jet_physical_field_images': [(c.wave_jet(x), c.wave_jet(y)) for x, y in jet_bridge],
        'pressure_identities': [(p, a.dw_delta_p[s] if p.name.startswith('d_w_')
                                else (a.delta_p_plus if s == 1 else a.delta_p_minus))
                               for p in slots for s in c.FACES
                               if p.name.endswith('plus' if s == 1 else 'minus')],
        'normalizations': ('retained_eta_sigma_rectangle', 'compact_support_interior_spatial_ibp',
                           'alpha_aligned_formal_kernel_signatures', 'inherited_on_shell_charts',
                           'inherited_profile_mixed_partials', 'inherited_fourier_profile_jet_identity'),
        'retained_grades': n.GRADES,
        'strong_jet_depth': 3, 'weak_jet_depth': 4,
    })


def build_material_velocity(a, alpha):
    # This source stage precedes and is independent of the carrier stage below.
    return {s: a.finalize(a.face_velocity_raw(a.build_face_source(
        alpha, s, 'DELTA_W', 'RHO4_CONSTANT', route='MATERIAL'))) for s in c.FACES}


def build_material_carrier(a, b, inputs, alpha, rho, mu_slot, slots):
    # A temporary-copy carrier-map ablation belongs inside this function; its
    # returned velocity is deliberately unused, so it cannot rebind source V.
    rows, _, provenance = n.face_factory(a, b, inputs, alpha, rho, 'MATERIAL', mu_slot)
    return n.pressure_coefficients(rows, slots), rows, provenance


def coefficient_table(comp, coeffs):
    result = {}
    for (row, p), value in coeffs.items():
        grades = comp.grades(value)
        face = 1 if p.name.endswith('plus') else -1
        for grade in n.GRADES:
            result[row, str(p), face, grade] = grades.get(grade, n.number(0))
    return result


def source_table(comp, sources):
    result = {}
    for face, terms in sources.items():
        grades = n.source_value(comp, terms)
        for grade in n.GRADES:
            result[face, grade] = grades.get(grade, n.number(0))
    return result


def closed_response(comp, inputs, coeffs, sources, kernels):
    """B(C,S): the imported template/kernel path, signatures 6/9/12 only."""
    output = {}
    for (row, p), coefficient in coeffs.items():
        coefficient = sp.cancel(coefficient / inputs.eps)
        face = 1 if p.name.endswith('plus') else -1
        jet = p.name.startswith('d_w_')
        if coefficient == 0:
            continue
        cg = [comp.grades(coefficient)] + [comp.grades(inputs.dx(coefficient, i)) for i in range(3)]
        for signature in (6, 9, 12):
            kernel = kernels[face][0][signature] * (kernels[face][1] if jet else 1)
            kg = comp.grades(kernel)
            for wave, source_coefficient in sources[face].items():
                sg = comp.grades(source_coefficient, 'Y')
                cf, weak = n.template(row, wave, signature, inputs)
                for cg_grade in n.GRADES:
                    bindings = {cf: cg[0].get(cg_grade, n.number(0))}
                    bindings.update({sp.diff(cf, c.X[i]): cg[i+1].get(cg_grade, n.number(0))
                                     for i in range(3)})
                    for block, expression in weak.items():
                        w = comp.scalar(expression, binding=bindings)
                        products = n.gmul({cg_grade: w}, n.gmul(sg, kg))
                        for grade, value in products.items():
                            key = (block, grade, signature, face)
                            output[key] = n.plus(output.get(key, n.number(0)),
                                n.times(n.times(value, n.measure(signature)), n.amplitude()))
        n.progress('closed_response_slot', row=row, slot=str(p))
    # An empty signature-0 column is a zero contraction over an empty set.
    for block, grade, signature in itertools.product(n.BLOCKS, n.GRADES, (0, 6, 9, 12)):
        for face in c.FACES:
            output.setdefault((block, grade, signature, face), n.number(0))
        output[block, grade, signature, 0] = n.plus(output[block, grade, signature, 1],
                                                               output[block, grade, signature, -1])
    return output


def emit_circuits(objects):
    """Print the actual arithmetic DAG, using references instead of expansion."""
    indices = {}
    nodes = []
    roots = {}
    for (name, _), table in objects.items():
        entries = []
        for key in sorted(table, key=str):
            root = table[key]
            stack = [root]
            while stack:
                node = stack[-1]
                if node in indices:
                    stack.pop()
                    continue
                missing = [x for x in node.args if isinstance(x, n.Node) and x not in indices]
                if missing:
                    stack.extend(missing)
                    continue
                args = [{'ref': indices[x]} if isinstance(x, n.Node) else {'literal': x}
                        for x in node.args]
                indices[node] = len(nodes)
                nodes.append({'op': node.op, 'args': args, 'degree': node.degree})
                stack.pop()
            entries.append((key, indices[root]))
        roots[name] = entries
    emit('N6RC_ARITHMETIC_DAG', {'nodes': nodes})
    for (name, _), table in objects.items():
        emit(name + '_NODES', {'columns': [key for key, _ in roots[name]],
             'root_ids': [index for _, index in roots[name]],
             'root_nodes': [nodes[index] for _, index in roots[name]]})


def run(args):
    import S11c_a_interface_geometry_sympy_audit as a
    import S11c_b_brane_operator_sympy_audit as b
    n.CASE = {'anchoring': args.anchoring, 'density': args.density}
    alpha, rho = args.anchoring, args.density
    n.progress('imports')
    fold, _ = c.load_model(str(ROOT / 'scripts/S11c_b_exports.py'), str(ROOT / 'scripts/S11c_c1_exports.py'))
    inputs = c.bind_inputs(fold)
    slots = tuple(inputs.a(prefix + label) for label in ('plus', 'minus')
                  for prefix in ('delta_p_', 'd_w_delta_p_'))
    frozen_relations(a, b, inputs, alpha, rho, slots)
    comp = n.Compiler(inputs)
    imported = n.flatten(c.expanded_rows(inputs.slab[alpha, rho]))
    e_coeff = n.pressure_coefficients(imported, slots)
    _, material, t, adv = n.constitutive(a, b, inputs, alpha, rho)
    mu_e = [inputs.mu[alpha, rho][1] / inputs.eps]
    mu_m = [e.subs(t, 1) for e in material]
    mu_slot = sp.Symbol('n6MaterialMu', real=True)
    c.DIMENSION_SCHEMA[mu_slot.name] = (-1, -2, 1)
    n.progress('native_faces')
    m_v = build_material_velocity(a, alpha)
    m_coeff, m_rows, m_prov = build_material_carrier(a, b, inputs, alpha, rho, mu_slot, slots)
    ev = {s: inputs.geometry['face_velocity'][alpha, s, 'DELTA_W'] for s in c.FACES}
    es = {s: n.source_terms(inputs, alpha, rho, s, mu_e, ev[s]) for s in c.FACES}
    ms = {s: n.source_terms(inputs, alpha, rho, s, mu_m, m_v[s]) for s in c.FACES}
    # Differences originate at the two interface operands, never from R_N6.
    dc = {key: e_coeff[key] - m_coeff[key] for key in e_coeff}
    ds = {s: {w: es[s].get(w, sp.S.Zero) - ms[s].get(w, sp.S.Zero)
              for w in sorted(set(es[s]) | set(ms[s]), key=str)} for s in c.FACES}
    emit('N6RC_PROVENANCE', {
        'source_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in
            (Path(__file__), Path(n.__file__), Path(a.__file__), Path(b.__file__), Path(c.__file__),
             ROOT / 'scripts/S11c_b_exports.py', ROOT / 'scripts/S11c_c1_exports.py',
             ROOT / '_measurements/S11c_c2_N6_route2_spec_astra.md',
             ROOT / 'directives/S11c_c2_SHARED_PHYSICS.md')},
        'operand_fingerprints': {name: n.sha(value) for name, value in
            (('C_E', e_coeff), ('C_M', m_coeff), ('es', es), ('ms', ms),
             ('imported_rows', imported), ('material_rows', m_rows), ('material_faces', m_prov),
             ('mu_E', mu_e), ('mu_M', mu_m), ('V_E', ev), ('V_M', m_v))},
        'builder_fingerprints': {fn.__name__: hashlib.sha256(inspect.getsource(fn).encode()).hexdigest()
            for fn in (n.constitutive, n.face_factory, n.source_terms, n.build_increment,
                       n.template, n.kernel_coefficients, n.pit, closed_response,
                       build_material_carrier, build_material_velocity)},
        'stage_inputs': {'C_E': ('inputs.slab', 'expanded_rows', 'pressure_coefficients'),
                         'C_M': ('MATERIAL', 'face_factory', 'pressure_coefficients'),
                         'es': ('inputs.mu/eps', "inputs.geometry[face_velocity]", 'source_terms'),
                         'ms': ('constitutive(t=1)', 'build_material_velocity', 'source_terms'),
                         'CARRIER_CHANNEL': ('C_E-C_M', 'ms', 'build_increment'),
                         'SOURCE_CHANNEL': ('C_M', 'es-ms', 'closed_response'),
                         'CROSS_CHANNEL': ('C_E-C_M', 'es-ms', 'closed_response'),
                         'R_N6': ('build_increment(C_E,es)', 'build_increment(C_M,ms)', 'residual')},
        'carrier_dependencies': {route: [(row, str(p), e.has(mu_slot), e.has(t), e.has(b.theta),
                                          tuple(e.has(u) for u in b.u))
                                        for (row, p), e in coeffs.items()]
                                 for route, coeffs in (('EULERIAN', e_coeff), ('MATERIAL', m_coeff))},
    })
    if rho == 'RHO4_CONSTANT':
        emit('N6RC_ADVECTION_ABSENCE', {'a_rho': adv,
             'density_gradient': tuple(b.total_derivative(b.density_pair(rho)[0], i,
                                       background_depth=3) for i in range(3)),
             'material_mu_tag_derivatives': [sp.diff(e, t) for e in material]})

    objects = {}
    def add(name, table):
        objects['N6RC_' + name, None] = table

    ce, cm = coefficient_table(comp, e_coeff), coefficient_table(comp, m_coeff)
    se, sm = source_table(comp, es), source_table(comp, ms)
    add('CARRIER_EULERIAN', ce)
    add('CARRIER_MATERIAL', cm)
    add('CARRIER_BRIDGE_RESIDUAL', n.residual(ce, cm))
    add('SOURCE_EULERIAN', se)
    add('SOURCE_MATERIAL', sm)
    add('SOURCE_BRIDGE_RESIDUAL', n.residual(se, sm))
    kernels = {s: n.kernel_coefficients(inputs, alpha, rho, s) for s in c.FACES}
    carrier, _ = n.build_increment(comp, inputs, dc, ms, kernels, slots)
    source = closed_response(comp, inputs, m_coeff, ds, kernels)
    cross = closed_response(comp, inputs, dc, ds, kernels)
    add('CARRIER_CHANNEL', carrier)
    add('SOURCE_CHANNEL', source)
    add('CROSS_CHANNEL', cross)
    E, edim = n.build_increment(comp, inputs, e_coeff, es, kernels, slots)
    M, mdim = n.build_increment(comp, inputs, m_coeff, ms, kernels, slots)
    R = n.residual(E, M)
    add('EULERIAN_OPERAND', E)
    add('MATERIAL_OPERAND', M)
    add('R_N6', R)
    split_sum = {k: n.nsum((carrier[k], source[k], cross[k])) for k in R}
    add('SPLIT_SUM', split_sum)
    add('SPLIT_CHECK', n.residual(split_sum, R))
    emit('N6RC_DIMENSIONS', {'eulerian': edim, 'material': mdim})
    emit_circuits(objects)
    n.pit(objects, inputs, comp, args.draws, args.seed)
    n.progress('case_finished')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--anchoring', choices=c.ANCHORINGS, required=True)
    parser.add_argument('--density', choices=c.DENSITIES, required=True)
    parser.add_argument('--draws', type=int, default=8)
    parser.add_argument('--seed', type=int, default=110602)
    args = parser.parse_args()
    if args.draws < 8:
        parser.error('at least eight valid draws per prime and cell')
    n.emit = emit
    try:
        run(args)
    except n.BlockSize:
        return 2
    finally:
        n.emit = BASE_EMIT
    return 0


if __name__ == '__main__':
    sys.exit(main())
