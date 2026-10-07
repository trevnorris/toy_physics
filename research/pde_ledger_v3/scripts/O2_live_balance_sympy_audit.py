#!/usr/bin/env python3
"""O2 conditional live balance; O2_SHARED_PHYSICS.md is the physics input.

OPEN_* applications are unevaluated functional ACTIONS. Their displayed tuple
is an operand/obligation register, NOT an exhaustive constitutive argument list.
The OpenDomain token admits arbitrary fields, jets, histories and nonlocality.
Derivatives of these actions are held: expanding them into a finite chain rule
would select a local, finite-state realization that the input does not supply.
Derivatives of the supplied radial profiles and graph are computed explicitly.

Signs: storage + outward spatial transport + outward exchange - applied forces.
The energy power-accounting action retains unresolved overlap between storage,
transport, stress work, relaxation, conversion, boundary work and supply. It is
not an additive constitutive decomposition of the named energy operands.

No controls or expected physical outcomes live here. The ablation harness
changes one marked construction line in an in-memory source copy per knife.
"""
from __future__ import annotations

import hashlib
from functools import lru_cache
from sympy.core.function import AppliedUndef
from sympy.printing.repr import ReprPrinter
from pathlib import Path
from types import MappingProxyType
import sys

import sympy as sp
from sympy.core.symbol import Str
from sympy.functions.elementary.piecewise import ExprCondPair

HERE = Path(__file__).resolve()
ROOT = HERE.parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / 'scripts'))
from ledger_fold import (load_model, check_consumer, assert_lookups_equal_manifest,
                         assert_delta_is_minimal, _RELATIONALS)

# No supplied LIVE object is identified with a computed row of this historical
# fold. In particular, its uniform transverse anchor is not mu_perp(r).
IMPORT_KEYS = ()
SPEC = ROOT / 'directives/O2_SHARED_PHYSICS.md'
FOLD_PATHS = tuple(ROOT / 'scripts' / p for p in (
    'S11c_b_exports.py', 'S11c_c1_exports.py', 'S11c_c2_exports.py'))
SECTION9 = (
    'BASIS', 'MEASURES', 'GEOMETRY', 'PROFILES', 'MATERIAL_VELOCITY',
    'MATERIAL_INPUT_DIFFERENTIALS', 'MOMENTUM_DENSITY', 'MOMENTUM_FLUX',
    'MOMENTUM_STORAGE', 'MOMENTUM_TRANSPORT', 'MATERIAL_COMPATIBILITY',
    'INTERNAL_FORCE', 'NATIVE_BULK_LOAD', 'NATIVE_HOLD_LOAD', 'MECHANICAL_LOAD',
    'CARRIED_MOMENTUM', 'SOURCE_PARTNERS', 'DRIVE', 'HOLD_INPLANE', 'HOLD_W',
    'HOLD_GRAPH_NORMAL', 'MASS_INPUT', 'MASS_RESIDUAL', 'MATERIAL_POWER_PAIRING',
    'FACE_POWER_PAIRING', 'ENERGY_STORAGE', 'ENERGY_TRANSPORT', 'ENERGY_POWER',
    'ENERGY_STEADY', 'COUPLED_INPUTS', 'MODEL_POINT', 'TRACE',
)
LOCAL_NAMES = ('IMPORT_MANIFEST', 'FOLD', 'F9', 'F3', 'F6', 'ROUNDTRIP',
               'EXPORT_CLOSURE', 'INPUT_DIGESTS', 'NAME_AUDIT')


class LosslessReprPrinter(ReprPrinter):
    """Preserve the undefined-function constructor metadata omitted by srepr.

    In particular, a real native immersion function must not revive as an
    assumption-free function. This changes encoding, not any constructed object
    or the strict structural D3 comparison.
    """
    def _print_FunctionClass(self, expr):
        if issubclass(expr, AppliedUndef) and expr._kwargs:
            return 'Function(%r, **%r)' % (expr.__name__, dict(sorted(expr._kwargs.items())))
        return super()._print_FunctionClass(expr)


def serialize(value):
    return LosslessReprPrinter().doprint(value)


def record(**items):
    return sp.Tuple(*(sp.Tuple(Str(k), v) for k, v in items.items()))


def text_tuple(*items):
    return sp.Tuple(*(Str(s) for s in items))


def action(role, *operands):
    # A standard SymPy undefined function makes the OPEN action serializable;
    # it supplies neither a constitutive equation nor a finite-state closure.
    return sp.Function('OPEN_' + role)(*operands)


def derivative(expr, coordinate):
    return sp.Derivative(expr, coordinate, evaluate=False)


def vector(items):
    return sp.ImmutableMatrix(list(items))


def bind_inputs(fold):
    # Deliberately empty, but exercised through the access-recording proxy.
    return {key: fold[key]['value'] for key in IMPORT_KEYS}


@lru_cache(maxsize=None)
def graph_kinematics(x, vr, xi):
    """Differentiate supplied graph/profile inputs only AFTER their selection."""
    r = sp.sqrt(sum(q*q for q in x))
    embedding = vector((*x, xi))
    tangent = embedding.jacobian(x)
    metric = sp.ImmutableMatrix(tangent.T * tangent)  # KNIFE_K4
    metric_inverse = sp.ImmutableMatrix(metric.inv(method='DM'))
    metric_determinant = sp.factor(metric.det())
    slope = vector(sp.diff(xi, q) for q in x)
    graph_normal = vector((*(-slope), sp.S.One)) / sp.sqrt(1 + slope.dot(slope))
    v = vector(vr*q/r for q in x)
    graph_velocity = tangent * v  # KNIFE_K2
    return (embedding, tangent, metric, metric_inverse, metric_determinant,
            graph_normal, v, graph_velocity)


def all_native_faces(face_set, face_label, contribution):
    """Formal sum over the ENTIRE OPEN O6 face set; Lambda binds its index.

    This is an unevaluated set-indexed sum, not one selected face, a chosen
    cardinality, or a finite-slab reduction. The same operator is used for work.
    """
    return action('SumOverAllNativeFaces', face_set, sp.Lambda(face_label, contribution))


def construct(_bound_inputs):
    # ACTION / ANSATZ inputs: supplied equations and named OPEN operands.
    x = sp.symbols('x1 x2 x3', real=True)
    t, w = sp.symbols('t w', real=True)
    ell, c0 = sp.symbols('ell c0', positive=True)
    GM = sp.Symbol('GM', real=True)
    r = sp.sqrt(sum(q*q for q in x))
    radial = {name: sp.Function(name)(r) for name in (
        'V_r', 'o2_rho_br_live', 'mu_perp', 'xi_w', 'h', 'delta', 'j_n', 'f')}
    vr, rho, mu, xi, h, delta, jn, f = (radial[k] for k in radial)
    domain = action('OpenDomain', *x, t)
    names = ('B_A13 I_br_live P_br_cons T_br_live T_br_cons N_br_live '
             'A_rot_live R_ref_strain_live M_perp R_br E_h_live J_map '
             'Pi_n H_core E_br_live J_E_live P_ref_relax_live '
             'P_convert_exchange_live P_boundary_live S_E_net P_E_supply '
             'C_ref S12_local_source_controller S12_boundary_domain '
             'S12_reaction_system S12_energy_reaction_supply '
             'mouth_core_data bulk_state nonpassive_reservoir_budget_obligation '
             'material_action_compatibility energy_accounting_overlap '
             'face_support_partition')
    op = dict(zip(names.split(), sp.symbols(names)))
    ref_stress = op['R_ref_strain_live']  # KNIFE_K5
    profile_packet = sp.Tuple(*radial.values())
    source = sp.Tuple(op['B_A13'], op['S12_local_source_controller'])
    boundary = sp.Tuple(op['S12_boundary_domain'], op['H_core'], op['mouth_core_data'])
    history = action('StateHistory', domain, profile_packet, source, boundary)
    state = sp.Tuple(domain, history, op['bulk_state'], profile_packet)
    # O1/O7 stay input obligations, not response equations or extra forces.
    constitutive_inputs = record(density=sp.Tuple(rho, op['R_br']),
                                 optical_stiffness=sp.Tuple(mu, op['M_perp']))
    # Compute graph geometry in ambient Cartesian (x1,x2,x3,w) coordinates.
    (embedding, tangent, metric, metric_inverse, metric_determinant,
     graph_normal, v, graph_velocity) = graph_kinematics(x, vr, xi)
    dual_tangent = metric_inverse * tangent.T
    geom = sp.Tuple(embedding, tangent, metric, metric_inverse, graph_normal)
    material = sp.Tuple(op['I_br_live'], op['P_br_cons'], op['N_br_live'],
                        op['A_rot_live'], op['J_map'], op['B_A13'],
                        op['material_action_compatibility'], constitutive_inputs)
    # An OPEN material response sees the FULL field/history domain. The explicit
    # profiles identify its supplied argument state, not a derivative cutoff.
    momentum_profiles = sp.Tuple(vr, rho, xi)  # KNIFE_K6
    momentum_replacements = dict(zip((vr, rho, xi), momentum_profiles))
    (m_embedding, m_tangent, m_metric, m_inverse, m_determinant,
     m_normal, m_v, momentum_velocity) = graph_kinematics(x, momentum_profiles[0], momentum_profiles[2])
    momentum_geometry = sp.Tuple(m_embedding, m_tangent, m_metric, m_inverse, m_normal)
    momentum_state = state.xreplace(momentum_replacements)
    momentum_material = material.xreplace(momentum_replacements)
    momentum_args = sp.Tuple(momentum_material, momentum_state, momentum_profiles,
                            momentum_velocity, momentum_geometry, source, boundary,
                            op['R_ref_strain_live'])
    density = vector(action('MomentumDensity_' + str(a), momentum_args) for a in range(4))
    flux = sp.ImmutableMatrix(4, 3, lambda a, i:
                             action('MomentumFlux_' + str(a) + '_' + str(i), momentum_args))
    storage = density.applyfunc(lambda p: derivative(p, t))
    transport = vector(sum(derivative(flux[a, i], x[i]) for i in range(3))
                       for a in range(4))  # KNIFE_K8
    # Explicit differentials of known input fields, alongside held derivatives
    # of the OPEN maps: no rho*V momentum map or advective flux law is inserted.
    momentum_input_differentials = record(
        profiles=momentum_profiles,
        gradients=sp.ImmutableMatrix(momentum_profiles).jacobian(x),
        velocity_gradient=momentum_velocity.jacobian(x),
        material_velocity_rate=momentum_velocity.jacobian(x) * m_v,
        measure=Str('d3x'))
    stress_state = sp.Tuple(op['T_br_live'], op['T_br_cons'], ref_stress,
                           op['N_br_live'], op['A_rot_live'], material, state,
                           geom, graph_velocity, source, boundary)
    internal = vector(action('InternalForce_' + str(a), stress_state) for a in range(4))
    compatibility = action('MaterialCompatibility', material, stress_state,
                           momentum_args, op['E_h_live'])
    # An arbitrary regular native immersion X_s(q,t) into (x1,x2,x3,w).
    # Chart orientation is chosen outward; no face is identified with the graph.
    # O6 owns its set, domains, geometry and reduction. No face count is chosen.
    s = sp.Symbol('o2_s_face')
    q = sp.symbols('o2_q_face_1 o2_q_face_2 o2_q_face_3', real=True)
    face_set = action('NativeFaceSet', op['J_map'], state, boundary)
    face_context = sp.Tuple(s, state, geom, op['N_br_live'], op['J_map'], boundary)
    native_embedding = vector(sp.Function('o2_X_face_' + str(a), real=True)(s, *q, t)
                              for a in range(4))
    native_tangent = native_embedding.jacobian(q)
    native_cofactor = vector((-1)**(a+3) * native_tangent.extract(
                              [b for b in range(4) if b != a], range(3)).det()
                             for a in range(4))
    native_area = sp.sqrt(native_cofactor.dot(native_cofactor))
    native_normal = native_cofactor / native_area  # KNIFE_K7
    native_domain = action('NativeFaceChartDomain', face_context)
    native_point_context = sp.Tuple(face_context, sp.Tuple(*q), native_embedding)
    amplitude = action('T_bulk_n_s_live', native_point_context)
    bulk_load = amplitude * native_normal
    hold = sp.Symbol('T_hold_s')
    # Single full load, with its restricted bulk part as an operand, never added
    # beside another full T_hold. External support has not been selected.
    native_hold = vector(action('FullFaceSupportLoad_' + str(a), hold,
                               sp.Tuple(*bulk_load), op['face_support_partition'],
                               native_point_context) for a in range(4))
    per_face_load = vector(action('FaceLoadReduction_' + str(a), op['J_map'],
                                  native_domain, sp.Lambda(q, sp.Tuple(native_area, native_hold)),
                                  face_context) for a in range(4))
    reduced_load = vector(all_native_faces(face_set, s, per_face_load[a]) for a in range(4))
    carried_v = v  # KNIFE_K3
    carried_w = action('CarriedMomentumW', op['Pi_n'], jn, op['J_map'],
                       op['N_br_live'], graph_velocity, geom, state, source,
                       boundary, Str('outward_native_relative_mass_current'),
                       Str('premise_3_local_material_velocity_at_each_transfer'))  # KNIFE_K13
    carry = vector((*list(jn * carried_v), carried_w))
    partners = vector(action('OutwardSourcePartner_' + str(a), op['Pi_n'], source,
                             op['S12_reaction_system'], op['J_map'], state, geom,
                             boundary) for a in range(4))
    body_entries = ()  # KNIFE_K12
    # Balance assembly is computed from signed physical accounting entries.
    entries = ((1, storage), (1, transport), (1, carry), (1, partners),
               (-1, internal), (-1, reduced_load)) + body_entries
    balance = sum((sign*entry for sign, entry in entries), sp.zeros(4, 1))
    balance = sp.ImmutableMatrix(balance)
    mass_density = rho  # KNIFE_K1
    mass_current = mass_density * v
    mass_divergence = sp.factor(sum(sp.diff(mass_current[i], x[i]) for i in range(3)))
    mass_equation = sp.Eq(mass_divergence, -jn, evaluate=False)
    mass_residual = mass_equation.lhs - mass_equation.rhs
    # Force/power on actual application velocities. The native face velocity
    # identification remains OPEN, and the same face context/map/area is used.
    face_velocity = vector(action('FaceApplicationVelocity_' + str(a),
                                  graph_velocity, native_point_context) for a in range(4))
    paired_face_velocity = face_velocity  # KNIFE_K9
    native_face_power = native_hold.dot(paired_face_velocity)
    per_face_power = action('FaceWorkReduction', op['J_map'], native_domain,
                            sp.Lambda(q, sp.Tuple(native_area, native_face_power)), face_context)
    mechanical_power = all_native_faces(face_set, s, per_face_power)
    material_pairing = internal.dot(graph_velocity)
    rotational_power = action('RotationalGeneralizedWork', op['A_rot_live'],
                              material, stress_state, state)
    energy_transport_channels = (op['J_E_live'],)  # KNIFE_K10
    material_work = action('MaterialStressNormalWork', material_pairing,
                           rotational_power, stress_state, op['E_br_live'],
                           sp.Tuple(*energy_transport_channels))
    relaxation = op['P_ref_relax_live']  # KNIFE_K11
    energy_context = sp.Tuple(state, material, stress_state, compatibility,
                              op['C_ref'], op['R_ref_strain_live'])
    energy_density = action('MaterialEnergyDensity', op['E_br_live'], energy_context)
    energy_flux = vector(sum(action('MaterialEnergyFlux_' + str(i), channel,
                                    energy_context) for channel in energy_transport_channels)
                         for i in range(3))
    energy_storage = derivative(energy_density, t)
    energy_transport = sum(derivative(energy_flux[i], x[i]) for i in range(3))
    # This joint action names the unresolved overlaps rather than adding named
    # powers that might describe the same work. The operator is NOT a closure.
    net_power = action('JointPowerAccounting',
        record(material_work=material_work, face_work=mechanical_power,
               relaxation=relaxation, conversion=op['P_convert_exchange_live'],
               boundary=op['P_boundary_live'], supplier=op['S_E_net'],
               budget=op['P_E_supply'], reactions=op['S12_energy_reaction_supply'],
               reference=op['C_ref'], overlap=op['energy_accounting_overlap'],
               storage_operand=energy_density, transport_operand=sp.Tuple(*energy_flux),
               nonpassive_obligation=op['nonpassive_reservoir_budget_obligation']),
        energy_context, source, boundary)
    energy_balance = sum((energy_storage, energy_transport, -net_power))
    optical = sp.Tuple(sp.Eq(action('c_gamma', r)**2, mu/rho, evaluate=False),
                       sp.Eq(action('c_gamma', r), c0*(1+delta), evaluate=False))
    grades = text_tuple('o2_rho_br_live', 'mu_perp', 'I_br_live', 'T_br_live', 'N_br_live',
                        'R_ref_strain_live', 'A_rot_live', 'P_ref_relax_live',
                        'Pi_n', 'j_n', 'T_hold_s', 'H_core', 'mouth_core_data',
                        'embedding_coefficients_source_longitudinal_field')
    restrictions = text_tuple(
        'far_field_isolated_spherical_mass_at_rest_lab_time',
        'Eulerian_steady_profiles_not_material_constancy',
        'linear_optical_waves_leading_eikonal',
        'radial_profiles_only_no_constitutive_isotropy_parity_or_stress_symmetry',
        'isotropic_optical_speed_same_shear_regime_only',
        'LAB_HELD_c_gamma_no_material_reference_or_holder_law',
        'O2_untruncated_optical_monomial_box_not_a_material_order_contract',
        'mass_law_relative_O_epsilon_qualification_for_induced_measure_claims',
        'bulk_f_first_order_recorded_domain_no_relation_to_epsilon',
        'no_finite_slab_or_sharp_sheet_material_reduction_selected',
        'S9_no_dissipation_frequency_independent_moduli_revisited_by_premise_1',
        'S9_flow_and_isotropic_strain_freezes_lifted',
        'S11b_rest_bulk_uniform_quadratic_breathing_slice_historical_only',
        'S11c_first_shape_background_jet_frozen_current_uniform_fold_historical_only',
        'L3_postulated_parent_L4_constant_L5_frozen_L6_static_held_mouth_historical_only',
        'homogeneous_linear_kinetic_continuity_anchor_not_flowing_kinetic_law',
        'uniform_transverse_anchor_not_live_mu_perp',
        'static_wall_tension_not_slab_tension_or_live_normal_response',
        'slab_kinematics_shear_projection_not_O6_O7_O1_closures',
        'rest_acoustics_perturbation_face_work_held_support_not_live_bulk_load',
        'real_fraction_order_work_and_n5_C_ref_record_not_live_energy_completion',
        'reference_comparison_EXPLORATORY_PAUSED_no_comparator_law_adopted')
    out = {
        'BASIS': record(coordinates=sp.Tuple(*x, w), time=t,
                        ambient_basis=sp.ImmutableMatrix(sp.eye(4)),
                        graph_tangent=tangent, graph_dual=dual_tangent,
                        graph_normal=graph_normal),
        'MEASURES': record(density_measure=Str('d3x'),
                           graph_area_factor=sp.sqrt(metric_determinant),
                           native_face_area=native_area, map=op['J_map'],
                           face_context=face_context, face_set=face_set,
                           native_chart_domain=native_domain,
                           chart_coordinates=sp.Tuple(*q),
                           orientation=Str('outward_oriented_regular_native_chart')),
        'GEOMETRY': record(embedding=embedding, metric=metric,
                           inverse=metric_inverse, determinant=metric_determinant,
                           identity=sp.Eq(xi, ell*h, evaluate=False),
                           native_embedding=native_embedding,
                           native_tangent=native_tangent,
                           native_cofactor=native_cofactor,
                           native_area=native_area, native_normal=native_normal,
                           native_regular_domain=sp.Ne(native_area**2, 0, evaluate=False),
                           native_normal_norm=sp.cancel(native_normal.dot(native_normal)),
                           native_normal_tangent_pairing=(native_tangent.T * native_normal).applyfunc(sp.expand)),
        'PROFILES': profile_packet,
        'MATERIAL_VELOCITY': graph_velocity,
        'MATERIAL_INPUT_DIFFERENTIALS': momentum_input_differentials,
        'MOMENTUM_DENSITY': density, 'MOMENTUM_FLUX': flux,
        'MOMENTUM_STORAGE': storage, 'MOMENTUM_TRANSPORT': transport,
        'MATERIAL_COMPATIBILITY': compatibility,
        'INTERNAL_FORCE': internal, 'NATIVE_BULK_LOAD': bulk_load,
        'NATIVE_HOLD_LOAD': native_hold, 'MECHANICAL_LOAD': reduced_load,
        'CARRIED_MOMENTUM': carry, 'SOURCE_PARTNERS': partners,
        'DRIVE': record(body_entries=sp.Tuple(*(sp.Tuple(a, b) for a, b in body_entries)),
                        local_source=source, boundary=boundary,
                        occurrences=sp.Tuple(internal, reduced_load, carry, partners),
                        representation=Str('dynamical_order_conversion_entry_provenance')),
        'HOLD_INPLANE': balance[:3, :], 'HOLD_W': balance[3],
        'HOLD_GRAPH_NORMAL': graph_normal.dot(balance),
        'MASS_INPUT': record(current=mass_current, divergence=mass_divergence,
                             outward_loss=jn, supplied_equation=mass_equation,
                             measure=Str('d3x'), density_input=constitutive_inputs),
        'MASS_RESIDUAL': mass_residual,
        'MATERIAL_POWER_PAIRING': record(velocity=graph_velocity, force=internal,
                                         pair=material_pairing, rotational=rotational_power,
                                         work=material_work),
        'FACE_POWER_PAIRING': record(velocity=paired_face_velocity, traction=native_hold,
                                     native_power=native_face_power,
                                     reduced_power=mechanical_power, map=op['J_map'],
                                     native_area=native_area),
        'ENERGY_STORAGE': energy_storage, 'ENERGY_TRANSPORT': energy_transport,
        'ENERGY_POWER': net_power, 'ENERGY_STEADY': energy_balance,
        'COUPLED_INPUTS': record(operands=sp.Tuple(*op.values()), optical=optical,
            embedding=sp.Tuple(op['E_h_live'], sp.Eq(xi, ell*h, evaluate=False),
                               Str('normal_equation_identity_and_count_unsettled')),
            unknown_grades_and_derivative_scales=grades,
            live_density_name=sp.Tuple(Str('ρ_br(r)'), Str('o2_rho_br_live')),
            constitutive_inputs=constitutive_inputs,
            open_action_semantics=Str('displayed_operands_not_closed_arguments_full_fields_jets_history_nonlocality'),
            native_reduction=Str('unevaluated_J_map_action_no_sheet_slab_choice'),
            supplier_obligation=sp.Tuple(op['S_E_net'], op['P_E_supply']),
            ownership=text_tuple('S12_conversion_map_reaction_supply', 'S14a_projected_flux_bridge',
                                 'S16_response_GM_matching_inference', 'S21_provenance_integration',
                                 'S1.5_S8_Q1_Q2_S22_material_holder_completion',
                                 'reference_relaxation_owner_unassigned',
                                 'nonpassive_live_successor_owner_unassigned')),
        'MODEL_POINT': record(restrictions=restrictions,
            premises=text_tuple('1_material_reference_optical_elastic_steady_relaxing',
                                '2_conversion_drive', '3_carried_local_material_velocity',
                                '4_shear_free_scalar_bulk_normal_load'),
            premise_status=Str('adopted substrate input to a conditional model (2026-10-06)'),
            epsilon=GM/(c0**2*r),
            optical_monomial_indices=sp.Tuple(*(sp.Tuple(a,b,c) for a in range(2)
                                               for b in range(3) for c in range(2))),
            recorded_grades=record(delta=sp.Symbol('epsilon'),
                                   slope_squared=sp.Symbol('epsilon'),
                                   velocity_over_c0=sp.sqrt(sp.Symbol('epsilon')))),
    }
    traces = {
        'BASIS MEASURES GEOMETRY PROFILES MATERIAL_VELOCITY': '§1,§3.1; C§1,2,5,6',
        'MATERIAL_INPUT_DIFFERENTIALS MOMENTUM_DENSITY MOMENTUM_FLUX MOMENTUM_STORAGE MOMENTUM_TRANSPORT MATERIAL_COMPATIBILITY': '§1,§3.2,§4; C§1,2,3,4,6',
        'INTERNAL_FORCE': '§3.2,§4; C§3,4',
        'NATIVE_BULK_LOAD NATIVE_HOLD_LOAD MECHANICAL_LOAD': '§2 premise 4,§3.3,§4,§5; C§1,6,7',
        'CARRIED_MOMENTUM SOURCE_PARTNERS': '§2 premise 3,§3.3,§5; C§1,6,7,10',
        'DRIVE': '§2 premise 2,§3.3,§4,§5; C§1,7,10',
        'HOLD_INPLANE HOLD_W HOLD_GRAPH_NORMAL': '§1–5; C§1–7,10; storage+outward transport+outward exchange-applied force',
        'MASS_INPUT MASS_RESIDUAL': '§1,§3.1,§7; C§6,9; supplied coordinate-measure law',
        'MATERIAL_POWER_PAIRING FACE_POWER_PAIRING ENERGY_STORAGE ENERGY_TRANSPORT ENERGY_POWER ENERGY_STEADY': '§2 premise 1,§6,§8.6; C§1,3,4,8; storage+outward transport-net power',
        'COUPLED_INPUTS MODEL_POINT TRACE': '§1–3,§7–10; C§1–10',
    }
    out['TRACE'] = sp.Tuple(*(sp.Tuple(Str(name), Str(trace))
                              for group, trace in traces.items() for name in group.split()))
    if set(out) != set(SECTION9):
        raise ValueError('section 9 object inventory')
    return out


def emit(name, value, local=False):
    prefix = 'PY_LOCAL_O2_' if local else 'PY_O2_'
    print(prefix + name + ': ' + serialize(value), flush=True)


def equal_three(left, right):
    """Total, three-valued comparison without added physical assumptions."""
    if left == right:
        return sp.true
    if isinstance(left, sp.MatrixBase) and isinstance(right, sp.MatrixBase):
        if left.shape != right.shape:
            return sp.false
        return aggregate_comparison(zip(left, right))
    if isinstance(left, (tuple, list, sp.Tuple)) and isinstance(right, (tuple, list, sp.Tuple)):
        if len(left) != len(right):
            return sp.false
        return aggregate_comparison(zip(left, right))
    if isinstance(left, dict) and isinstance(right, dict):
        if left.keys() != right.keys():
            return sp.false
        return aggregate_comparison((left[k], right[k]) for k in left)
    if isinstance(left, sp.core.relational.Relational) and isinstance(right, sp.core.relational.Relational):
        if left.func == right.func:
            return aggregate_comparison(zip(left.args, right.args))
        return Str('UNKNOWN')
    if isinstance(left, (str, Str, bool)) and isinstance(right, type(left)):
        return sp.sympify(left == right)
    if isinstance(left, sp.Expr) and isinstance(right, sp.Expr):
        result = (left-right).equals(0)
        return Str('UNKNOWN') if result is None else sp.sympify(result)
    return Str('UNKNOWN')


def aggregate_comparison(pairs):
    outcomes = [equal_three(a, b) for a, b in pairs]
    if any(value is sp.false for value in outcomes):
        return sp.false
    return sp.true if all(value is sp.true for value in outcomes) else Str('UNKNOWN')


def publish(fold, objects, emitted):
    complete = set(emitted) == set(SECTION9)
    emit('F6', sp.sympify(complete), local=True)
    if not complete:
        raise ValueError('section 9 publication coverage')
    roots = {'B_hold_live_x' + str(i+1): objects['HOLD_INPLANE'][i] for i in range(3)}
    roots.update(B_hold_live_w=objects['HOLD_W'],
                 B_hold_live_graph_normal=objects['HOLD_GRAPH_NORMAL'],
                 B_E_steady=objects['ENERGY_STEADY'])
    delta, routes, write_keys = {}, [], []
    for key, value in roots.items():
        route, write_key, comparison = 'F9A_ABSENT', key, Str('NOT_COMPARED_ABSENT')
        extra = {}
        if key in fold:
            prior = fold[key]['value']
            comparison = equal_three(prior, value)
            route = 'F9B_EQUAL' if comparison is sp.true else 'F9C_NEW'
            write_key = key if comparison is sp.true else 'o2_' + key
            extra = {'f9_operands': (serialize(prior), serialize(value)),
                     'f9_comparison': serialize(comparison)}
        if write_key in delta or (write_key != key and write_key in fold):
            raise ValueError('F9 routed-key collision: ' + write_key)
        delta[write_key] = dict(value=value, display=sp.sstr(value),
            value_kind='COMPUTED_OBJECT', **{'class': 'DERIVED'}, step='O2',
            route=route, evidence_tag='PY_O2_' + ('ENERGY_STEADY' if key == 'B_E_steady'
                                               else 'HOLD_GRAPH_NORMAL' if key.endswith('graph_normal')
                                               else 'HOLD_W' if key.endswith('_w') else 'HOLD_INPLANE'),
            **extra)
        routes.append(sp.Tuple(Str(key), Str(write_key), Str(route), comparison))
        write_keys.append(write_key)
    emit('F9', sp.Tuple(*routes), local=True)
    emit('F3', sp.Tuple(*(sp.Tuple(Str(k), text_tuple(*row['f9_operands']))
                         for k, row in delta.items() if 'f9_operands' in row)), local=True)
    combined = dict(fold) | delta
    closure = check_consumer(combined, write_keys)['closure']
    own_closure = set(closure) & set(delta)
    assert_delta_is_minimal(delta, own_closure)
    digests = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
               for p in (HERE, SPEC, *FOLD_PATHS, ROOT/'scripts/ledger_fold.py')}
    lines = ['# Generated O2 own-rows delta; no live profile is bound from the historical fold.',
             'from types import MappingProxyType', 'import sympy as sp',
             'from sympy.core.symbol import Str',
             'from sympy.functions.elementary.piecewise import ExprCondPair',
             '_RELATIONALS = {',
             "'Equality': lambda a,b: sp.Eq(a,b,evaluate=False),",
             "'Unequality': lambda a,b: sp.Ne(a,b,evaluate=False),",
             "'StrictGreaterThan': lambda a,b: sp.Gt(a,b,evaluate=False),",
             "'StrictLessThan': lambda a,b: sp.Lt(a,b,evaluate=False),",
             "'GreaterThan': lambda a,b: sp.Ge(a,b,evaluate=False),",
             "'LessThan': lambda a,b: sp.Le(a,b,evaluate=False),", '}',
             'def _restore(source):',
             "    return eval(source, {'__builtins__': {}, **vars(sp), 'Str': Str, 'ExprCondPair': ExprCondPair, **_RELATIONALS})",
             'IMPORT_KEYS = ' + repr(IMPORT_KEYS),
             'BUILD_INPUT_DIGESTS = MappingProxyType(' + repr(digests) + ')',
             'EXPORT_ROOTS = ' + repr(tuple(write_keys)), '_LEDGER = {']
    for key, row in delta.items():
        fields = [repr(k) + ': ' + ('_restore(' + repr(serialize(v)) + ')'
                                   if k == 'value' else repr(v)) for k, v in row.items()]
        lines.append(repr(key) + ': {' + ', '.join(fields) + '},')
    lines += ['}', 'LEDGER = MappingProxyType({k: MappingProxyType(v) for k,v in _LEDGER.items()})',
              'del _LEDGER', '']
    code = '\n'.join(lines)
    namespace = {}
    exec(compile(code, str(ROOT/'scripts/O2_exports.py'), 'exec'), namespace)
    restored = namespace['LEDGER']
    roundtrip = sp.Tuple(*(sp.sympify(delta[k]['value'] == restored[k]['value']) for k in delta))
    emit('ROUNDTRIP', roundtrip, local=True)
    if not all(item is sp.true for item in roundtrip):
        raise ValueError('D3 serialization roundtrip')
    restored_closure = check_consumer(dict(fold) | dict(restored), write_keys)['closure']
    assert_delta_is_minimal(restored, set(restored_closure) & set(restored))
    emit('EXPORT_CLOSURE', text_tuple(*sorted(restored_closure)), local=True)
    emit('INPUT_DIGESTS', sp.Tuple(*(sp.Tuple(Str(k), Str(v)) for k,v in digests.items())), local=True)
    # The only file write; after coverage, F9, closure, minimality and D3 guards.
    (ROOT/'scripts/O2_exports.py').write_text(code)


def name_audit(objects, fold):
    """Census every emitted physical symbol/function name against upstream keys.

    Metadata-only inspection: no upstream value becomes a construction operand.
    The five shared coordinate names have the same meanings in the fold.
    """
    names = set()
    for value in objects.values():
        names.update(atom.name for atom in value.atoms(sp.Symbol))
        names.update(atom.func.__name__ for atom in value.atoms(AppliedUndef))
    overlap = sorted(names.intersection(fold))
    shared_coordinates = {'x1', 'x2', 'x3', 't', 'w'}
    unexpected = sorted(set(overlap) - shared_coordinates)
    emit('NAME_AUDIT', record(names=text_tuple(*sorted(names)),
         upstream_matches=sp.Tuple(*(sp.Tuple(Str(k), Str(fold[k].get('class', '')),
                                  Str(fold[k].get('description', ''))) for k in overlap)),
         other_matches=text_tuple(*unexpected)), local=True)
    if unexpected:
        raise ValueError('unresolved upstream name collision: ' + ', '.join(unexpected))


def run(publish_delta=True, fold=None):
    if fold is None:
        fold, audit = load_model(*(str(p) for p in FOLD_PATHS))
    else:
        audit = {'source_row_counts': [(str(p), Str('shared_fold')) for p in FOLD_PATHS]}
    closure = check_consumer(fold, IMPORT_KEYS)
    witness = assert_lookups_equal_manifest(bind_inputs, fold, IMPORT_KEYS)
    objects = construct(witness['result'])
    emit('LOCAL_NAMES', text_tuple(*('PY_LOCAL_O2_' + n for n in LOCAL_NAMES)))
    emit('IMPORT_MANIFEST', text_tuple(*sorted(witness['lookups'])), local=True)
    emit('FOLD', record(rows=sp.Integer(len(fold)),
         import_closure=text_tuple(*sorted(closure['closure'])),
         files=text_tuple(*(str(p.relative_to(ROOT)) for p in FOLD_PATHS))), local=True)
    emitted = []
    for name in SECTION9:
        emit(name, objects[name])
        emitted.append(name)
    name_audit(objects, fold)
    if publish_delta:
        publish(fold, objects, emitted)
    return objects


if __name__ == '__main__':
    run()
