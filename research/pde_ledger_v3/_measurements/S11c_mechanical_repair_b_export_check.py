#!/usr/bin/env python3
"""Four-case face-work, action-orientation, and baseline preservation operands."""
import hashlib
import json
from pathlib import Path
import resource
import sys
import time

STARTED = time.monotonic()
ROOT = Path(__file__).resolve().parents[1]
BASE = Path('/tmp/s11c-mechanical-repair-20260912')
sys.path.insert(0, str(ROOT / 'scripts'))
import sympy as sp
from ledger_fold import _restore
from S11c_inertia_artifact_audit import export_data, locate
from S11c_mechanical_sign_audit import Output, a, b, d, delta
import S11b_interface_coupling_law_sympy_audit as legacy


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mechanical(body, slot='EXPANDED'):
    return sp.ImmutableMatrix([*b.named_tuple_row(b.named_tuple_row(body, 'U_BODY_BALANCE'), slot),
                               b.named_tuple_row(b.named_tuple_row(body, 'E_W_BALANCE'), slot)])


def case_rows(root):
    return {tuple(map(str, axes)): b.named_tuple_row(payload, 'VALUE') for axes, payload in root}


def run():
    run = BASE / 'b_checks'
    run.mkdir(exist_ok=True)
    sources = [Path(__file__), ROOT / 'scripts/S11c_b_exports.py',
               ROOT / 'scripts/S11c_b_brane_operator_sympy_audit.py',
               ROOT / 'scripts/S11c_a_exports.py', ROOT / 'scripts/S11c_a_interface_geometry_sympy_audit.py',
               ROOT / 'scripts/S11b_interface_coupling_law_sympy_audit.py',
               ROOT / 'scripts/S11c_d_mixing_scattering_sympy_audit.py',
               ROOT / '_measurements/S11c_mechanical_sign_audit.py',
               BASE / 'baseline/scripts/S11c_b_exports.py']
    pins = {str(path): digest(path) for path in sources}
    new_values, export_pins, _ = export_data(ROOT / 'scripts/S11c_b_exports.py')
    old_values, _, _ = export_data(BASE / 'baseline/scripts/S11c_b_exports.py')
    new, old = case_rows(_restore(new_values['slab_operator'])), case_rows(_restore(old_values['slab_operator']))
    origins = case_rows(_restore(new_values['slab_operator_term_origins']))
    old_origins = case_rows(_restore(old_values['slab_operator_term_origins']))
    basis = {str(branch): b.named_tuple_row(payload, 'VALUE') for branch, payload in _restore(new_values['energy_basis_variable'])}
    for energy in basis.values():
        for row in energy[1:]:
            if str(row[0]) in ('W_BG', 'MU_R_BG'):
                b.SYMBOL_DIMENSIONS[row[3]] = row[5] - b.dimension_of(row[2])
    out = Output()
    out.emit('s11cMechanicalRepairFullSourcePins', pins)
    out.emit('s11cMechanicalRepairFullExportPinResidual',
             {name: int(digest(locate(name)) != value) for name, value in export_pins.items()}, residual=True)
    changed = sorted(key for key in new_values.keys() & old_values.keys() if new_values[key] != old_values[key])
    out.emit('s11cMechanicalRepairFullChangedValueSerializations', changed)
    old_c1, _, c1_imports = export_data(BASE / 'baseline/scripts/S11c_c1_exports.py')
    out.emit('s11cMechanicalRepairFullChangedC1DirectInputs', sorted(set(changed) & set(c1_imports)))
    legacy_model = legacy.derive_model()
    legacy_coordinates = (legacy.delta_p_plus, legacy.A_plus_affinity, legacy.delta_p_minus, legacy.A_minus_affinity)
    legacy_internal = legacy_model['thickness_uneliminated'].subs(dict.fromkeys(legacy_coordinates, sp.S.Zero))
    legacy_work = sp.expand(legacy_model['thickness_uneliminated']-legacy_internal)
    legacy_live_bindings = {symbol: b.INCOMING_LEDGER[symbol.name]['value']
                            for symbol in legacy_work.free_symbols
                            if symbol.name in b.INCOMING_LEDGER and symbol not in legacy_coordinates}
    out.emit('s11cMechanicalRepairFullS11bLiveCoordinateBindings',
             tuple((sp.core.symbol.Str(sp.srepr(source)), target)
                   for source, target in sorted(legacy_live_bindings.items(), key=lambda item: sp.default_sort_key(item[0]))))
    vector_units = lambda path: tuple(b.DIM_BODY_U if path[0] < 3 else b.DIM_ENERGY)
    for case, body in new.items():
        branch, representative = case
        prefix = 's11cMechanicalRepairFull' + ''.join(word.title().replace('_', '') for word in case)
        def emit(name, value, **kwargs):
            out.emit(prefix + name, value, **kwargs)
        before, after = mechanical(old[case]), mechanical(body)
        bundle = b.named_tuple_row(body, 'FACE_FLUX_BOUNDARY_OPERANDS')
        evolution = b.named_tuple_row(b.named_tuple_row(body, 'THETA_BALANCE'), 'EVOLUTION_TERM_ORIGINS')
        mu = b.named_tuple_row(body, 'MU_THETA_FACE_BINDING')[1]
        physical = b.retained_grade(b.face_generalized_force_rows(bundle, branch, representative, evolution, mu))
        force = sp.ImmutableMatrix([*physical['U'], physical['E_W']])
        recorded = b.named_tuple_row(body, 'FACE_GENERALIZED_FORCE_ROWS')
        recorded_force = sp.ImmutableMatrix([*b.named_tuple_row(recorded, 'U'), b.named_tuple_row(recorded, 'E_W')])
        controls = {b.INCOMING_LEDGER[name]['value']: sp.S.Zero for name in
                    ('delta_p_plus', 'delta_p_minus', 'd_w_delta_p_plus', 'd_w_delta_p_minus', 'Lambda_X_0')}
        internal = after.xreplace(controls)
        increment = (after - internal).applyfunc(sp.expand)
        energy = sp.Add(*(row[4] if str(row[0]) in ('W_BG', 'MU_R_BG') else row[3] for row in basis[branch][1:]))
        action_anchor = sp.expand(-sp.diff(energy, b.k_W, b.e_W, 2) / b.epsilon).coeff(b.epsilon, 1)
        row_anchor = sp.expand(sp.diff(internal[3], b.k_W, b.e_W)).coeff(b.epsilon, 1)
        multiplier = sp.cancel(row_anchor / action_anchor)
        emit('StiffnessOrientationOperands', (action_anchor, row_anchor, multiplier))
        reciprocal_bases = {power.base for expression in (*after, *force, action_anchor)
                            for power in expression.atoms(sp.Pow) if power.exp.is_negative}
        emit('ReciprocalDomainOperands', tuple(sorted(reciprocal_bases, key=sp.default_sort_key)), heavy=True)
        emit('NormalizationDivisor', action_anchor)
        emit('PhysicalForce', force, units=vector_units, heavy=True)
        emit('ForceExtractionResidual', delta(force, recorded_force), units=vector_units, residual=True)
        emit('BeforeRows', before, units=vector_units, heavy=True)
        emit('AfterRows', after, units=vector_units, heavy=True)
        emit('AssembledLoad', increment, units=vector_units, heavy=True)
        emit('ActionNormalizedLoad', multiplier*force, units=vector_units, heavy=True)
        emit('ActionLoadResidual', delta(increment, multiplier * force), units=vector_units, residual=True)
        emit('NonFacePreservationResidual', delta(internal, before.xreplace(controls)), units=vector_units, residual=True)
        old_recorded = b.named_tuple_row(old[case], 'FACE_GENERALIZED_FORCE_ROWS')
        old_force = sp.ImmutableMatrix([*b.named_tuple_row(old_recorded, 'U'), b.named_tuple_row(old_recorded, 'E_W')])
        emit('PhysicalForcePreservationResidual', delta(recorded_force, old_force), units=vector_units, residual=True)
        old_mass = b.named_tuple_row(b.named_tuple_row(old[case], 'THETA_BALANCE'), 'EXPANDED')
        mass = b.named_tuple_row(b.named_tuple_row(body, 'THETA_BALANCE'), 'EXPANDED')
        emit('MassPreservationResidual', delta(mass, old_mass), units=tuple(b.DIM_RHOBR-b.DIM_T), residual=True)
        emit('ChemicalPreservationResidual', delta(mu, b.named_tuple_row(old[case], 'MU_THETA_FACE_BINDING')[1]),
             units=tuple(b.DIM_ENERGY), residual=True)
        for label in ('THETA_BALANCE', 'MU_THETA_FACE_BINDING'):
            current, previous = b.named_tuple_row(body, label), b.named_tuple_row(old[case], label)
            emit(label.title().replace('_', '') + 'ProvenanceIdentity',
                 (digest_text(current), digest_text(previous), int(current != previous)))
        for label in ('KINETIC', 'BULK_ENERGY'):
            current = b.named_tuple_row(origins[case], label)
            previous = b.named_tuple_row(old_origins[case], label)
            # Hashes record preservation of the entire structured provenance;
            # physical zero rows below retain their operand-derived dimensions.
            emit(label.title().replace('_', '') + 'ProvenanceIdentity',
                 (digest_text(current), digest_text(previous), int(current != previous)))
        time_coordinate = sp.Symbol('s11cMechanicalRepairTime', real=True)
        fields = tuple(sp.Function('s11cMechanicalRepairField'+str(i))(time_coordinate) for i in range(4))
        rho = b.density_pair(representative)[1]
        kinetic = b.epsilon**2*(rho*sum(sp.diff(field, time_coordinate)**2 for field in fields[:3])
                               + b.mu_W*sp.diff(b.W_bg*fields[3], time_coordinate)**2)/2
        accelerations = (*b.u_tt, b.e_tt)
        mapping = {sp.diff(field, time_coordinate, 2): acceleration for field, acceleration in zip(fields, accelerations)}
        action_kinetic = b.retained_grade(sp.ImmutableMatrix([
            -sp.diff(sp.diff(kinetic, sp.diff(field, time_coordinate)), time_coordinate).xreplace(mapping)/b.epsilon
            for field in fields]))
        incoming_kinetic = b.named_tuple_row(origins[case], 'KINETIC')
        incoming_kinetic = sp.ImmutableMatrix([*incoming_kinetic[0], incoming_kinetic[1]])
        inertia_anchor = sp.expand(sp.diff(action_kinetic[3], b.mu_W, b.e_tt)).coeff(b.epsilon, 1)
        inertia_row = sp.expand(sp.diff(internal[3], b.mu_W, b.e_tt)).coeff(b.epsilon, 1)
        emit('KineticActionRows', action_kinetic, units=vector_units)
        emit('IncomingKineticRows', incoming_kinetic, units=vector_units)
        emit('NormalizedKineticRows', multiplier*action_kinetic, units=vector_units)
        emit('InertiaOrientationOperands', (inertia_anchor, inertia_row, sp.cancel(inertia_row/inertia_anchor)))
        emit('InertiaStiffnessOrientationResidual', sp.cancel(inertia_row/inertia_anchor-multiplier), residual=True)
        emit('KineticNormalizationResidual', delta(incoming_kinetic, multiplier*action_kinetic), units=vector_units, residual=True)
        source_origin = b.named_tuple_row(b.named_tuple_row(origins[case], 'FACE_VIRTUAL_WORK'), 'ROWS')
        origin_vector = sp.ImmutableMatrix([*b.named_tuple_row(b.named_tuple_row(source_origin, 'U'), 'EXPANDED'),
                                           b.named_tuple_row(b.named_tuple_row(source_origin, 'E_W'), 'EXPANDED')])
        emit('OriginLoad', origin_vector, units=vector_units, heavy=True)
        emit('OriginActionLoadResidual', delta(origin_vector, multiplier * force), units=vector_units, residual=True)
        for slot in ('LOCAL', 'EXPANDED'):
            emit(slot.title() + 'LoadDeltaResidual',
                 delta(mechanical(body, slot)-mechanical(old[case], slot), increment-(before-before.xreplace(controls))),
                 units=vector_units, residual=True)
        impermeable = {b.INCOMING_LEDGER['Lambda_A_0']['value']: sp.S.Zero,
                       b.INCOMING_LEDGER['Lambda_V_0']['value']: sp.S.Zero}
        emit('ImpermeableInputControl', tuple(sp.Eq(k, v, evaluate=False) for k, v in impermeable.items()))
        emit('ImpermeableActionLoad', (multiplier * force).xreplace(impermeable), units=vector_units, heavy=True)
        emit('ImpermeableAssembledLoad', increment.xreplace(impermeable), units=vector_units, heavy=True)
        emit('ImpermeableActionLoadResidual', delta(increment.xreplace(impermeable), (multiplier * force).xreplace(impermeable)),
             units=vector_units, residual=True)
        emit('ZeroSourceControl', tuple(sp.Eq(k, v, evaluate=False) for k, v in controls.items()))
        emit('ZeroSourcePhysicalForce', force.xreplace(controls).applyfunc(sp.expand), units=vector_units)
        work = b.named_tuple_row(physical['SOURCE_OPERANDS'], 'VIRTUAL_WORK_SHAPE_DERIV')
        face_sum = sp.zeros(4, 1)
        for face_label, density in work[0]:
            face = {'UPPER': 1, 'LOWER': -1}[str(face_label)]
            face_name = 'Plus' if face == 1 else 'Minus'
            face_force = sp.ImmutableMatrix([sp.diff(density, test) for test in (*b.delta_v_u, b.delta_v_e_W)])
            face_sum += face_force
            emit(face_name + 'VirtualWork', density, heavy=True)
            emit(face_name + 'GeneralizedForce', face_force, units=vector_units, heavy=True)
            emit(face_name + 'WorkReconstructionResidual',
                 sp.cancel(density-(face_force.T*sp.ImmutableMatrix([*b.delta_v_u, b.delta_v_e_W]))[0]),
                 units=tuple(b.DIM_ENERGY), residual=True)
            source = a.build_face_source(branch, int(face), 'DELTA_W', representative)
            pressure = b.uniformize(a.finalize(a.pressure_trace_raw(source)) / b.epsilon)
            affinity = b.uniformize(b.bind_mu_theta_operand(a.finalize(a.affinity_raw(source)) / b.epsilon, branch, mu))
            own = legacy_coordinates[:2] if face == 1 else legacy_coordinates[2:]
            other = legacy_coordinates[2:] if face == 1 else legacy_coordinates[:2]
            legacy_face = legacy_work.subs(dict.fromkeys(other, sp.S.Zero))
            legacy_face = legacy_face.xreplace(legacy_live_bindings)
            legacy_bound = b.uniformize(legacy_face.subs(dict(zip(own, (pressure, affinity))), simultaneous=True))
            uniform_load = b.uniformize(multiplier * face_force[3] / b.epsilon)
            emit(face_name + 'UniformS11bBindingOperands', (pressure, affinity), heavy=True)
            emit(face_name + 'UniformLoad', uniform_load, heavy=True)
            emit(face_name + 'UniformS11bLoad', legacy_bound, heavy=True)
            emit(face_name + 'UniformS11bLoadResidual', sp.cancel(uniform_load-legacy_bound),
                 units=tuple(b.DIM_ENERGY), residual=True)
        emit('FaceSumForceResidual', delta(force, face_sum), units=vector_units, residual=True)
    out.emit('s11cMechanicalRepairFullSourceStabilityResidual',
             {name: int(digest(Path(name)) != value) for name, value in pins.items()}, residual=True)
    report = {'sourcePins': pins, 'records': out.records,
              'wallSeconds': time.monotonic()-STARTED, 'peakRssKiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              'changedValueSerializations': changed, 'addedKeys': sorted(new_values.keys()-old_values.keys()),
              'removedKeys': sorted(old_values.keys()-new_values.keys())}
    (run / 'inventory.json').write_text(json.dumps(report, indent=2) + '\n')


def digest_text(value):
    return hashlib.sha256(sp.srepr(value).encode()).hexdigest()


if __name__ == '__main__':
    run()
