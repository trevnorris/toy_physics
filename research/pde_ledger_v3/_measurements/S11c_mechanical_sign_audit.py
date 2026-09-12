#!/usr/bin/env python3
"""Pinned native virtual-work extraction and mechanical-row sign operands.

Print computed objects and literal residuals. Interpretation belongs to the
companion report. This instrument neither repairs nor regenerates an operator.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import resource
import sys
import time

STARTED = time.monotonic()
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import sympy as sp
from sympy.core.symbol import Str
from ledger_fold import _restore
from S11c_inertia_artifact_audit import export_data, locate
import S11c_a_interface_geometry_sympy_audit as a
import S11c_b_brane_operator_sympy_audit as b
import S11c_d_mixing_scattering_sympy_audit as d

SOURCE_PATHS = (
    '_measurements/S11c_mechanical_sign_audit.py',
    '_measurements/S11c_inertia_artifact_audit.py',
    'scripts/ledger_fold.py', 'scripts/S11c_d_output_codec.py',
    'scripts/S11b_interface_coupling_law_sympy_audit.py',
    'scripts/S11c_a_interface_geometry_sympy_audit.py',
    'scripts/S11c_b_brane_operator_sympy_audit.py',
    'scripts/S11c_c1_bulk_closure_sympy_audit.py',
    'scripts/S11c_c2_selfenergy_fold_sympy_audit.py',
    'scripts/S11c_d_mixing_scattering_sympy_audit.py',
    *(f'scripts/S11{stage}_exports.py' for stage in ('b', 'c_a', 'c_b', 'c_c1', 'c_c2')),
    *(f'directives/S11{stage}_SHARED_PHYSICS.md' for stage in ('b', 'c_a', 'c_b', 'c_c1', 'c_c2', 'c_d')),
)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def delta(left, right):
    return b.map_object(d.tree_difference(b.casify(left), b.casify(right)), sp.cancel)


class Output:
    def __init__(self):
        self.records = []
        self.keys = set()

    def emit(self, key, value, *, units=None, residual=False, heavy=False):
        if key in self.keys or key in b.INCOMING_LEDGER:
            raise ValueError(('write-key collision', key))
        self.keys.add(key)
        value = d.cas(value)
        metadata = []
        scalar_values = []
        for path, leaf in d.leaves(value):
            if isinstance(leaf, Str):
                continue
            measured = leaf.lhs if isinstance(leaf, sp.Equality) and leaf.rhs == 0 else leaf
            dim = tuple(b.dimension_of(measured)) if leaf != 0 and not (leaf.is_number and units is not None) else (
                units(path) if callable(units) else units if units is not None else tuple(b.DIM_ZERO))
            metadata.append({'path':list(path), 'dimensionLTM':list(map(str, dim)),
                             'epsEtaSigmaOrder':str(b.multigrade(leaf))})
            scalar_values.append(leaf)
        representation = 'carrierSha256AndNumericPit' if heavy else 'sympySrepr'
        payload = sp.srepr(d.carrier_fingerprint(value) if heavy else value)
        record = {'writeKey':key, 'representation':representation, 'object':payload,
                  'metadata':metadata}
        print(json.dumps(record, separators=(',', ':')), flush=True)
        self.records.append({'writeKey':key, 'scalarLeaves':len(scalar_values),
            'payloadBytes':len(payload.encode()), 'residual':residual,
            'zeroLeaves':sum(leaf == 0 for leaf in scalar_values),
            'nonzeroLeaves':sum(leaf != 0 for leaf in scalar_values),
            'nonfiniteLeaves':sum(leaf.has(sp.nan, sp.zoo, sp.oo, -sp.oo) for leaf in scalar_values),
            'unknownDimensions':sum('nan' in entry['dimensionLTM'] for entry in metadata)})


def native_case(out, axes, payload, basis):
    branch, representative = map(str, axes)
    suffix = ''.join(word.title().replace('_', '') for word in (branch, representative))
    def emit(name, value, **kwargs):
        out.emit('s11cMechanicalSign'+suffix+name, value, **kwargs)
    body = b.named_tuple_row(payload, 'VALUE')
    bundle = b.named_tuple_row(body, 'FACE_FLUX_BOUNDARY_OPERANDS')
    evolution = b.named_tuple_row(b.named_tuple_row(body, 'THETA_BALANCE'), 'EVOLUTION_TERM_ORIGINS')
    mu = b.named_tuple_row(body, 'MU_THETA_FACE_BINDING')[1]
    source = b.selected_substrate_axes(bundle, 'virtual_work_shape_deriv',
                                     (branch, 'DELTA_W', 'DELTA_W', representative))
    replay_a = a.finalize(a.virtual_work_cases(branch, 'DELTA_W', 'DELTA_W', representative))
    emit('NativeGeometryWork', replay_a, heavy=True)
    emit('ConsumedGeometryWork', source, heavy=True)
    emit('GeometryReplayResidual', delta(replay_a, source), units=tuple(b.DIM_ENERGY), residual=True)
    replay = b.face_generalized_force_rows(bundle, branch, representative, evolution, mu)
    replay = b.retained_grade(replay)
    stored = b.named_tuple_row(body, 'FACE_GENERALIZED_FORCE_ROWS')
    q = sp.ImmutableMatrix([*replay['U'], replay['E_W']])
    stored_q = sp.ImmutableMatrix([*b.named_tuple_row(stored, 'U'), b.named_tuple_row(stored, 'E_W')])
    row_units = [tuple(b.DIM_ENERGY-b.DIM_L)]*3+[tuple(b.DIM_ENERGY)]
    vector_units = lambda path:row_units[path[0]]
    emit('ExtractedExternalForce', q, units=vector_units, heavy=True)
    emit('StoredExternalForce', stored_q, units=vector_units, heavy=True)
    emit('ForceReplayResidual', delta(q, stored_q), units=vector_units, residual=True)
    work = b.named_tuple_row(replay['SOURCE_OPERANDS'], 'VIRTUAL_WORK_SHAPE_DERIV')
    tests = sp.ImmutableMatrix([*b.delta_v_u, b.delta_v_e_W])
    emit('BoundFaceWork', work, heavy=True)
    emit('WorkReconstructionResidual', sp.cancel(work[1]-(q.T*tests)[0]),
         units=tuple(b.DIM_ENERGY), residual=True)
    emit('FaceSumResidual', sp.cancel(work[1]-sum(row[1] for row in work[0])),
         units=tuple(b.DIM_ENERGY), residual=True)

    mechanical = sp.ImmutableMatrix([
        *b.named_tuple_row(b.named_tuple_row(body, 'U_BODY_BALANCE'), 'EXPANDED'),
        b.named_tuple_row(b.named_tuple_row(body, 'E_W_BALANCE'), 'EXPANDED')])
    control_names = ('delta_p_plus', 'delta_p_minus', 'd_w_delta_p_plus', 'd_w_delta_p_minus', 'Lambda_X_0')
    controls = {b.INCOMING_LEDGER[name]['value']:sp.S.Zero for name in control_names}
    emit('ZeroFaceSourceControl', tuple(sp.Eq(k,v,evaluate=False) for k,v in controls.items()))
    controlled = mechanical.xreplace(controls)
    face_increment = (mechanical-controlled).applyfunc(sp.expand)
    emit('AssembledMechanicalRows', mechanical, units=vector_units, heavy=True)
    emit('ControlledMechanicalRows', controlled, units=vector_units, heavy=True)
    emit('AssembledFaceIncrement', face_increment, units=vector_units, heavy=True)
    emit('ControlledExternalForce', q.xreplace(controls).applyfunc(sp.expand), units=vector_units)
    emit('AssemblyMinusExternalForceResidual', (face_increment-q).applyfunc(sp.cancel),
         units=vector_units, residual=True)

    # The action is T-U. Its Euler convention is derived before applying the
    # same normalization to the independently supplied external virtual work.
    energy = sp.Add(*(row[4] if str(row[0]) in ('W_BG', 'MU_R_BG') else row[3] for row in basis[1:]))
    stiffness_action = sp.expand(-sp.diff(energy, b.k_W, b.e_W, 2)/b.epsilon).coeff(b.epsilon, 1)
    stiffness_row = sp.expand(sp.diff(controlled[3], b.k_W, b.e_W)).coeff(b.epsilon, 1)
    stiffness_scale = sp.cancel(stiffness_row/stiffness_action)
    rho = b.density_pair(representative)[1]
    velocities = (*b.u_t, b.e_t)
    accelerations = (*b.u_tt, b.e_tt)
    time_coordinate = sp.Symbol('s11cMechanicalSignTime', real=True)
    fields = tuple(sp.Function('s11cMechanicalSignField'+str(i))(time_coordinate) for i in range(4))
    time_jets = {sp.diff(field,time_coordinate,order):jet
                 for order,jets in ((1,velocities),(2,accelerations)) for field,jet in zip(fields,jets)}
    kinetic_time = b.epsilon**2*(rho*sum(sp.diff(f,time_coordinate)**2 for f in fields[:3])+
                    b.mu_W*sp.diff(b.W_bg*fields[3],time_coordinate)**2)/2
    kinetic = kinetic_time.xreplace(time_jets)
    kinetic_euler = b.retained_grade(sp.ImmutableMatrix([
        -sp.diff(sp.diff(kinetic_time,sp.diff(field,time_coordinate)),time_coordinate).xreplace(time_jets)/b.epsilon
        for field in fields]))
    native_kinetic_u, native_kinetic_e = b.kinetic_balance_from_energy(rho)
    native_kinetic = b.retained_grade(sp.ImmutableMatrix([*native_kinetic_u, native_kinetic_e]))
    inertia_row = sp.expand(sp.diff(controlled[3], b.mu_W, b.e_tt)).coeff(b.epsilon, 1)
    inertia_action = sp.expand(sp.diff(kinetic_euler[3], b.mu_W, b.e_tt)).coeff(b.epsilon, 1)
    inertia_scale = sp.cancel(inertia_row/inertia_action)
    emit('StoredEnergy', energy, heavy=True)
    emit('KineticEnergy', b.retained_grade(kinetic))
    emit('StiffnessActionCoefficient', sp.expand(stiffness_action))
    emit('StiffnessRowCoefficient', sp.expand(stiffness_row))
    emit('ActionToRowScale', stiffness_scale)
    emit('KineticActionRows', kinetic_euler, units=vector_units)
    emit('NativeKineticRows', native_kinetic, units=vector_units)
    emit('InertiaActionCoefficient', inertia_action)
    emit('InertiaRowCoefficient', inertia_row)
    emit('KineticActionToRowScale', inertia_scale)
    emit('OrientationAnchorResidual', sp.cancel(stiffness_scale-inertia_scale), residual=True)
    emit('KineticNormalizationResidual', (native_kinetic-stiffness_scale*kinetic_euler).applyfunc(sp.cancel),
         units=vector_units, residual=True)
    lhs_load = stiffness_scale*q
    emit('ActionNormalizedFaceLoad', lhs_load, units=vector_units, heavy=True)
    emit('AssemblyMinusActionLoadResidual', (face_increment-lhs_load).applyfunc(sp.expand),
         units=vector_units, residual=True, heavy=True)
    emit('AssemblyPlusActionLoadResidual', (face_increment+lhs_load).applyfunc(sp.cancel),
         units=vector_units, residual=True)
    denominator_bases = {power.base for expression in (*mechanical,*q,stiffness_action,inertia_action)
                         for power in expression.atoms(sp.Pow) if power.exp.is_negative}
    emit('OriginalReciprocalBases', tuple(sorted(denominator_bases,key=sp.default_sort_key)), heavy=True)
    emit('NormalizationDivisors', (stiffness_action,inertia_action))
    emit('LiveFaceParameterSymbols', tuple(sorted(q.free_symbols, key=sp.default_sort_key)))


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', choices=('LAB_HELD_RHO4_CONSTANT', 'ALL'), default='ALL')
    parser.add_argument('--inventory', type=Path, required=True)
    args = parser.parse_args()
    source_pins = {name:sha(ROOT/name) for name in SOURCE_PATHS}
    out = Output()
    out.emit('s11cMechanicalSignSourcePins', source_pins)
    values, pins, _ = export_data(ROOT/'scripts/S11c_b_exports.py')
    slab = _restore(values['slab_operator'])
    basis = {str(branch):b.named_tuple_row(payload, 'VALUE')
             for branch,payload in _restore(values['energy_basis_variable'])}
    # New-invariant coefficients carry their dimensions in the consumed basis.
    for body in basis.values():
        for row in body[1:]:
            if str(row[0]) in ('W_BG', 'MU_R_BG'):
                b.SYMBOL_DIMENSIONS[row[3]] = row[5]-b.dimension_of(row[2])
    pin_residuals = {name:int(sha(locate(name)) != expected) for name,expected in pins.items()}
    out.emit('s11cMechanicalSignConsumedExportPinResidual', pin_residuals, residual=True)
    for axes,payload in slab:
        if args.case != 'ALL' and '_'.join(map(str, axes)) != args.case:
            continue
        native_case(out, axes, payload, basis[str(axes[0])])
    # Read the actual inherited S11b closed thickness row and its separately
    # emitted bulk contribution; no homogeneous model is re-entered by hand.
    thickness = b.INCOMING_LEDGER['thickness_eom']['value'].lhs
    bulk = b.INCOMING_LEDGER['bulk_force_on_thickness']['value']
    dimensional_aliases = {}
    for symbol in thickness.free_symbols | bulk.free_symbols:
        if symbol in b.SYMBOL_DIMENSIONS:
            continue
        if symbol.name in ('omega', 'k', 'q_out', 'eta'):
            dimension = {'omega':-b.DIM_T,'k':-b.DIM_L,'q_out':-b.DIM_L,'eta':b.DIM_ZERO}[symbol.name]
        else:
            candidates = {tuple(dim) for atom,dim in b.SYMBOL_DIMENSIONS.items()
                          if isinstance(atom,sp.Symbol) and atom.name == symbol.name}
            if len(candidates) != 1:
                raise ValueError(('legacy symbol dimension', sp.srepr(symbol), candidates))
            dimension = sp.ImmutableMatrix(next(iter(candidates)))
        b.SYMBOL_DIMENSIONS[symbol] = dimension
        dimensional_aliases[symbol.name] = (symbol, dimension)
    out.emit('s11cMechanicalSignS11bDimensionBindings', dimensional_aliases)
    out.emit('s11cMechanicalSignS11bThicknessRow', thickness, heavy=True)
    out.emit('s11cMechanicalSignS11bBulkContribution', bulk, heavy=True)
    out.emit('s11cMechanicalSignS11bThicknessWithoutBulk', sp.cancel(thickness-bulk))
    out.emit('s11cMechanicalSignS11bStiffnessCoefficient', sp.cancel(sp.diff(thickness, b.k_W, b.e_W)))
    out.emit('s11cMechanicalSignS11bInertiaCoefficient', sp.cancel(sp.diff(thickness, b.mu_W, b.e_W)))
    dependencies = {}
    for stage in ('c1', 'c2'):
        _, export_pins, imports = export_data(ROOT/f'scripts/S11c_{stage}_exports.py')
        dependencies[stage] = {'directImportKeys':imports, 'buildInputPins':export_pins,
            'pinResidual':{name:int(sha(locate(name)) != expected) for name,expected in export_pins.items()}}
    out.emit('s11cMechanicalSignDownstreamDependencies', dependencies)
    out.emit('s11cMechanicalSignSourceStabilityResidual',
             {name:int(sha(ROOT/name) != expected) for name,expected in source_pins.items()}, residual=True)
    resources = {'wallSeconds':time.monotonic()-STARTED,
                 'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    out.emit('s11cMechanicalSignWallSeconds', sp.Float(resources['wallSeconds']), units=tuple(b.DIM_T))
    out.emit('s11cMechanicalSignPeakRssKiB', resources['peakRssKiB'])
    args.inventory.write_text(json.dumps({'sourcePins':source_pins, 'case':args.case,
        'resources':resources, 'records':out.records}, indent=2)+'\n')


if __name__ == '__main__':
    run()
