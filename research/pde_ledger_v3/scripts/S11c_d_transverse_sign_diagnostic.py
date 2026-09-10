#!/usr/bin/env python3
"""Focused sign trace authorized by the user; no production rows are modified.

Print computed operands and residuals. Interpret them in the measurement report.
This homogeneous conservative diagnostic is not the S11c-d scattering engine.
The upstream/open rows are used for provenance tracing only.
"""
from __future__ import annotations

import argparse
import hashlib
import inspect
import pickle
import re
import resource
import time
from pathlib import Path

import sympy as sp
from sympy.core.symbol import Str

import S11c_d_mixing_scattering_sympy_audit as engine
from ledger_fold import load_model, assert_lookups_equal_manifest

ROOT = Path(__file__).resolve().parent.parent
IMPORT_KEYS = (
    'energy_basis_variable', 'slab_operator', 'transverse_dispersion',
    'energy_reexpressed_density', 'rho_br', 'mu_R', 'mu_S', 'epsilon_shape',
    'eta_bg', 'sigma_W', 'W_0', 'mu_W', 'omega', 'k',
    's11cc2ClosedSlabOperator',
)
INPUT_PATHS = tuple(ROOT / p for p in (
    'scripts/S11c_b_exports.py', 'scripts/S11c_c1_exports.py',
    'scripts/S11c_c2_exports.py', 'scripts/S11c_d_mixing_scattering_sympy_audit.py',
    'scripts/S11c_b_brane_operator_sympy_audit.py',
    'scripts/S11c_c2_selfenergy_fold_sympy_audit.py',
    'scripts/S11b_interface_coupling_law_sympy_audit.py',
    'scripts/ledger_fold.py', 'directives/S11b_SHARED_PHYSICS.md',
    'directives/S11c_b_SHARED_PHYSICS.md',
)) + (Path(__file__).resolve(),)
named = engine.named


def bind(fold):
    return {key: dict(fold[key]) for key in IMPORT_KEYS}


class Dimensions(engine.DimensionAnalysis):
    """Supplied units for this small conservative slice; inherited arithmetic."""

    def __init__(self, rows, x, t, q, h, standing):
        self.zero = (sp.S.Zero,) * 3
        self.known, self.unknown, self.constraints, self.solution = {}, {}, set(), {}
        units = {
            'rho_br': (-3, 0, 1), 'mu_R': (-1, -2, 1), 'mu_S': (-1, -2, 1),
            'W_0': (1, 0, 0), 'mu_W': (-3, 0, 1), 'omega': (0, -1, 0),
            'k': (-1, 0, 0), 'B_rho_3': (-1, -2, 1), 'C': (-2, -2, 1),
            'k_W': (-3, -2, 1), 'epsilon_shape': (0, 0, 0),
            'eta_bg': (0, 0, 0), 'sigma_W': (0, 0, 0),
        }
        for row in rows.values():
            for a in row['value'].atoms(sp.Symbol):
                if a.name in units:
                    self.known[a] = tuple(map(sp.Integer, units[a.name]))
        self.known.update({x: (1, 0, 0), t: (0, 1, 0), q.func: (1, 0, 0),
                           h.func: self.zero, standing.func: (1, 0, 0)})


class Records:
    def __init__(self, units, epsilon):
        self.units, self.epsilon, self.lines = units, epsilon, {}

    def emit(self, name, value, *, zero_unit=None):
        if name in self.lines:
            raise ValueError(('duplicate diagnostic key', name))
        self.lines[name] = inspect.currentframe().f_back.f_lineno
        body = engine.cas(value)
        metadata = []
        for path, leaf in engine.leaves(body):
            if isinstance(leaf, Str):
                continue
            unit = self.units.measure(leaf)
            support = tuple((int(g[0]), 0, 0) for g, coefficient in
                            sp.Poly(leaf, self.epsilon).terms() if coefficient != 0)
            metadata.append((path, {
                'DIMENSION_L_T_M': unit if unit is not None else
                    (zero_unit.get(path) if isinstance(zero_unit, dict) else zero_unit),
                'MULTIGRADE': support,
                'EPSILON_LAMBDA_SUPPORT': tuple((g[0], 0) for g in support),
            }))
        engine.emit('SIGN_' + name, {'VALUE': body, 'METADATA': metadata})


def value_cases(row):
    return {tuple(map(str, key)) if isinstance(key, sp.Tuple) else str(key):
            named(payload, 'VALUE') for key, payload in row['value']}


def jet_ansatz(expression, x, t, q, h, direction, polarization, *, thickness=False):
    """Map every displacement/scalar jet by differentiation of a field ansatz."""
    mapping = {}
    for a in expression.atoms(sp.Symbol):
        name = re.sub(r'^grad_theta_([123])$', r'theta_d\1', a.name)
        match = re.fullmatch(r'(u_[123]|theta|e_W)((?:_t{1,2})?(?:_?d[123])*)', name)
        if not match:
            continue
        field, suffix = match.groups()
        if thickness:
            base = {'theta': -h, 'e_W': h}.get(field, sp.S.Zero)
        else:
            base = q if field == 'u_' + str(polarization + 1) else sp.S.Zero
        spatial = [int(v) - 1 for v in re.findall(r'd([123])', suffix)]
        for i in spatial:
            base = sp.diff(base, x) if i == direction and not thickness else sp.S.Zero
        if '_tt' in suffix:
            base = sp.diff(base, t, 2)
        elif '_t' in suffix:
            base = sp.diff(base, t)
        mapping[a] = base
    return sp.expand(expression.xreplace(mapping))


def euler_balance(lagrangian, field, x, t, epsilon):
    """Euler variation with positive time-momentum derivative normalization."""
    divergence = sum(sp.diff(sp.diff(lagrangian, sp.diff(field, coordinate)), coordinate)
                     for coordinate in (t, x) if coordinate in field.args)
    return sp.expand((divergence - sp.diff(lagrangian, field)) / epsilon)


def closed_slice(value, x, t, q, direction, polarization):
    functions = tuple(sp.Function('s11cc2Field' + s) for s in
                      ('u1', 'u2', 'u3', 'theta', 'eW'))
    mapping = {f: (lambda *args, j=j: q.func(args[direction], args[-1])
                   if j == polarization else sp.S.Zero) for j, f in enumerate(functions)}
    projected = engine.dag_substitute(value, mapping)
    coords = {a: (t if a.name == 's11cc2Time' else x) for a in projected.free_symbols
              if a.name in ('s11cc2Time', 's11cc2X' + str(direction + 1))}
    return engine.map_leaves(projected.xreplace(coords), sp.expand)


def run(rows, inputs, *, cached=False, started=None):
    started = time.monotonic() if started is None else started
    eps = rows['epsilon_shape']['value']
    rho, mu_s = (rows[k]['value'] for k in ('rho_br', 'mu_S'))
    omega, k = (rows[n]['value'] for n in ('omega', 'k'))
    x, t = sp.symbols('s11cdSignPosition s11cdSignTime', real=True)
    q = sp.Function('s11cdSignDisplacement', real=True)(x, t)
    h = sp.Function('s11cdSignThickness', real=True)(t)
    standing = sp.Function('s11cdSignStandingAmplitude', real=True)(t)
    units = Dimensions(rows, x, t, q, h, standing)
    records = Records(units, eps)
    engine.emit('SIGN_INPUT_DIGESTS', inputs)
    engine.emit('SIGN_CACHED_DEVELOPMENT', sp.sympify(cached))
    engine.emit('SIGN_IMPORT_PROVENANCE', {key: rows[key]['step'] for key in IMPORT_KEYS})
    cut = {rows[n]['value']: sp.S.Zero for n in ('eta_bg', 'sigma_W')}
    # Apply the uniform control to source operands, never to a computed result.
    source_energy = value_cases(rows['energy_basis_variable'])
    source_open = value_cases(rows['slab_operator'])
    source_closed = value_cases(rows['s11cc2ClosedSlabOperator'])
    energies = {case: sp.Add(*(row[4] if str(row[0]) in ('W_BG', 'MU_R_BG') else row[3]
                               for row in body[1:])).xreplace(cut)
                for case, body in source_energy.items()}
    legacy_eq = rows['transverse_dispersion']['value'][0]
    legacy_bindings = {a: rows[a.name]['value'] for a in legacy_eq.free_symbols
                       if a.name in ('rho_br', 'mu_R', 'mu_S', 'omega', 'k')}
    records.emit('S11B_SYMBOL_ALIGNMENT_OPERANDS', tuple(legacy_bindings.items()))
    kinetic = eps**2 * rho * sp.diff(q, t)**2 / 2
    records.emit('KINETIC_ACTION_DENSITY', kinetic)
    for case, operator in source_open.items():
        suffix = '_'.join(case)
        open_u = named(named(operator, 'U_BODY_BALANCE'), 'EXPANDED').xreplace(cut)
        open_e = named(named(operator, 'E_W_BALANCE'), 'EXPANDED').xreplace(cut)
        closed = engine.memo_xreplace(source_closed[case], cut)
        for direction in range(3):
            for polarization in range(3):
                if direction == polarization:
                    continue  # Enumerate the six supplied transverse ansatz choices.
                label = suffix + f'_D{direction + 1}_P{polarization + 1}'
                energy = jet_ansatz(energies[case[0]], x, t, q, h, direction, polarization)
                action = kinetic - energy
                balance = euler_balance(action, q, x, t, eps)
                stored = sp.ImmutableMatrix([jet_ansatz(v, x, t, q, h, direction, polarization)
                                            for v in open_u])
                closed_rows = closed_slice(closed, x, t, q, direction, polarization)
                closed_u = sp.ImmutableMatrix(named(closed_rows, 'U'))
                scalar = stored[polarization]
                records.emit('STORED_ENERGY_' + label, energy)
                records.emit('ACTION_BALANCE_' + label, balance)
                records.emit('OPEN_U_' + label, stored, zero_unit=(-2, -2, 1))
                records.emit('CLOSED_U_' + label, closed_u, zero_unit=(-2, -2, 1))
                records.emit('CLOSED_MINUS_OPEN_' + label, closed_u - stored,
                             zero_unit=(-2, -2, 1))
                records.emit('CLOSED_SCALAR_ROWS_' + label,
                             (named(closed_rows, 'THETA'), named(closed_rows, 'E_W')),
                             zero_unit={(0,): (-3, -1, 1), (1,): (-1, -2, 1)})
                records.emit('ROW_MINUS_ACTION_' + label, scalar - balance,
                             zero_unit=(-2, -2, 1))
                records.emit('ROW_PLUS_ACTION_' + label, scalar + balance,
                             zero_unit=(-2, -2, 1))
                jets = (sp.diff(q, t, 2), sp.diff(q, x, 2))
                coefficients = sp.ImmutableMatrix([[sp.diff(row, v, eps) for v in jets]
                                                   for row in (scalar, balance)])
                records.emit('INERTIA_STIFFNESS_COEFFICIENTS_' + label, coefficients)
                scale = sp.cancel(coefficients[0, 1] / coefficients[1, 1])
                records.emit('STIFFNESS_ALIGNED_ROW_RESIDUAL_' + label,
                             scalar - scale * balance, zero_unit=(-2, -2, 1))
                # Compute both Fourier conventions directly from their ansatz.
                for time_sign in (-1, 1):
                    phase = sp.exp(sp.I * (k * x + time_sign * omega * t))
                    symbols = tuple(sp.cancel(row.subs(q, phase).doit() / (eps * phase))
                                    for row in (scalar, balance))
                    records.emit('HARMONIC_OPERANDS_' + label + '_T' + str(time_sign), symbols)
                    records.emit('FREQUENCY_SQUARED_' + label + '_T' + str(time_sign),
                                 tuple(sp.solve(row, omega**2)[0] for row in symbols))
                if (direction, polarization) == (0, 1):
                    sum_action_balance = euler_balance(kinetic + energy, q, x, t, eps)
                    records.emit('SUM_ACTION_BALANCE_' + suffix, sum_action_balance)
                    records.emit('STORED_PLUS_SUM_ACTION_RESIDUAL_' + suffix,
                                 scalar + sum_action_balance, zero_unit=(-2, -2, 1))
                    # S11b has a different mu_S energy normalization. The pure
                    # curl cut isolates the time/stiffness sign from that issue.
                    legacy_full = (legacy_eq.lhs - legacy_eq.rhs).xreplace(legacy_bindings)
                    records.emit('S11B_S11CB_STIFFNESS_OPERANDS_' + suffix,
                                 (legacy_full.subs(omega, 0), symbols[1].subs(omega, 0)))
                    records.emit('S11CB_MINUS_S11B_STIFFNESS_' + suffix,
                                 symbols[1].subs(omega, 0) - legacy_full.subs(omega, 0))
                    legacy = legacy_full.subs(mu_s, 0)
                    records.emit('S11B_PURE_CURL_OPERAND_' + suffix, legacy)
                    records.emit('S11B_ACTION_PURE_CURL_RESIDUAL_' + suffix,
                                 legacy - symbols[1].subs(mu_s, 0), zero_unit=(-3, -2, 1))
                    records.emit('S11B_STORED_PURE_CURL_RESIDUAL_' + suffix,
                                 symbols[0].subs(mu_s, 0) - legacy, zero_unit=(-3, -2, 1))
                    # Spatially averaged energy per volume of a real standing
                    # wave; integration uses dimensionless phase, not a typed mean.
                    phase_coordinate = sp.Symbol('s11cdSignSpatialPhase', real=True)
                    units.known[phase_coordinate] = units.zero
                    wave = standing * sp.cos(k * x)
                    total_energy = (kinetic + energy).subs(q, wave).doit()
                    mean_energy = sp.simplify(sp.integrate(total_energy.subs(x, phase_coordinate / k),
                                             (phase_coordinate, 0, 2 * sp.pi)) / (2 * sp.pi))
                    wave_rows = [sp.simplify(row.subs(q, wave).doit() / sp.cos(k * x))
                                 for row in (scalar, balance)]
                    accelerations = [sp.solve(row, sp.diff(standing, t, 2))[0] for row in wave_rows]
                    rates = [sp.factor(sp.diff(mean_energy, t).subs(sp.diff(standing, t, 2), acc))
                             for acc in accelerations]
                    records.emit('STANDING_WAVE_ENERGY_' + suffix, mean_energy)
                    records.emit('STANDING_WAVE_ACCELERATIONS_' + suffix, accelerations)
                    records.emit('STANDING_WAVE_ENERGY_RATES_' + suffix, rates,
                                 zero_unit=(-1, -3, 1))
        # Independently enforce the supplied no-flux k=0 virtual constraint in
        # the energy before varying. The open pressure slots are external data.
        thickness_energy = jet_ansatz(energies[case[0]], x, t, q, h, 0, 0, thickness=True)
        thickness_kinetic = eps**2 * rows['mu_W']['value'] * rows['W_0']['value']**2 * sp.diff(h, t)**2 / 2
        thickness_balance = euler_balance(thickness_kinetic - thickness_energy, h, x, t, eps)
        pressures = {a: sp.S.Zero for a in open_e.free_symbols if a.name in
                     ('delta_p_plus', 'delta_p_minus', 'Lambda_A_0', 'Lambda_V_0', 'Lambda_X_0')}
        thickness_stored = jet_ansatz(open_e.xreplace(pressures), x, t, q, h, 0, 0, thickness=True)
        records.emit('THICKNESS_ENERGY_' + suffix, thickness_energy)
        records.emit('THICKNESS_ACTION_BALANCE_' + suffix, thickness_balance)
        records.emit('THICKNESS_STORED_BALANCE_' + suffix, thickness_stored)
        records.emit('THICKNESS_ROW_MINUS_ACTION_' + suffix, thickness_stored - thickness_balance,
                     zero_unit=units.measure(thickness_balance))
    engine.emit('SIGN_DIMENSION_CONSTRAINT_RESIDUALS', sorted(units.constraints, key=sp.default_sort_key))
    engine.emit('SIGN_DIMENSION_UNRESOLVED', units.unknown)
    engine.emit('SIGN_EMISSION_LINES', records.lines)
    engine.emit('SIGN_RESOURCE_MEASUREMENTS', (time.monotonic() - started,
                                              resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    if units.constraints or units.unknown:
        raise ValueError('unresolved diagnostic dimensions')


if __name__ == '__main__':
    started = time.monotonic()
    parser = argparse.ArgumentParser()
    parser.add_argument('--development-cache', type=Path)
    args = parser.parse_args()
    digests = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in INPUT_PATHS}
    if args.development_cache:
        # This optional local cache is explicitly marked, never a final run.
        with args.development_cache.open('rb') as stream:
            rows = pickle.load(stream)
    else:
        fold, audit = load_model(ROOT / 'scripts/S11c_b_exports.py', ROOT / 'scripts/S11c_c1_exports.py',
                                 ROOT / 'scripts/S11c_c2_exports.py')
        witness = assert_lookups_equal_manifest(bind, fold, IMPORT_KEYS)
        rows = witness['result']
        engine.emit('SIGN_IMPORT_FOLD', audit)
        engine.emit('SIGN_IMPORT_LOOKUPS', sorted(witness['lookups']))
    run(rows, digests, cached=args.development_cache is not None, started=started)
