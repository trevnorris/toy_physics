#!/usr/bin/env python3
"""Compact S9 native-constructor check; never runs an audit or writes exports.

Execute only selected original SymPy declarations/functions. Compare the MAIN
action and two native modal routes with the mathematical contract, symbolically.
The handwritten reference is a reviewed translation, not a Lean-certified Python
interpreter. Wolfram connection is source inspection, explicitly not a fresh run.
"""
import ast
import hashlib
import json
from pathlib import Path
import sympy as sp

BASE = Path(__file__).resolve().parents[1]
PY = BASE / 'scripts/S9_light_requires_shear_sympy_audit.py'
WL = BASE / 'mathematica/S9_light_requires_shear_mathematica_audit.wl'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    names = set('t x y z coordinates spatial_coordinates u_functions a_symbols '
                'b_symbols omega omegaSquared k_input rho_br mu_R identity3 '
                'main_action position input_wavevector phase_exponent phase plane_wave'.split())
    functions = {'curl_of', 'construct_curl_action', 'euler_lagrange_residual',
                 'route_a', 'route_b'}
    selected, assigned = [], set()
    for node in ast.parse(PY.read_text()).body:
        if isinstance(node, ast.Assign):
            targets = {item.id for target in node.targets for item in ast.walk(target)
                       if isinstance(item, ast.Name)}
            if targets <= names:
                assigned |= targets
                selected.append(node)
        elif isinstance(node, ast.FunctionDef) and node.name in functions:
            assigned.add(node.name)
            selected.append(node)
    assert assigned == names | functions, assigned ^ (names | functions)
    env = {'sp': sp}
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(PY), 'exec'), env)
    rho, mu, s = (env[n] for n in ('rho_br', 'mu_R', 'omegaSquared'))
    k = sp.Matrix(env['k_input'])
    jet = sp.Matrix(4, 3, lambda j, i: sp.Symbol(f'J{j}{i}', real=True))
    substitutions = {sp.diff(u, c): jet[j, i] for j, c in enumerate(env['coordinates'])
                     for i, u in enumerate(env['u_functions'])}
    curl = sp.Matrix([jet[2, 2] - jet[3, 1], jet[3, 0] - jet[1, 2],
                      jet[1, 1] - jet[2, 0]])
    expected_action = rho / 2 * sum(jet[0, i]**2 for i in range(3)) - mu / 2 * curl.dot(curl)
    expected_operator = rho * s * sp.eye(3) - mu * (k.dot(k) * sp.eye(3) - k * k.T)
    action = env['main_action']
    checks = []

    def check(name, residual, should_vanish=True):
        entries = list(residual) if isinstance(residual, sp.MatrixBase) else [residual]
        simplified = [sp.factor(e) for e in entries]
        vanishes = all(e == 0 for e in simplified)
        assert vanishes == should_vanish, (name, simplified)
        checks.append({'name': name, 'expected_zero': should_vanish,
                       'observed_zero': vanishes, 'residual': str(simplified)})

    check('native_MAIN_action_equals_Lean_jet_reference',
          action.xreplace(substitutions) - expected_action)
    check('native_route_A_equals_contract_operator', env['route_a'](action)[2] - expected_operator)
    check('native_paired_route_B_equals_contract_operator', env['route_b'](action) - expected_operator)
    check('native_squared_frequency_determinant',
          env['route_a'](action)[2].det() - rho * s * (rho * s - mu * k.dot(k))**2)
    wrong_action = env['construct_curl_action'](rho * sp.eye(3), mu, stiffness_sign=1)
    check('wrong_shear_sign_action_rejected', wrong_action.xreplace(substitutions) - expected_action, False)
    check('wrong_shear_sign_operator_rejected', env['route_a'](wrong_action)[2] - expected_operator, False)

    wl = WL.read_text()
    wl_anchors = [
        'coordinates = {t, x, y, z};',
        'velocityVector = D[fieldVector, t];',
        'curlVector = Curl[fieldVector, spaceCoordinates];',
        'mainLagrangian = rhoBr velocityVector.velocityVector/2 - muR curlVector.curlVector/2;',
        'wavePhase = Exp[I (waveVector.spaceCoordinates - omega t)];',
    ]
    for anchor in wl_anchors:
        assert wl.count(anchor) == 1, anchor
    result = {
        'status': 'PASS',
        'scope': 'S9 MAIN action and modal operator; exact symbolic SymPy checks, Wolfram source anchors',
        'limits': 'References are handwritten and independently reviewed. No Lean execution of CAS, no fresh Wolfram run, no export/comparator certification.',
        'sympy_version': sp.__version__, 'checks': checks,
        'selected_original_source_ranges': [[n.lineno, n.end_lineno] for n in selected],
        'wolfram_source_anchors': wl_anchors,
        'source_sha256': {str(p.relative_to(BASE)): digest(p) for p in (
            PY, WL, Path(__file__), BASE / 'lean/s9/S9Pilot/Action.lean',
            BASE / 'lean/s9/S9Pilot/PlaneWave.lean')},
    }
    (BASE / '_measurements/S9_lean_source_checks.json').write_text(
        json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
