#!/usr/bin/env python3
"""H1: execute the native D3 MAIN constructor and two compact modal routes.

Select original definitions and the three needed upstream scalar records, never
run an audit or write exports. References are handwritten translations of Lean;
hashes record provenance, and independent review supplies statement fidelity.
Wolfram is source-inspected and pinned, not re-executed by this instrument.
"""
import ast
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import sympy as sp
from sympy.core.symbol import Str

BASE = Path(__file__).resolve().parents[1]
PY = BASE / 'scripts/S11_stray_longitudinal_sympy_audit.py'
WL = BASE / 'mathematica/S11_stray_longitudinal_mathematica_audit.wl'
UPSTREAM = BASE / 'scripts/S10_exports.py'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    upstream_tree = ast.parse(UPSTREAM.read_text())
    upstream_env = {'sp': sp, 'Str': Str}
    restore_nodes = [n for n in upstream_tree.body if
                     (isinstance(n, ast.FunctionDef) and n.name == '_restore') or
                     (isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and
                      t.id == '_RELATIONALS' for t in n.targets))]
    assert len(restore_nodes) == 2
    exec(compile(ast.Module(restore_nodes, []), str(UPSTREAM), 'exec'), upstream_env)
    ledger = next(n.value for n in upstream_tree.body if isinstance(n, ast.Assign)
                  and any(isinstance(t, ast.Name) and t.id == '_LEDGER' for t in n.targets))
    incoming = {}
    upstream_ranges = []
    for key, record in zip(ledger.keys, ledger.values):
        if isinstance(key, ast.Constant) and key.value in {'rho_br', 'mu_R', 'omegaSquared'}:
            value = next(v for k, v in zip(record.keys, record.values)
                         if isinstance(k, ast.Constant) and k.value == 'value')
            incoming[key.value] = {'value': eval(compile(ast.Expression(value), str(UPSTREAM), 'eval'),
                                                upstream_env)}
            upstream_ranges.append([value.lineno, value.end_lineno])
    assert {'rho_br', 'mu_R'} <= incoming.keys()
    variables = set('CLASS_TAGS DECLARED_SYMBOLS rho_br mu_R B_comp mu_br beta s s_rho c_s0 '
                    'omegaSquared phase t X_ALL K_ALL A_ALL'.split())
    definitions = set('register_symbol declared_symbol derivative_placeholders stiffness_densities '
                      'package_build coefficient_ordering route_a_matrix period_average '
                      'route_b_matrix Term PackageBuild'.split())
    selected, found = [], set()
    for node in ast.parse(PY.read_text()).body:
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)) and node.name in definitions:
            selected.append(node)
            found.add(node.name)
        elif isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            names = {part.id for target in targets for part in ast.walk(target) if isinstance(part, ast.Name)}
            if names and names <= variables:
                selected.append(node)
                found |= names
    assert found == variables | definitions, found ^ (variables | definitions)
    env = {'sp': sp, 'dataclass': dataclass, 'INCOMING_LEDGER': incoming, '__name__': __name__}
    exec(compile(ast.Module(selected, []), str(PY), 'exec'), env)
    g, velocity = env['derivative_placeholders'](3)
    sentinel = sp.Symbol('UNUSED_INVARIANT_CENSUS_INPUT', real=True)
    build = env['package_build']('MAIN', 3, {'PD_DENSITY_PLACEHOLDER': sentinel})
    assert sentinel not in build.lagrangian.free_symbols
    rho, mu, B, z = [env[n] for n in ['rho_br', 'mu_R', 'B_comp', 'omegaSquared']]
    k = sp.Matrix(env['K_ALL'][:3]); a = sp.Matrix(env['A_ALL'][:3]); K = k.dot(k)
    curl = sp.Matrix([g[1,2]-g[2,1], g[2,0]-g[0,2], g[0,1]-g[1,0]])
    density = rho/2*sum(v**2 for v in velocity)-mu/2*curl.dot(curl)-B/2*sp.trace(g)**2
    operator = (rho*z-mu*K)*sp.eye(3)+(mu-B)*k*k.T
    checks = []

    def check(name, residual, expected_zero=True):
        entries = list(residual) if isinstance(residual, sp.MatrixBase) else [residual]
        values = [sp.factor(v) for v in entries]
        zero = all(v == 0 for v in values)
        assert zero == expected_zero, (name, values)
        checks.append({'name': name, 'expected_zero': expected_zero, 'observed_zero': zero,
                       'residual': str(values)})

    check('native_MAIN_action', build.lagrangian-density)
    ma, phase = env['route_a_matrix'](build.lagrangian, 3, g, velocity, tuple(k), tuple(a))
    mb, averaged = env['route_b_matrix'](build.lagrangian, 3, g, velocity, tuple(k), tuple(a))
    check('native_A_is_negative_Lean_operator', ma+operator)
    check('native_B_is_half_Lean_operator', 2*mb-operator)
    check('native_phase_average_normalization', averaged-(a.T*operator*a)[0]/4)
    check('native_B_determinant', mb.det()-(rho*z-mu*K)**2*(rho*z-B*K)/8)
    wrong = build.lagrangian+B*sp.trace(g)**2
    wrong_a = env['route_a_matrix'](wrong, 3, g, velocity, tuple(k), tuple(a))[0]
    check('reversed_compression_sign_detected_in_action', wrong-density, False)
    check('reversed_compression_sign_detected_in_operator', wrong_a+operator, False)
    check('missing_phase_half_detected', mb-operator, False)
    wl = WL.read_text()
    anchors = [
        '"MAIN" | "XKIN_ANISO", {{muR/2, "curl"}, {bComp/2, "div"}}',
        'divDensity[gradient_] := Tr[gradient]^2;',
        '{<|"Factor" -> rhoBr/2, "DensityJet" -> Total[velocityJet^2]|>}',
        'lagrangianJet = Total[kineticTermsJet] - Total[stiffnessTermsJet];',
        'D[D[lagrangianJet, velocityJet[[component]]] /. jetRules,',
        'Integrate[Expand[planeLagrangian], {phaseVariable, 0, 2 Pi}]/(2 Pi)',
        'D[averagedLagrangian, amplitudes[[row]], amplitudes[[column]]]',
    ]
    for anchor in anchors:
        assert wl.count(anchor) == 1, anchor
    result = {'status': 'PASS', 'scope': 'H1 compact D3 MAIN action/operator correspondence',
              'limits': 'Original selected SymPy code executed; Wolfram source inspection only. '
                        'References are handwritten; hashes do not certify translation or future edits. '
                        'No invariant census, full audit, export/comparator or S11c calculation executed.',
              'sympy_version': sp.__version__, 'checks': checks,
              'native_symbol_bindings': {n: sp.srepr(env[n]) for n in ['rho_br','mu_R','B_comp','omegaSquared']},
              'selected_original_source_ranges': [[n.lineno,n.end_lineno] for n in selected],
              'selected_upstream_scalar_ranges': upstream_ranges,
              'unused_PD_sentinel_absent_from_MAIN': True, 'route_A_stripped_phase': str(phase),
              'wolfram_literal_anchors': anchors,
              'source_sha256': {str(p.relative_to(BASE)):sha(p) for p in [PY,WL,UPSTREAM,Path(__file__),
                  BASE/'lean/s11/S11Homogeneous/Action.lean']}}
    (BASE/'_measurements/S11_lean_source_checks.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__ == '__main__':
    main()
