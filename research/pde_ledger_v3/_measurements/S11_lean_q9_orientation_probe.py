#!/usr/bin/env python3
"""Read-only D2 counterexample for the next Lean contract's source reconnaissance.

Execute five original Q9 helpers with real coordinate symbols. No production
audit, registry, export, Wolfram or S11c execution, and no native source edits.
"""
import ast
from datetime import datetime, timezone
import hashlib
from itertools import combinations_with_replacement
import json
from pathlib import Path
import sympy as sp

BASE = Path(__file__).resolve().parents[1]
SOURCE = BASE / 'scripts/S11_stray_longitudinal_sympy_audit.py'
NAMES = {'compute_q9', 'q9_vector', 'q9_row_to_poly',
         'matrix_from_rows', 'derivative_placeholders'}


def main():
    nodes = [n for n in ast.parse(SOURCE.read_text()).body
             if isinstance(n, ast.FunctionDef) and n.name in NAMES]
    assert {n.name for n in nodes} == NAMES
    coordinates = tuple(sp.Symbol(f'g_{i}', real=True) for i in range(1, 26))
    env = {'sp': sp, 'QG_ALL': coordinates,
           'combinations_with_replacement': combinations_with_replacement}
    exec(compile(ast.Module(nodes, []), str(SOURCE), 'exec'), env)
    native = env['compute_q9'](2)
    variables = coordinates[:4]
    G = sp.Matrix(2, 2, variables)
    J = sp.Matrix([[0, 1], [-1, 0]])
    delta = J * G - G * J
    R = sp.Matrix([[sp.Rational(3, 5), sp.Rational(-4, 5)],
                   [sp.Rational(4, 5), sp.Rational(3, 5)]])
    assert R.T * R == sp.eye(2) and R.det() == 1

    def infinitesimal(poly):
        return sp.expand(sum(sp.diff(poly, v) * delta[i // 2, i % 2]
                             for i, v in enumerate(variables)))

    def rotation_residual(poly):
        return sp.expand(poly.xreplace(dict(zip(variables, list(R * G * R.T)))) - poly)

    basis = native['V1_BASIS']
    polynomials = native['V1_POLYS']
    pairs = tuple(combinations_with_replacement(range(4), 2))
    trace_square = sp.trace(G) ** 2
    trace_row = sp.Matrix([env['q9_vector'](trace_square, variables, pairs)])
    rows = sp.Matrix([env['q9_vector'](infinitesimal(variables[p] * variables[q]),
                                     variables, pairs) for p, q in pairs])
    corrected = sp.Matrix.hstack(*rows.T.nullspace()).T
    corrected_polys = [env['q9_row_to_poly'](corrected.row(i), native['MONOMIAL_ORDERING'])
                       for i in range(corrected.rows)]
    assert infinitesimal(trace_square) == rotation_residual(trace_square) == 0
    assert basis.rank() == 4 and basis.col_join(trace_row).rank() == 5
    assert any(rotation_residual(poly) != 0 for poly in polynomials)
    assert corrected.rows == 4
    assert all(infinitesimal(poly) == rotation_residual(poly) == 0 for poly in corrected_polys)
    sample = dict(zip(variables, [1, 0, 0, 0]))
    witness = polynomials[0]
    before = witness.subs(sample)
    after = (witness + rotation_residual(witness)).subs(sample)
    assert before == 1 and after == sp.Rational(481, 625)
    record = {
        'status': 'CONFIRMED_SOURCE_FIDELITY_DISCREPANCY',
        'checked_utc': datetime.now(timezone.utc).isoformat(),
        'scope': 'Read-only native D2 Q9 counterexample; no complete census claim or engine repair',
        'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'instrument_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'selected_original_functions': sorted(NAMES),
        'coordinate_binding': 'g_1,...,g_25 as real symbols; registry initialization not executed',
        'native_SO_dimension': int(native['V1_DIM']),
        'native_O_dimension': int(native['V2_DIM']),
        'native_SO_basis': list(map(str, polynomials)),
        'native_infinitesimal_residuals': [str(infinitesimal(p)) for p in polynomials],
        'native_finite_rotation_residuals': [str(sp.factor(rotation_residual(p))) for p in polynomials],
        'proper_rotation': str(R), 'witness_G': 'diag(1,0)',
        'witness_polynomial': str(witness), 'before': str(before), 'after': str(after),
        'trace_square_rotation_residual': str(rotation_residual(trace_square)),
        'native_basis_rank': int(basis.rank()),
        'rank_after_adding_trace_square': int(basis.col_join(trace_row).rank()),
        'transposed_control_dimension': corrected.rows,
        'transposed_control_polynomials': list(map(str, corrected_polys)),
        'transposed_control_rotation_residuals': [str(rotation_residual(p)) for p in corrected_polys],
        'limits': 'A counterexample disproves the emitted invariant span. Passing one rotation does not prove full group invariance. D3-D5 and downstream effects are not measured here.',
        'initial_probe_note': 'An earlier inline assertion compared an unexpanded zero expression structurally; expanding residuals fixed that instrument error. It was not counted as mathematical failure.',
        'homogeneous_contract': 'H1-H4 excluded Q9 and checked MAIN independence from its sentinel; no homogeneous theorem is reopened.',
    }
    (BASE / '_measurements/S11_lean_q9_orientation_probe.json').write_text(
        json.dumps(record, indent=2) + '\n')


if __name__ == '__main__':
    main()
