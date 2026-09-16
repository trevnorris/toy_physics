#!/usr/bin/env python3
"""Compact D2 Q9 fidelity check: actual native spans, not transcript replication.

Run five original helpers only; no registry, production audit, or S11c jobs.
The symbolic span identities below are CAS evidence. The full-group theorem
and completeness are supplied separately by S11Invariants in Lean.
"""
import ast
from datetime import datetime, timezone
import hashlib
from itertools import combinations_with_replacement
import json
from pathlib import Path
import sympy as sp

BASE = Path(__file__).resolve().parents[1]
SOURCE = BASE/'scripts/S11_stray_longitudinal_sympy_audit.py'
WOLFRAM = BASE/'mathematica/S11_stray_longitudinal_mathematica_audit.wl'
NAMES = {'compute_q9', 'q9_vector', 'q9_row_to_poly',
         'matrix_from_rows', 'derivative_placeholders'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    before = sha(SOURCE)
    nodes = [n for n in ast.parse(SOURCE.read_text()).body
             if isinstance(n, ast.FunctionDef) and n.name in NAMES]
    assert {n.name for n in nodes} == NAMES
    coordinates = tuple(sp.Symbol(f'g_{i}', real=True) for i in range(1, 26))
    env = {'sp': sp, 'QG_ALL': coordinates,
           'combinations_with_replacement': combinations_with_replacement}
    exec(compile(ast.Module(nodes, []), str(SOURCE), 'exec'), env)
    native = env['compute_q9'](2)
    variables = coordinates[:4]
    a,b,c,d = variables
    G = sp.Matrix([[a,b],[c,d]])
    pairs = tuple(combinations_with_replacement(range(4), 2))
    t,s,x,y = a+d,b-c,a-d,b+c
    forms = [t*t,s*s,x*x+y*y,t*s]
    def rows(polys):
        return sp.Matrix([env['q9_vector'](p, variables, pairs) for p in polys])
    assert native['MONOMIAL_ORDERING'] == sp.Tuple(*(variables[i]*variables[j] for i,j in pairs))
    groups = {}
    for name,key,indices,count in [('SO','V1',[0,1,2,3],4),
                                   ('O','V2',[0,1,2],3),('odd','V6',[3],1)]:
        actual = native[key+'_BASIS']
        expected = rows([forms[i] for i in indices])
        assert actual.rows == actual.rank() == expected.rank() == count
        assert actual.col_join(expected).rank() == count
        # Row-space identity is a complete comparison of the two finite spans.
        assert actual.rref()[0] == expected.rref()[0]
        groups[name] = {'dimension':count,'native_basis':str(actual),
                        'basis_polynomials':[str(env['q9_row_to_poly'](actual.row(i),native['MONOMIAL_ORDERING']))
                                             for i in range(actual.rows)],
                        'contract_polynomials':[str(sp.expand(forms[i])) for i in indices],
                        'stacked_rank':count,'same_rref':True}
    # Deliberately undo just the orientation in memory: the same counts must
    # not pass a span check. This is an instrument control, not a Lean mutation.
    fixed = 'lie_equations.extend(matrix_from_rows(action_rows, len(monomials)).T.tolist())'
    assert SOURCE.read_text().count(fixed) == 1
    mutant_nodes = [n for n in ast.parse(SOURCE.read_text().replace(fixed,
                    'lie_equations.extend(action_rows)')).body
                    if isinstance(n, ast.FunctionDef) and n.name in NAMES]
    mutant_env = dict(env)
    exec(compile(ast.Module(mutant_nodes, []), 'wrong_orientation_control', 'exec'), mutant_env)
    wrong = mutant_env['compute_q9'](2)
    assert wrong['V1_DIM'] == 4 and wrong['V2_DIM'] == 3
    assert wrong['V1_BASIS'].col_join(rows(forms)).rank() > 4
    R = sp.Matrix([[sp.Rational(3,5),sp.Rational(-4,5)],
                   [sp.Rational(4,5),sp.Rational(3,5)]])
    p = a*a+b*c+d*d
    residual = sp.expand(p.xreplace(dict(zip(variables,list(R*G*R.T))))-p)
    assert residual.subs(dict(zip(variables,[1,0,0,0]))) == sp.Rational(-144,625)
    # Compact source-level Wolfram orientation link, not a Wolfram execution.
    wl = WOLFRAM.read_text()
    assert 'Transpose[actionRows]' in wl
    assert before == sha(SOURCE)
    record = {'status':'PASS','checked_utc':datetime.now(timezone.utc).isoformat(),
              'source_sha256':before,'wolfram_source_sha256':sha(WOLFRAM),
              'instrument_sha256':sha(Path(__file__)), 'selected_functions':sorted(NAMES),
              'coordinate_convention':'row-major real g_1,...,g_4; G_ij = partial_i u_j; no EL quotient',
              'native_spans':groups,
              'orientation_control':{'wrong_counts':[4,3],'wrong_span_rejected':True,
                                     'counterexample_residual':'-144/625'},
              'wolfram_link':'source inspection: Transpose[actionRows]; no engine execution claimed',
              'limits':'D2 Q9 V1/V2/V6 spans only. No D3-D5 census, V5/EL classification, production rerun, comparator/export or S11c work.'}
    (BASE/'_measurements/S11_lean_invariant_source_checks.json').write_text(json.dumps(record,indent=2)+'\n')


if __name__ == '__main__':
    main()
