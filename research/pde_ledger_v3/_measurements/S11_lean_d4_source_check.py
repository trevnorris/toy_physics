#!/usr/bin/env python3
"""D4.3: compact actual D4 Q9 span identification; no production driver/export.

Lean establishes the full-group theorem. This native SymPy check identifies
its four forms with actual Q9 V1/V2 and the one-dimensional odd V6 space, before EL quotients.
Wolfram evidence here is explicitly source inspection, not an engine execution.
"""
import ast
from datetime import datetime, timezone
import hashlib
from itertools import combinations_with_replacement, permutations
import json
from pathlib import Path
import sympy as sp

BASE=Path(__file__).resolve().parents[1]
SOURCE=BASE/'scripts/S11_stray_longitudinal_sympy_audit.py'
WOLFRAM=BASE/'mathematica/S11_stray_longitudinal_mathematica_audit.wl'
NAMES={'compute_q9','q9_vector','q9_row_to_poly','matrix_from_rows','derivative_placeholders'}
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def native_environment(text):
 nodes=[n for n in ast.parse(text).body if isinstance(n,ast.FunctionDef) and n.name in NAMES]
 assert {n.name for n in nodes}==NAMES
 env={'sp':sp,'QG_ALL':tuple(sp.Symbol(f'g_{i}',real=True) for i in range(1,26)),
      'combinations_with_replacement':combinations_with_replacement}
 exec(compile(ast.Module(nodes,[]),str(SOURCE),'exec'),env)
 return env

def main():
 before={str(p.relative_to(BASE)):sha(p) for p in [SOURCE,WOLFRAM]}
 text=SOURCE.read_text();env=native_environment(text);native=env['compute_q9'](4)
 variables=env['QG_ALL'][:16];pairs=tuple(combinations_with_replacement(range(16),2));G=sp.Matrix(4,4,variables)
 P=(G[0,1]-G[1,0])*(G[2,3]-G[3,2])-(G[0,2]-G[2,0])*(G[1,3]-G[3,1])+(G[0,3]-G[3,0])*(G[1,2]-G[2,1])
 forms=[sp.trace(G)**2,sp.trace(G*G),sp.trace(G*G.T),P]
 rows=lambda fs: sp.Matrix([env['q9_vector'](f,variables,pairs) for f in fs]) if fs else sp.zeros(0,136)
 expected=rows(forms)
 assert native['MONOMIAL_ORDERING']==sp.Tuple(*(variables[i]*variables[j] for i,j in pairs))
 groups={}
 for group,key,target,count in [('SO','V1',expected,4),('O','V2',rows(forms[:3]),3),('odd','V6',rows([P]),1)]:
  actual=native[key+'_BASIS']
  assert actual.cols==136 and actual.rows==actual.rank()==target.rank()==count
  assert actual.col_join(target).rank()==count
  assert actual.rref()[0]==target.rref()[0]
  groups[group]={'dimension':count,'native_basis':str(actual),
   'native_polynomials':[str(env['q9_row_to_poly'](actual.row(i),native['MONOMIAL_ORDERING'])) for i in range(actual.rows)],
   'same_rref':True,'stacked_rank':count}
 reflection=sp.diag(-1,1,1,1);reflected=reflection*G*reflection.T
 substitutions=dict(zip(variables,list(reflected)))
 assert all(sp.expand(forms[i].xreplace(substitutions)-forms[i])==0 for i in range(3))
 assert sp.expand(P.xreplace(substitutions)+P)==0
 assert sp.expand(native['PD_POLY']-P)==0
 epsilon=sum(sp.LeviCivita(*q)*G[q[0],q[1]]*G[q[2],q[3]] for q in permutations(range(4)))
 assert sp.expand(epsilon-2*P)==0
 actual=native['V1_BASIS'];reflection_rows=rows([sp.expand(f.xreplace(substitutions)) for f in native['V1_POLYS']])
 assert native['V6_OPERATOR'].T*actual==reflection_rows
 # A wrong four-dimensional span must not pass by having the right count.
 wrong=rows([forms[0],forms[1],forms[2],G[0,0]**2])
 assert wrong.rank()==4 and wrong.col_join(expected).rank()>4
 # Deliberately restore the old generator orientation only in memory.
 fixed='lie_equations.extend(matrix_from_rows(action_rows, len(monomials)).T.tolist())'
 assert text.count(fixed)==1
 bad=native_environment(text.replace(fixed,'lie_equations.extend(action_rows)'))['compute_q9'](4)
 assert bad['V1_BASIS'].col_join(expected).rank()>4
 assert 'Transpose[actionRows]' in WOLFRAM.read_text()
 assert before=={str(p.relative_to(BASE)):sha(p) for p in [SOURCE,WOLFRAM]}
 record={'status':'PASS','checked_utc':datetime.now(timezone.utc).isoformat(),
  'source_sha256':before,'instrument_sha256':sha(Path(__file__)),
  'selected_functions':sorted(NAMES),'coordinate_convention':'all real row-major g_1,...,g_16; G_ij = partial_i u_j; no EL quotient',
  'contract_polynomials':[str(sp.expand(f)) for f in forms],'native_spans':groups,
  'reflection_operator_matches_native_basis':True,'PD_equals_contract_P':True,'epsilon_contraction_equals_2P':True,
  'controls':{'same_count_wrong_span_rejected':True,'wrong_orientation_counts':[int(bad['V1_DIM']),int(bad['V2_DIM']),int(bad['V6_DIM'])],
              'wrong_orientation_span_rejected':True},
  'wolfram_link':'source inspection of Transpose[actionRows] only; no engine execution claimed',
  'limits':'D4 Q9 V1/V2/V6 only; no V5/EL, divergence, dynamics, D5, production rerun/export/comparator or S11c work'}
 (BASE/'_measurements/S11_lean_d4_source_checks.json').write_text(json.dumps(record,indent=2)+'\n')

if __name__=='__main__':main()
