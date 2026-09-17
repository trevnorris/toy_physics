#!/usr/bin/env python3
"""D4B.4: actual native odd density, momenta and null operator; no production run."""
import ast
from datetime import datetime,timezone
import hashlib,json
from itertools import combinations_with_replacement
from pathlib import Path
import sympy as sp
BASE=Path(__file__).resolve().parents[1]
SOURCE=BASE/'scripts/S11_stray_longitudinal_sympy_audit.py'
WOLFRAM=BASE/'mathematica/S11_stray_longitudinal_mathematica_audit.wl'
NAMES={'compute_q9','q9_vector','q9_row_to_poly','matrix_from_rows','derivative_placeholders',
       'u_functions','coordinate_substitution','to_coordinate','euler_lagrange_from_placeholders','q9_v5'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def zero(e):
 if isinstance(e,sp.MatrixBase):return all(zero(x) for x in e)
 return sp.expand(e.doit())==0
def environment(source):
 nodes=[n for n in ast.parse(source).body if isinstance(n,ast.FunctionDef) and n.name in NAMES]
 assert {n.name for n in nodes}==NAMES
 env={'sp':sp,'combinations_with_replacement':combinations_with_replacement,
      'QG_ALL':tuple(sp.Symbol(f'g_{i}',real=True) for i in range(1,26)),
      'X_ALL':tuple(sp.Symbol(f'x{i}',real=True) for i in range(1,6)),
      't':sp.Symbol('t',real=True)}
 exec(compile(ast.Module(nodes,[]),str(SOURCE),'exec'),env)
 return env

def main():
 before={str(p.relative_to(BASE)):sha(p) for p in [SOURCE,WOLFRAM]}
 text=SOURCE.read_text();env=environment(text);native=env['compute_q9'](4)
 G,velocity=env['derivative_placeholders'](4);F=G-G.T
 P=F[0,1]*F[2,3]-F[0,2]*F[1,3]+F[0,3]*F[1,2]
 M=sp.Matrix([[0,F[2,3],-F[1,3],F[1,2]],[-F[2,3],0,F[0,3],-F[0,2]],
              [F[1,3],-F[0,3],0,F[0,1]],[-F[1,2],F[0,2],-F[0,1],0]])
 beta=sp.Symbol('beta',real=True)
 replace=dict(zip(native['QG_VARIABLES'],list(G)))
 actualP=sp.expand(native['PD_POLY'].xreplace(replace))
 assert native['V6_DIM']==1 and zero(actualP-P)
 L=-beta*actualP/2
 momenta=sp.Matrix(4,4,lambda i,j:sp.diff(L,G[i,j]))
 assert zero(momenta+beta*M/2) and any(not zero(x) for x in momenta)
 assert all(zero(sp.diff(L,v)) for v in velocity)
 xs,u=env['u_functions'](4);coord=env['to_coordinate'](M,4,G,velocity)
 K=sp.Matrix([sum(u[j]*coord[i,j] for j in range(4))/2 for i in range(4)])
 actual_coord_P=env['to_coordinate'](actualP,4,G,velocity)
 divK=sum(sp.diff(K[i],xs[i]) for i in range(4))
 assert zero(divK-actual_coord_P)
 assert all(zero(sum(sp.diff(coord[i,j],xs[i]) for i in range(4))) for j in range(4))
 assert zero(sum(G[i,j]*M[i,j] for i in range(4) for j in range(4))-2*P)
 nativeEL=sp.Matrix(env['euler_lagrange_from_placeholders'](L,4,G,velocity));assert zero(nativeEL)
 # Identify just the native odd combination of the existing V1/V5 objects.
 weights=native['V1_BASIS'].T.gauss_jordan_solve(native['V6_BASIS'].row(0).T)[0]
 v5=[sp.Matrix(row) for row in env['q9_v5'](4,native)]
 odd_v5=sum((weights[i]*v5[i] for i in range(len(v5))),sp.zeros(4,1));assert zero(odd_v5)
 # Zero bulk alone cannot identify sign or scale: nonzero momenta must detect both.
 assert not zero(momenta-beta*M/2)
 assert not zero(momenta+beta*M)
 assert not zero(2*divK-actual_coord_P)
 assert not zero(actualP)
 badL=-(actualP+G[0,0]**2)/2
 badEL=sp.Matrix(env['euler_lagrange_from_placeholders'](badL,4,G,velocity))
 assert not zero(badEL) and zero(badEL-sp.Matrix([-sp.diff(u[0],xs[0],2),0,0,0]))
 old='system.append(sp.expand((time_coord + space_coord).xreplace(coord_sub)))'
 assert text.count(old)==1
 broken=environment(text.replace(old,'system.append(sp.S.Zero)'))
 assert zero(sp.Matrix(broken['euler_lagrange_from_placeholders'](badL,4,G,velocity)))
 # This would disagree with the explicit non-null target, so the native failure is detected.
 wrong_index=sum(sp.diff(coord[i,0],xs[0]) for i in range(4))
 assert not zero(wrong_index)
 witness={G[i,j]:int((i,j) in [(0,1),(2,3)]) for i in range(4) for j in range(4)}
 assert P.subs(witness)==1 and momenta[0,1].subs(witness).subs(beta,2)==-1
 assert 'Transpose[actionRows]' in WOLFRAM.read_text()
 assert before=={str(p.relative_to(BASE)):sha(p) for p in [SOURCE,WOLFRAM]}
 record={'status':'PASS','checked_utc':datetime.now(timezone.utc).isoformat(),'instrument_sha256':sha(Path(__file__)),
  'source_sha256':before,'selected_definitions':sorted(NAMES),
  'identities':{'actual_PD_equals_reviewed_P':True,'actual_momenta_equal_minus_beta_over_two_M':True,'time_momentum_zero':True,
    'dual_divergence_zero_all_components':True,'contraction_equals_2P':True,'explicit_current_divergence_equals_P':True,
    'actual_native_EL_zero':True,'actual_native_V5_odd_combination_zero':True},
  'objects':{'P':str(sp.expand(P)),'dual_matrix':str(M),'native_momenta':str(momenta),'odd_weights_in_native_V1':str(weights),
    'native_odd_V5':str(odd_v5),'native_EL':str(nativeEL),'current':str(K),
    'normalization':'L=-beta P_D/2, P_D=P; Lean EL=-native EL. Nonzero momenta check the otherwise invisible sign/factor.'},
  'controls':{'wrong_momentum_sign_rejected':True,'wrong_momentum_factor_rejected':True,'wrong_current_factor_rejected':True,
    'zero_density_claim_rejected':True,'wrong_derivative_index_rejected':True,'added_G00_square_nonnull':str(badEL),
    'native_always_zero_operator_detected_on_nonnull_density':True,'density_witness':1,'momentum_witness_beta_2':-1},
  'wolfram_link':'source inspection of generator transpose only; no Wolfram execution',
  'limits':'D4 odd family only; selected native helper calls and exact symbolic translation. No full D4 bulk census, variable beta, D5, spectra/interfaces, production driver, exports, S11c or systematic CAS bridge.'}
 (BASE/'_measurements/S11_lean_d4_odd_source_checks.json').write_text(json.dumps(record,indent=2)+'\n')
if __name__=='__main__':main()
