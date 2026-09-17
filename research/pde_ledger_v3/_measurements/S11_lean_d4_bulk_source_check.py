#!/usr/bin/env python3
"""D4C.4: compact native D4 V5 and bulk-operator identification; no production runs."""
import ast
from datetime import datetime, timezone
import hashlib
from itertools import combinations_with_replacement
import json
from pathlib import Path
import sympy as sp
BASE=Path(__file__).resolve().parents[1]
SOURCE=BASE/'scripts/S11_stray_longitudinal_sympy_audit.py'
WOLFRAM=BASE/'mathematica/S11_stray_longitudinal_mathematica_audit.wl'
REPORT=BASE/'_measurements/S11_lean_d4_bulk_source_checks.json'
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

def objects(env,native):
    g,velocity=env['derivative_placeholders'](4)
    F=g-g.T;P=F[0,1]*F[2,3]-F[0,2]*F[1,3]+F[0,3]*F[1,2]
    forms=[sp.trace(g)**2,sp.trace(g*g),sp.trace(g*g.T),P]
    pairs=tuple(combinations_with_replacement(range(16),2))
    variables=tuple(g);target=sp.Matrix([env['q9_vector'](f,variables,pairs) for f in forms])
    native_basis=native['V1_BASIS'];assert native_basis.shape==(4,136)
    basis_map=sp.Matrix.vstack(*[target.T.gauss_jordan_solve(native_basis.row(i).T)[0].T for i in range(4)])
    assert basis_map.det()!=0 and zero(basis_map*target-native_basis)
    coeff=sp.Matrix(sp.symbols('a b c beta',real=True));a,b,c,beta=coeff
    native_weights=basis_map.T.inv()*coeff
    qg=native['QG_VARIABLES'];replacement={qg[i]:variables[i] for i in range(16)}
    native_forms=[sp.expand(p.xreplace(replacement)) for p in native['V1_POLYS']]
    actual_Q=sp.expand(sum(native_weights[i]*native_forms[i] for i in range(4)))
    xs,u=env['u_functions'](4)
    div=sum(sp.diff(u[i],xs[i]) for i in range(4))
    grad_div=sp.Matrix([sp.diff(div,xs[i]) for i in range(4)])
    lap=sp.Matrix([sum(sp.diff(u[i],x,2) for x in xs) for i in range(4)])
    lean_el=(a+b)*grad_div+c*lap
    v5=[sp.Matrix(row) for row in env['q9_v5'](4,native)]
    combined_v5=sum((native_weights[i]*v5[i] for i in range(4)),sp.zeros(4,1))
    native_el=sp.Matrix(env['euler_lagrange_from_placeholders'](-actual_Q/2,4,g,velocity))
    identities={'full_native_span':zero(basis_map*target-native_basis),
      'actual_density_normalization':zero(actual_Q-sum(coeff[i]*forms[i] for i in range(4))),
      'V5_all_basis_combination':zero(combined_v5-2*lean_el),
      'actual_L_native_EL_equals_negative_Lean_EL':zero(native_el+lean_el)}
    per_basis=[]
    for i in range(4):
        row=basis_map.row(i)
        per_basis.append(zero(v5[i]-2*((row[0]+row[1])*grad_div+row[2]*lap)))
    identities['V5_every_actual_basis_element']=all(per_basis)
    k=sp.Matrix(sp.symbols('k1 k2 k3 k4',real=True));amp=sp.Matrix(sp.symbols('p1 p2 p3 p4',real=True))
    replacements={}
    for i in range(4):
        for j in range(4):
            for ell in range(4):
                replacements[sp.diff(u[i],xs[j],xs[ell])]=-k[j]*k[ell]*amp[i]
    M=-c*k.dot(k)*sp.eye(4)-(a+b)*k*k.T
    identities['native_EL_modal_sign']=zero(native_el.doit().xreplace(replacements)+M*amp)
    homogeneous=-c*k.dot(k)*sp.eye(4)+(c-(a+b+c))*k*k.T
    identities['homogeneous_mu_c_B_a_b_c']=zero(M-homogeneous)
    momenta=sp.Matrix(4,4,lambda i,j:sp.diff(-actual_Q/2,g[i,j]))
    dual=sp.Matrix(4,4,lambda i,j:sp.diff(P,g[i,j]))
    expected=-(a*sp.trace(g)*sp.eye(4)+b*g.T+c*g)-beta*dual/2
    identities['actual_full_momentum']=zero(momenta-expected)
    actualP=sp.expand(native['PD_POLY'].xreplace(replacement))
    identities['native_odd_normalization']=native['V6_DIM']==1 and zero(actualP-P)
    identities['zero_time_momentum']=all(sp.diff(-actual_Q/2,q)==0 for q in velocity)
    response=sp.Matrix([[0,0,1,0],[1,1,0,0]])
    identities['response_rank_two_nullity_two']=response.rank()==2 and len(response.nullspace())==2
    current=sp.Matrix([sum(u[i]*sp.diff(u[j],xs[j])-u[j]*sp.diff(u[i],xs[j]) for j in range(4)) for i in range(4)])
    G=sp.Matrix(4,4,lambda i,j:sp.diff(u[j],xs[i]))
    identities['even_current_divergence']=zero(sum(sp.diff(current[i],xs[i]) for i in range(4))-div**2+sp.trace(G*G))
    witness={q:0 for q in g};witness.update({g[0,1]:1,g[2,3]:1})
    assert P.subs(witness)==1
    assert momenta[0,1].subs({a:0,b:0,c:0,beta:2}).subs(witness)==-1
    assert not zero(momenta-(-(a*sp.trace(g)*sp.eye(4)+b*g.T+c*g)-beta*dual))
    assert not zero(momenta-(-(a*sp.trace(g)*sp.eye(4)+b*g.T+c*g)+beta*dual/2))
    # Target-side sensitivity controls supplement native sign/factor mutations.
    assert not zero(M-(-c*k.dot(k)*sp.eye(4)+(c-a)*k*k.T))
    assert zero(lean_el.subs({a:1,b:-1,c:0}))
    assert not zero(lean_el.subs({a:1,b:1,c:0}))
    return identities,{'native_basis_to_trace_forms':str(basis_map),
      'general_native_basis_weights':str(native_weights),'native_V5':list(map(str,v5)),
      'per_basis_matches':per_basis,'Lean_EL':str(lean_el),'native_EL_of_L':str(native_el),
      'modal_matrix':str(M),'normalization':'Q=a trace^2+b trace(GG)+c trace(GG^T)+beta P, L=-Q/2; native EL=-Lean EL.',
      'actual_momenta':str(momenta),'response_map':str(response),'nullspace':list(map(str,response.nullspace()))}

def main():
    before={str(p.relative_to(BASE)):sha(p) for p in [SOURCE,WOLFRAM]}
    source=SOURCE.read_text();env=environment(source);native=env['compute_q9'](4)
    identities,details=objects(env,native)
    assert all(identities.values()),identities
    old='system.append(sp.expand((time_coord + space_coord).xreplace(coord_sub)))'
    assert source.count(old)==1
    controls={}
    for name,new in [('native_EL_sign','system.append(sp.expand((-time_coord - space_coord).xreplace(coord_sub)))'),
                     ('native_EL_factor','system.append(sp.expand((2*time_coord + 2*space_coord).xreplace(coord_sub)))')]:
        bad,_=objects(environment(source.replace(old,new)),native)
        assert not bad['V5_every_actual_basis_element'] and not bad['actual_L_native_EL_equals_negative_Lean_EL'],bad
        controls[name]=bad
    wl=WOLFRAM.read_text()
    anchors=['eulerOperatorForGradientDensity[polynomial_, gradientVariables_, frame_] := Module[','Transpose[actionRows]']
    assert all(a in wl for a in anchors)
    assert before=={str(p.relative_to(BASE)):sha(p) for p in [SOURCE,WOLFRAM]}
    record={'status':'PASS','checked_utc':datetime.now(timezone.utc).isoformat(),
      'instrument_sha256':sha(Path(__file__)),'source_sha256':before,'selected_definitions':sorted(NAMES),
      'identities':identities,'objects':details,'native_mutations':controls,
      'target_controls':{'wrong_B_map_rejected':True,'null_even_positive':True,'nonnull_even_positive':True,'odd_density_nonzero_positive':True,'odd_momentum_nonzero_positive':True,'odd_momentum_factor_rejected':True,'odd_momentum_sign_rejected':True},
      'wolfram_source_only_anchors':anchors,
      'limits':'Selected original D4 Q9/V5/EL helper execution; no production drivers, exports, S11c or Wolfram execution. Exact symbolic translation checks, not a kernel-certified CAS bridge.'}
    REPORT.write_text(json.dumps(record,indent=2)+'\n')
if __name__=='__main__':main()
