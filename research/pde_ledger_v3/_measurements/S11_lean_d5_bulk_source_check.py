#!/usr/bin/env python3
"""D5B.4: actual D5 Q9/V5/action comparison with selected helpers, no production driver."""
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
REPORT=BASE/'_measurements/S11_lean_d5_bulk_source_checks.json'
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
    g,velocity=env['derivative_placeholders'](5)
    forms=[sp.trace(g)**2,sp.trace(g*g),sp.trace(g*g.T)]
    pairs=tuple(combinations_with_replacement(range(25),2))
    variables=tuple(g);target=sp.Matrix([env['q9_vector'](f,variables,pairs) for f in forms])
    basis=native['V1_BASIS'];assert basis.shape==(3,325)
    basis_map=sp.Matrix.vstack(*[target.T.gauss_jordan_solve(basis.row(i).T)[0].T for i in range(3)])
    assert basis_map.det()!=0 and zero(basis_map*target-basis)
    coeff=sp.Matrix(sp.symbols('a b c',real=True));a,b,c=coeff
    weights=basis_map.T.inv()*coeff
    qg=native['QG_VARIABLES'];replacement={qg[i]:variables[i] for i in range(25)}
    native_forms=[sp.expand(p.xreplace(replacement)) for p in native['V1_POLYS']]
    actual_Q=sp.expand(sum(weights[i]*native_forms[i] for i in range(3)))
    xs,u=env['u_functions'](5)
    div=sum(sp.diff(u[i],xs[i]) for i in range(5))
    grad_div=sp.Matrix([sp.diff(div,xs[i]) for i in range(5)])
    lap=sp.Matrix([sum(sp.diff(u[i],x,2) for x in xs) for i in range(5)])
    lean_el=(a+b)*grad_div+c*lap
    v5=[sp.Matrix(row) for row in env['q9_v5'](5,native)]
    combined_v5=sum((weights[i]*v5[i] for i in range(3)),sp.zeros(5,1))
    native_el=sp.Matrix(env['euler_lagrange_from_placeholders'](-actual_Q/2,5,g,velocity))
    identities={'full_native_span':zero(basis_map*target-basis),
      'actual_density_normalization':zero(actual_Q-sum(coeff[i]*forms[i] for i in range(3))),
      'V5_all_basis_combination':zero(combined_v5-2*lean_el),
      'actual_L_native_EL_equals_negative_Lean_EL':zero(native_el+lean_el)}
    per_basis=[zero(v5[i]-2*((basis_map[i,0]+basis_map[i,1])*grad_div+basis_map[i,2]*lap)) for i in range(3)]
    identities['V5_every_actual_basis_element']=all(per_basis)
    k=sp.Matrix(sp.symbols('k1:6',real=True));amp=sp.Matrix(sp.symbols('p1:6',real=True))
    replacements={sp.diff(u[i],xs[j],xs[ell]):-k[j]*k[ell]*amp[i]
                  for i in range(5) for j in range(5) for ell in range(5)}
    M=-c*k.dot(k)*sp.eye(5)-(a+b)*k*k.T
    identities['native_EL_modal_sign']=zero(native_el.doit().xreplace(replacements)+M*amp)
    homogeneous=-c*k.dot(k)*sp.eye(5)+(c-(a+b+c))*k*k.T
    identities['homogeneous_mu_c_B_a_b_c']=zero(M-homogeneous)
    momenta=sp.Matrix(5,5,lambda i,j:sp.diff(-actual_Q/2,g[i,j]))
    expected=-(a*sp.trace(g)*sp.eye(5)+b*g.T+c*g)
    identities['actual_full_momentum']=zero(momenta-expected)
    identities['zero_time_momentum']=all(sp.diff(-actual_Q/2,q)==0 for q in velocity)
    response=sp.Matrix([[0,0,1],[1,1,0]])
    identities['response_rank_two_nullity_one']=response.rank()==2 and len(response.nullspace())==1
    current=sp.Matrix([sum(u[i]*sp.diff(u[j],xs[j])-u[j]*sp.diff(u[i],xs[j]) for j in range(5)) for i in range(5)])
    G=sp.Matrix(5,5,lambda i,j:sp.diff(u[j],xs[i]))
    identities['current_divergence']=zero(sum(sp.diff(current[i],xs[i]) for i in range(5))-div**2+sp.trace(G*G))
    identities['no_odd_density']=native['V6_DIM']==0 and native['PD_POLY']==0
    witness={q:0 for q in g};witness.update({g[0,0]:1,g[1,1]:1})
    null={a:1,b:-1,c:0}
    wrong_momentum=-(a*sp.trace(g)*sp.eye(5)+b*g+c*g)
    controls={
      'momentum_transpose_rejected':not zero(momenta-wrong_momentum),
      'wrong_B_map_rejected':not zero(M-(-c*k.dot(k)*sp.eye(5)+(c-a)*k*k.T)),
      'null_positive':zero(lean_el.subs(null)),
      'nonnull_positive':not zero(lean_el.subs({a:1,b:1,c:0})),
      'null_density_nonzero':(-actual_Q/2).subs(null).subs(witness)==-1,
      'momentum_nonzero':momenta[0,0].subs({a:1,b:0,c:0}).subs(witness)==-2,
      'current_factor_rejected':not zero(2*(sp.trace(g)**2-sp.trace(g*g))-(sp.trace(g)**2-sp.trace(g*g))),
      'fifth_transverse':M.subs({a:2,b:3,c:7,**{k[i]:int(i==4) for i in range(5)}})[0,0]==-7,
      'fifth_longitudinal':M.subs({a:2,b:3,c:7,**{k[i]:int(i==4) for i in range(5)}})[4,4]==-12}
    assert all(controls.values()),controls
    return identities,controls,{'native_basis_to_trace_forms':str(basis_map),
      'general_native_basis_weights':str(weights),'native_V5':list(map(str,v5)),
      'per_basis_matches':per_basis,'Lean_EL':str(lean_el),'native_EL_of_L':str(native_el),
      'modal_matrix':str(M),'actual_momenta':str(momenta),'response_map':str(response),
      'nullspace':list(map(str,response.nullspace())),
      'normalization':'Q=a trace^2+b trace(GG)+c trace(GG^T), L=-Q/2; native EL=-Lean EL.'}

def main():
    before={str(p.relative_to(BASE)):sha(p) for p in [SOURCE,WOLFRAM]}
    source=SOURCE.read_text();env=environment(source);native=env['compute_q9'](5)
    identities,target_controls,details=objects(env,native)
    assert all(identities.values()),identities
    old='system.append(sp.expand((time_coord + space_coord).xreplace(coord_sub)))'
    assert source.count(old)==1
    controls={}
    for name,new in [('native_EL_sign','system.append(sp.expand((-time_coord - space_coord).xreplace(coord_sub)))'),
                     ('native_EL_factor','system.append(sp.expand((2*time_coord + 2*space_coord).xreplace(coord_sub)))')]:
        bad,_,_=objects(environment(source.replace(old,new)),native)
        assert not bad['V5_every_actual_basis_element'] and not bad['actual_L_native_EL_equals_negative_Lean_EL'],bad
        controls[name]=bad
    anchors=['eulerOperatorForGradientDensity[polynomial_, gradientVariables_, frame_] := Module[','Transpose[actionRows]']
    assert all(a in WOLFRAM.read_text() for a in anchors)
    assert before=={str(p.relative_to(BASE)):sha(p) for p in [SOURCE,WOLFRAM]}
    record={'status':'PASS','checked_utc':datetime.now(timezone.utc).isoformat(),
      'instrument_sha256':sha(Path(__file__)),'source_sha256':before,'selected_definitions':sorted(NAMES),
      'identities':identities,'objects':details,'native_mutations':controls,'target_controls':target_controls,
      'wolfram_source_only_anchors':anchors,
      'limits':'Selected original D5 Q9/V5/EL helpers only; no production driver, exports, S11c or Wolfram execution. Translation checked outside the Lean kernel.'}
    REPORT.write_text(json.dumps(record,indent=2)+'\n')
if __name__=='__main__':main()
