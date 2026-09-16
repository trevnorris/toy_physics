#!/usr/bin/env python3
"""K4: compact native D3 V5 and bulk-operator identification; no production runs."""
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
REPORT=BASE/'_measurements/S11_lean_d3_bulk_source_checks.json'
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

def objects(env):
    native=env['compute_q9'](3);g,velocity=env['derivative_placeholders'](3)
    forms=[sp.trace(g)**2,sp.trace(g*g),sp.trace(g*g.T)]
    pairs=tuple(combinations_with_replacement(range(9),2))
    variables=tuple(g);target=sp.Matrix([env['q9_vector'](f,variables,pairs) for f in forms])
    native_basis=native['V1_BASIS'];assert native_basis.shape==(3,45)
    basis_map=sp.Matrix.vstack(*[target.T.gauss_jordan_solve(native_basis.row(i).T)[0].T for i in range(3)])
    assert basis_map.det()!=0 and zero(basis_map*target-native_basis)
    coeff=sp.Matrix(sp.symbols('a b c',real=True));a,b,c=coeff
    native_weights=basis_map.T.inv()*coeff
    qg=native['QG_VARIABLES'];replacement={qg[i]:variables[i] for i in range(9)}
    native_forms=[sp.expand(p.xreplace(replacement)) for p in native['V1_POLYS']]
    actual_Q=sp.expand(sum(native_weights[i]*native_forms[i] for i in range(3)))
    xs,u=env['u_functions'](3)
    div=sum(sp.diff(u[i],xs[i]) for i in range(3))
    grad_div=sp.Matrix([sp.diff(div,xs[i]) for i in range(3)])
    lap=sp.Matrix([sum(sp.diff(u[i],x,2) for x in xs) for i in range(3)])
    lean_el=(a+b)*grad_div+c*lap
    v5=[sp.Matrix(row) for row in env['q9_v5'](3,native)]
    combined_v5=sum((native_weights[i]*v5[i] for i in range(3)),sp.zeros(3,1))
    native_el=sp.Matrix(env['euler_lagrange_from_placeholders'](-actual_Q/2,3,g,velocity))
    identities={'full_native_span':zero(basis_map*target-native_basis),
      'actual_density_normalization':zero(actual_Q-sum(coeff[i]*forms[i] for i in range(3))),
      'V5_all_basis_combination':zero(combined_v5-2*lean_el),
      'actual_L_native_EL_equals_negative_Lean_EL':zero(native_el+lean_el)}
    per_basis=[]
    for i in range(3):
        row=basis_map.row(i)
        per_basis.append(zero(v5[i]-2*((row[0]+row[1])*grad_div+row[2]*lap)))
    identities['V5_every_actual_basis_element']=all(per_basis)
    k=sp.Matrix(sp.symbols('k1 k2 k3',real=True));amp=sp.Matrix(sp.symbols('p1 p2 p3',real=True))
    replacements={}
    for i in range(3):
        for j in range(3):
            for ell in range(3):
                replacements[sp.diff(u[i],xs[j],xs[ell])]=-k[j]*k[ell]*amp[i]
    M=-c*k.dot(k)*sp.eye(3)-(a+b)*k*k.T
    identities['native_EL_modal_sign']=zero(native_el.doit().xreplace(replacements)+M*amp)
    homogeneous=-c*k.dot(k)*sp.eye(3)+(c-(a+b+c))*k*k.T
    identities['homogeneous_mu_c_B_a_b_c']=zero(M-homogeneous)
    # Target-side sensitivity controls supplement native sign/factor mutations.
    assert not zero(M-(-c*k.dot(k)*sp.eye(3)+(c-a)*k*k.T))
    assert zero(lean_el.subs({a:1,b:-1,c:0}))
    assert not zero(lean_el.subs({a:1,b:1,c:0}))
    return identities,{'native_basis_to_trace_forms':str(basis_map),
      'general_native_basis_weights':str(native_weights),'native_V5':list(map(str,v5)),
      'per_basis_matches':per_basis,'Lean_EL':str(lean_el),'native_EL_of_L':str(native_el),
      'modal_matrix':str(M),'normalization':'Q=a trace^2+b trace(GG)+c trace(GG^T), L=-Q/2; native EL=-Lean EL.'}

def main():
    before={str(p.relative_to(BASE)):sha(p) for p in [SOURCE,WOLFRAM]}
    source=SOURCE.read_text();identities,details=objects(environment(source))
    assert all(identities.values()),identities
    old='system.append(sp.expand((time_coord + space_coord).xreplace(coord_sub)))'
    assert source.count(old)==1
    controls={}
    for name,new in [('native_EL_sign','system.append(sp.expand((-time_coord - space_coord).xreplace(coord_sub)))'),
                     ('native_EL_factor','system.append(sp.expand((2*time_coord + 2*space_coord).xreplace(coord_sub)))')]:
        bad,_=objects(environment(source.replace(old,new)))
        assert not bad['V5_every_actual_basis_element'] and not bad['actual_L_native_EL_equals_negative_Lean_EL'],bad
        controls[name]=bad
    wl=WOLFRAM.read_text()
    anchors=['eulerOperatorForGradientDensity[polynomial_, gradientVariables_, frame_] := Module[','Transpose[actionRows]']
    assert all(a in wl for a in anchors)
    assert before=={str(p.relative_to(BASE)):sha(p) for p in [SOURCE,WOLFRAM]}
    record={'status':'PASS','checked_utc':datetime.now(timezone.utc).isoformat(),
      'instrument_sha256':sha(Path(__file__)),'source_sha256':before,'selected_definitions':sorted(NAMES),
      'identities':identities,'objects':details,'native_mutations':controls,
      'target_controls':{'wrong_B_map_rejected':True,'null_1_minus1_0':True,'nonnull_1_1_0':True},
      'wolfram_source_only_anchors':anchors,
      'limits':'Selected original D3 Q9/V5/EL helper execution; no production drivers, exports, S11c or Wolfram execution. Exact symbolic translation checks, not a kernel-certified CAS bridge.'}
    REPORT.write_text(json.dumps(record,indent=2)+'\n')
if __name__=='__main__':main()
