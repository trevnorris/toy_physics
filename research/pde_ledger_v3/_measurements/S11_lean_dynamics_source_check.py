#!/usr/bin/env python3
"""E4: compact D2 native density, variational and modal identification.

Execute selected original SymPy helpers only, with their symbolic inputs; no
registry, emitters, production audit, exports, or S11c invocation. Translation
comparisons are symbolic CAS evidence, not kernel-certified parsing.
"""
import ast
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
from itertools import combinations_with_replacement
import json
from pathlib import Path
import sympy as sp

BASE = Path(__file__).resolve().parents[1]
SOURCE = BASE/'scripts/S11_stray_longitudinal_sympy_audit.py'
WOLFRAM = BASE/'mathematica/S11_stray_longitudinal_mathematica_audit.wl'
REPORT = BASE/'_measurements/S11_lean_dynamics_source_checks.json'
NAMES = {'Term','PackageBuild','compute_q9','q9_vector','q9_row_to_poly','matrix_from_rows',
         'derivative_placeholders','u_functions','stiffness_densities','coordinate_substitution',
         'to_coordinate','euler_lagrange_from_placeholders','route_a_matrix','route_b_matrix',
         'period_average','coefficient_ordering','package_build','q9_v5'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def native_environment(source):
    nodes = [n for n in ast.parse(source).body
             if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in NAMES]
    assert {n.name for n in nodes} == NAMES
    env = {'sp':sp,'dataclass':dataclass,'combinations_with_replacement':combinations_with_replacement,
           'QG_ALL':tuple(sp.Symbol(f'g_{i}',real=True) for i in range(1,26)),
           'X_ALL':tuple(sp.Symbol(f'x{i}',real=True) for i in range(1,6))}
    for name in ['beta','t','phase','omegaSquared']:
        env[name]=sp.Symbol(name,real=True)
    for name in ['rho_br','mu_R','B_comp','mu_br','s','s_rho','c_s0']:
        env[name]=sp.Symbol(name,positive=True)
    exec(compile(ast.Module(nodes,[]),str(SOURCE),'exec'),env)
    return env


def zero(expr):
    if isinstance(expr,sp.MatrixBase):
        return all(zero(x) for x in expr)
    return sp.expand(expr.doit()) == 0


def objects(env):
    native=env['compute_q9'](2)
    g,v=env['derivative_placeholders'](2)
    beta=env['beta']; k=sp.Matrix(sp.symbols('k1 k2',real=True)); a=sp.Matrix(sp.symbols('a1 a2',real=True))
    pd=native['PD_DENSITY_PLACEHOLDER']
    delta=sp.expand(env['package_build']('XFORM_EXTRA',2,native).lagrangian-
                    env['package_build']('MAIN',2,native).lagrangian)
    xs,u=env['u_functions'](2); x,y=xs
    expected_pd=(g[0,0]+g[1,1])*(g[0,1]-g[1,0])
    div=sp.diff(u[0],x)+sp.diff(u[1],y)
    skew=sp.diff(u[1],x)-sp.diff(u[0],y)
    lean_el=beta/2*sp.Matrix([sp.diff(skew,x)-sp.diff(div,y),
                             sp.diff(div,x)+sp.diff(skew,y)])
    native_el=sp.Matrix(env['euler_lagrange_from_placeholders'](delta,2,g,v))
    # V5 is presented in the native V1 basis. Use its actual RREF coordinates
    # to select the odd combination; do not assume it occupies a particular row.
    basis=native['V1_BASIS']; pivots=basis.rref()[1]
    variables=native['QG_VARIABLES']; pairs=tuple(combinations_with_replacement(range(4),2))
    pdrow=env['q9_vector'](native['PD_POLY'],variables,pairs)
    coefficients=[pdrow[p] for p in pivots]
    assert zero(sum((coefficients[i]*native['V1_POLYS'][i] for i in range(len(coefficients))),sp.S.Zero)-native['PD_POLY'])
    v5=env['q9_v5'](2,native)
    v5odd=sum((coefficients[i]*sp.Matrix(v5[i]) for i in range(len(coefficients))),sp.zeros(2,1))
    turn=sp.Matrix([-k[1],k[0]])
    lean_matrix=-beta/2*(k*turn.T+turn*k.T)
    route_a,cosine=env['route_a_matrix'](delta,2,g,v,tuple(k),tuple(a))
    route_b,average=env['route_b_matrix'](delta,2,g,v,tuple(k),tuple(a))
    lean_modal=-beta/2*(k.dot(a))*(turn.dot(a))
    identities={
        'actual_PD_equals_oddPairing':zero(pd-expected_pd),
        'actual_action_increment_sign_and_half':zero(delta+beta/2*expected_pd),
        'native_EL_equals_negative_Lean_EL':zero(native_el+lean_el),
        'native_V5_odd_combination_scale':zero(native_el+beta/2*v5odd),
        'route_A_equals_negative_Lean_operator':zero(route_a+lean_matrix),
        'route_B_equals_half_Lean_operator':zero(route_b-lean_matrix/2),
        'averaged_density_equals_half_modalAction':zero(average-lean_modal/2),
        'cosine_factor':cosine==sp.cos(env['phase']),
    }
    return identities,{'native_PD':str(pd),'actual_delta':str(delta),
        'native_EL':str(native_el),'Lean_EL_convention':str(lean_el),
        'V5_odd_coefficients':list(map(str,coefficients)),'V5_odd':str(v5odd),
        'Lean_modal_matrix':str(lean_matrix),'route_A':str(route_a),'route_B':str(route_b),
        'average':str(average)},(beta,k,lean_matrix)


def main():
    before={str(p.relative_to(BASE)):sha(p) for p in [SOURCE,WOLFRAM]}
    source=SOURCE.read_text();env=native_environment(source)
    identities,detail,(beta,k,M)=objects(env)
    assert all(identities.values()),identities
    controls={}
    for name,old,new,required in [
        ('wrong_action_sign','Term(beta / 2, "P_D", pd_density)',
         'Term(-beta / 2, "P_D", pd_density)','actual_action_increment_sign_and_half'),
        ('missing_action_half','Term(beta / 2, "P_D", pd_density)',
         'Term(beta, "P_D", pd_density)','route_A_equals_negative_Lean_operator')]:
        assert source.count(old)==1
        bad,_,_=objects(native_environment(source.replace(old,new)))
        assert not bad[required],(name,bad)
        controls[name]={'required_failed_identity':required,'identity_results':bad}
    # Explicit admissible and exceptional loci; the exhaustive statement is Lean's.
    positives={}
    for b,kk in [(0,(0,0)),(0,(1,0)),(2,(0,0)),(2,(1,0)),(-2,(3,4))]:
        subs={beta:b,k[0]:kk[0],k[1]:kk[1]};value=M.subs(subs)
        r=sp.Matrix([-kk[1],kk[0]]);q=sp.Matrix(kk)
        cross=(r.T*value*q)[0]
        assert (cross!=0)==(b!=0 and kk!=(0,0))
        positives[str((b,kk))]={'matrix':str(value),'cross_pairing':str(cross)}
    # Source-only Wolfram connection: record the actual normalization/route
    # sites. These guards do not pretend to execute or certify Wolfram output.
    wl=WOLFRAM.read_text()
    anchors=[
        '"XFORM_EXTRA", {{muR/2, "curl"}, {bComp/2, "div"}, {beta/2, "pd"}}',
        'lagrangianJet = Total[kineticTermsJet] - Total[stiffnessTermsJet];',
        'eulerOperatorForGradientDensity[polynomial_, gradientVariables_, frame_] := Module[',
        'Integrate[Expand[planeLagrangian], {phaseVariable, 0, 2 Pi}]/(2 Pi)']
    assert all(anchor in wl for anchor in anchors)
    assert before=={str(p.relative_to(BASE)):sha(p) for p in [SOURCE,WOLFRAM]}
    record={'status':'PASS','checked_utc':datetime.now(timezone.utc).isoformat(),
        'source_sha256':before,'instrument_sha256':sha(Path(__file__)),
        'selected_definitions':sorted(NAMES),'identities':identities,'objects':detail,
        'wolfram_source_only_anchors':anchors,
        'instrument_mutations':controls,'positive_locus_checks':positives,
        'translation':'G_ij=partial_i u_j; Lean spacetime (t,x1,x2), native function arguments (x1,x2,t); beta constant real. Lean EL=-div(momentum), native EL=+div(momentum). Route A=-M, Route B=M/2, average=modalAction/2.',
        'boundary':'Selected SymPy helper execution and exact symbolic comparisons. No Wolfram execution or production export/comparator clearance; no S11c files written.'}
    REPORT.write_text(json.dumps(record,indent=2)+'\n')


if __name__=='__main__':
    main()
