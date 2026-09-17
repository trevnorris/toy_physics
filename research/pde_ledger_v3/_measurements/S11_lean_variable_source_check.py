#!/usr/bin/env python3
"""VC4: exact prescribed-profile identities for the previously identified local families.
Selected original Q9/coordinate helpers only; no production or S11c execution.
The variable-coefficient divergence is derived here, not attributed to the old
constant-coefficient native EL helper.
"""
from datetime import datetime,timezone
import hashlib,json,runpy
from pathlib import Path
import sympy as sp
B=Path(__file__).resolve().parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def zero(e):
 if isinstance(e,sp.MatrixBase):return all(zero(x) for x in e)
 return sp.expand(e.doit())==0

def main():
 sources=[B/'scripts/S11_stray_longitudinal_sympy_audit.py',B/'_measurements/S11_lean_d4_odd_source_check.py',
  B/'_measurements/S11_lean_d3_bulk_source_checks.json',B/'_measurements/S11_lean_d4_odd_source_checks.json']
 before={str(p.relative_to(B)):sha(p) for p in sources}
 helper=runpy.run_path(str(sources[1]));env=helper['environment'](sources[0].read_text())
 identities={};controls={};details={}
 for D in [3,4]:
  native=env['compute_q9'](D);G,V=env['derivative_placeholders'](D);xs,u=env['u_functions'](D)
  convert=lambda z:env['to_coordinate'](z,D,G,V)
  if D==3:
   a,b,c=sp.symbols('vc_a vc_b vc_c',real=True);coeff=[a,b,c]
   forms=sp.Matrix([sp.trace(G)**2,sp.trace(G*G),sp.trace(G*G.T)])
   from itertools import combinations_with_replacement
   pairs=tuple(combinations_with_replacement(range(9),2))
   target=sp.Matrix([env['q9_vector'](f,tuple(G),pairs) for f in forms])
   mapping=sp.Matrix.vstack(*[target.T.gauss_jordan_solve(native['V1_BASIS'].row(i).T)[0].T for i in range(3)])
   weights=mapping.T.inv()*sp.Matrix(coeff)
   polys=[p.xreplace(dict(zip(native['QG_VARIABLES'],list(G)))) for p in native['V1_POLYS']]
   Q=sp.expand(sum(weights[i]*polys[i] for i in range(3)))
   assert zero(Q-sum(coeff[i]*forms[i] for i in range(3)))
   funcs=[sp.Function(n)(*xs) for n in ['vc_a','vc_b','vc_c']];sub=dict(zip(coeff,funcs));af,bf,cf=funcs
   p=convert(sp.Matrix(D,D,lambda i,j:sp.diff(-Q/2,G[i,j]))).subs(sub)
   EL=sp.Matrix([-sum(sp.diff(p[i,j],xs[i]) for i in range(D)) for j in range(D)])
   div=sum(sp.diff(u[i],xs[i]) for i in range(D))
   correction=sp.Matrix([sp.diff(af,xs[j])*div+sum(sp.diff(bf,xs[i])*sp.diff(u[i],xs[j])+sp.diff(cf,xs[i])*sp.diff(u[j],xs[i]) for i in range(D)) for j in range(D)])
   constant=sp.Matrix([(af+bf)*sp.diff(div,xs[j])+cf*sum(sp.diff(u[j],x,2) for x in xs) for j in range(D)])
   J=sp.Matrix([u[i]*div-sum(u[j]*sp.diff(u[i],xs[j]) for j in range(D)) for i in range(D)])
   nullQ=convert(Q.subs({b:-a,c:0})).subs(a,af)
   weighted=sum(sp.diff(af*J[i],xs[i]) for i in range(D))
   gradpair=sum(sp.diff(af,xs[i])*J[i] for i in range(D))
   identities.update(D3_native_density=True,D3_EL_product=zero(EL-constant-correction),D3_weighted_null=zero(nullQ-weighted+gradpair))
   controls.update(D3_omit_gradient=not zero(EL-constant),D3_gradient_sign=not zero(EL-constant+correction),D3_omit_weighted_current_correction=not zero(nullQ-weighted))
   # VC4 explicitly requires sensitivity to derivative-index transposition.
   # Keep this separate from sign/omission controls, which need not detect it.
   wrong_b_correction=sp.Matrix([sp.diff(af,xs[j])*div+sum(sp.diff(bf,xs[i])*sp.diff(u[j],xs[i])+sp.diff(cf,xs[i])*sp.diff(u[j],xs[i]) for i in range(D)) for j in range(D)])
   index_witness={u[0]:xs[1],u[1]:0,u[2]:0,af:0,bf:xs[0],cf:0}
   correct_index_value=EL.subs(index_witness).doit()
   wrong_index_value=(constant+wrong_b_correction).subs(index_witness).doit()
   assert correct_index_value==sp.Matrix([0,1,0]) and wrong_index_value==sp.zeros(3,1)
   controls['D3_wrong_derivative_index']=not zero(EL-constant-wrong_b_correction)
   details['D3_index_witness']={'actual_EL':str(correct_index_value),'transposed_b_gradient_claim':str(wrong_index_value)}
   witness={u[0]:0,u[1]:xs[1],u[2]:0,af:xs[0],bf:-xs[0],cf:0}
   got=EL.subs(witness).doit();assert got==sp.Matrix([1,0,0])
   details['D3_variable_null_witness']=str(got)
  else:
   F=G-G.T;P=F[0,1]*F[2,3]-F[0,2]*F[1,3]+F[0,3]*F[1,2]
   actualP=native['PD_POLY'].xreplace(dict(zip(native['QG_VARIABLES'],list(G))));assert zero(actualP-P)
   beta=sp.Function('vc_beta')(*xs)
   M=sp.Matrix(D,D,lambda i,j:sp.diff(actualP,G[i,j]));Mc=convert(M)
   p=-beta*Mc/2
   EL=sp.Matrix([-sum(sp.diff(p[i,j],xs[i]) for i in range(D)) for j in range(D)])
   correct=sp.Matrix([sum(sp.diff(beta,xs[i])*Mc[i,j] for i in range(D))/2 for j in range(D)])
   K=sp.Matrix([sum(u[j]*Mc[i,j] for j in range(D))/2 for i in range(D)])
   Pc=convert(actualP);weighted=sum(sp.diff(beta*K[i],xs[i]) for i in range(D));gp=sum(sp.diff(beta,xs[i])*K[i] for i in range(D))
   identities.update(D4_native_PD=True,D4_EL_gradient=zero(EL-correct),D4_weighted_current=zero(beta*Pc-weighted+gp))
   controls.update(D4_zero_bulk_claim=not zero(EL),D4_factor=not zero(EL-2*correct),D4_sign=not zero(EL+correct),D4_omit_weighted_correction=not zero(beta*Pc-weighted))
   got=EL.subs({u[0]:0,u[1]:0,u[2]:0,u[3]:xs[2],beta:xs[0]}).doit();assert got==sp.Matrix([0,sp.Rational(1,2),0,0])
   details['D4_variable_beta_witness']=str(got)
 # Nonzero actual normal-slice integrals: p-=2, p+=5, test 1-x^2 on [-1,1].
 x=sp.Symbol('x',real=True);h=1-x*x
 interface=sp.integrate(2*sp.diff(h,x),(x,-1,0))+sp.integrate(5*sp.diff(h,x),(x,0,1))
 assert interface==-3
 identities['normal_slice_integral_jump']=True
 controls.update(interface_omission=interface!=0,interface_reversed_sign=interface!=3)
 details['normal_slice_jump']=str(interface)
 assert all(identities.values()) and all(controls.values())
 assert before=={str(p.relative_to(B)):sha(p) for p in sources}
 r={'status':'PASS','checked_utc':datetime.now(timezone.utc).isoformat(),'instrument_sha256':sha(Path(__file__)),
 'source_sha256':before,'identities':identities,'controls':controls,'witnesses':details,
 'limits':['Selected native Q9 densities and coordinate helpers only. Variable-profile total divergence derived by this instrument; old constant-coefficient EL helper is not used for this extension.',
 'Only the reviewed D3 family and D4 odd density. No full S11c operator identification, production driver, pinned export or Wolfram execution.',
 'Normal-slice interface witness; no multidimensional trace theorem or distribution product is claimed.']}
 (B/'_measurements/S11_lean_variable_source_checks.json').write_text(json.dumps(r,indent=2)+'\n')
if __name__=='__main__':main()
