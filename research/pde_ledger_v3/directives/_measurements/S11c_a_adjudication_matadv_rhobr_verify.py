import sys; sys.path.insert(0,"research/pde_ledger_v3/scripts")
import sympy as sp
from sympy import Symbol,Function,Integer,Rational,Mul,Add,Pow,I,Float,Dummy
from sympy.core.function import AppliedUndef
import S11c_a_cross_engine_comparator as C
NS=dict(Symbol=Symbol,Function=Function,Integer=Integer,Rational=Rational,Mul=Mul,Add=Add,Pow=Pow,I=I,Float=Float,Dummy=Dummy,Tuple=sp.Tuple,oo=sp.oo,pi=sp.pi)
def parse(s):
    v=eval(s,{"__builtins__":{}},NS); return Add(*v) if isinstance(v,tuple) else v
HEADS=("delta_j_bulk","mu_theta","delta_rho","theta","currentXPerturbation","currentWPerturbation","currentWWave")
def collapse(e): return e.xreplace({a:Symbol(a.func.__name__) for a in e.atoms(AppliedUndef) if any(a.func.__name__.startswith(h) for h in HEADS)})
raw=None
import re as _re, os as _os
_OUT=_os.path.expanduser("~/.s11_build/comparator_run.out")
_grab=False
for _l in open(_OUT):
    if _l.startswith("CASE family=PROJECTION_STATIC_OPERAND key=(BRANCH=MATERIAL_ADVECTED, DENSITY=RHOBR_CONSTANT"): _grab=True; continue
    if _grab and _l.startswith("A_minus_B = "): raw="A_minus_B = "+_l[len("A_minus_B = "):]; break

amb=[l[len("A_minus_B = "):] for l in raw.splitlines() if l.startswith("A_minus_B = ")][0]
red=C.combine_bound_integrals(collapse(parse(amb)))
red=sp.factor_terms(sp.expand(red))
print("MAT×RHOBR PROJECTION_STATIC_OPERAND — reduced residual (collapse + integral linearity, NO identity sub):")
print(sp.srepr(red)[:1400])
print()
print("atoms:", sorted({a.name for a in red.atoms(Symbol)}))
# Now bind PY symbol to WL's full density-time (symbol := rho_br/W_bg*theta_t + advection) and check ->0.
W0,rho_br,theta_t,eta,w1,sw=map(Symbol,['W_0','rho_br','theta_t','eta_bg','w1_profile','sigma_W'])
u=[Symbol(f'u_{i}_t') for i in (1,2,3)]; d=[Symbol(f'w1_profile_d{i}') for i in (1,2,3)]
Wbg=Symbol('W_bg'); drho=Symbol('delta_rho_4D_bulk_t')
G_t=sum(u[i]*d[i] for i in range(3))
# spec-built full advected density-time (RHOBR: rho_bg=rho_br/W_bg): d/dt[ (rho_br/W_bg(chi)) (1+theta) ] first order
# = (rho_br/W_bg) theta_t  +  d/dt(rho_br/W_bg(chi));  chi=x-eps u => d/dt(1/W_bg(chi)) = -(-1/W_bg^2) grad W_bg . u = (grad W_bg . u)/W_bg^2
# grad W_bg = sigma_W * dw1 (per §2a),  1/W_bg ~ (1-eta w1)/W0 to 1st order
full = (rho_br/Wbg)*theta_t + rho_br*sw*G_t/Wbg**2
test = red.xreplace({drho: full})
test = test.replace(Pow(Wbg,Integer(-1)),(1-eta*w1)/W0).xreplace({Wbg:W0*(1+eta*w1)})
# first shape order truncation via BI-atomization
BI=Function('BoundIntegral')
def t1(e):
    nodes=list(e.atoms(BI)); reps={n:Symbol(f'_B{i}') for i,n in enumerate(nodes)}; inv={v:k for k,v in reps.items()}
    ee=e.xreplace(reps)
    if ee.has(eta): ee=sp.series(ee,eta,0,2).removeO()
    if ee.has(sw): ee=sp.series(ee,sw,0,2).removeO()
    return sp.expand(ee).xreplace(inv)
test=C.combine_bound_integrals(sp.expand(test)); test=t1(test)
print("\nresidual after binding delta_rho_4D_bulk_t := (rho_br/W_bg)theta_t + rho_br*sigma_W*G_t/W_bg^2, 1st shape order:")
print("   simplified ->", sp.simplify(test))
