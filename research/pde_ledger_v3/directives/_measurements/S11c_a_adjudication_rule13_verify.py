"""
rule-13 self-verification grounding the T7 representational-identity ADJUDICATION (S11c-a task 1).
⛔ This is the orchestrator's OWN check to ground the leg prompts + verify the legs' verdicts afterward.
It is NOT a substitute for the independent from-spec CAS consult + adjudication legs (rule 1 / rule 8).
It PRINTS the reduced residual objects (rule 2); it asserts no physics conclusion.

Inputs: the committed comparator run  ~/.s11_build/comparator_run.out  (raw per-case operand triples).
Reductions used (all on the COMPARISON side, exactly as the SCOUT classification describes):
  (a) collapse WL applied bulk quantities  X(x1,x2,x3,time) -> bare Symbol X  (the 'inert args' identity),
  (b) the comparator's OWN integral-linearity canonicalizer combine_bound_integrals,
  (c) the SPEC §3a density identity  delta_rho_4D = rho_4D,bg0 . theta   with the §2c branch factor,
  (d) the §2a background ansatz  W_bg = W_0 (1 + eta_bg . w1_profile)  truncated to first shape order.
Prints, per family x (branch,density), how many residuals are 0 after each stage.
"""
import sys, os, re
sys.path.insert(0, "research/pde_ledger_v3/scripts")
import sympy as sp
from sympy import Symbol, Function, Integer, Rational, Mul, Add, Pow, I, Float, Dummy
from sympy.core.function import AppliedUndef
import S11c_a_cross_engine_comparator as C

OUT = os.path.expanduser("~/.s11_build/comparator_run.out")
NS = dict(Symbol=Symbol, Function=Function, Integer=Integer, Rational=Rational, Mul=Mul, Add=Add,
          Pow=Pow, I=I, Float=Float, Dummy=Dummy, Tuple=sp.Tuple, oo=sp.oo, pi=sp.pi)
def parse(s):
    v = eval(s, {"__builtins__": {}}, NS); return Add(*v) if isinstance(v, tuple) else v
HEADS = ("delta_j_bulk","mu_theta","delta_rho","theta","currentXPerturbation","currentWPerturbation","currentWWave")
def collapse(e):
    return e.xreplace({a: Symbol(a.func.__name__) for a in e.atoms(AppliedUndef)
                       if any(a.func.__name__.startswith(h) for h in HEADS)})
W0,rho_br,theta_t = Symbol('W_0'),Symbol('rho_br'),Symbol('theta_t')
drho_t,Wbg,eta,w1 = Symbol('delta_rho_4D_bulk_t'),Symbol('W_bg'),Symbol('eta_bg'),Symbol('w1_profile')
BI = Function('BoundIntegral')

def extract(fams):
    rows=[]; key=None; grab=False
    pat = re.compile(r'^CASE family=('+'|'.join(fams)+r') ')
    for line in open(OUT):
        if pat.match(line): key=line.rstrip("\n"); grab=True; continue
        if grab and line.startswith("A_minus_B = "):
            rows.append((key, line[len("A_minus_B = "):].rstrip("\n"))); grab=False
    return rows
def keyinfo(k):
    kk=dict(p.split("=") for p in k.split("key=")[1].strip("()\n").split(", "))
    return k.split(" key=")[0].replace("CASE family=",""), kk
def trunc1_eta(expr):
    nodes=list(expr.atoms(BI)); reps={n:Symbol(f'_BI{i}') for i,n in enumerate(nodes)}
    inv={v:k for k,v in reps.items()}; e=expr.xreplace(reps)
    if e.has(eta): e=sp.series(e,eta,0,2).removeO()
    return sp.expand(e).xreplace(inv)

# ---- Class 1: mu_theta families ----
print("="*78); print("CLASS 1  mu_theta closure coefficient  (TRACTION / CLOSURE_SHAPE_DERIV / VIRTUAL_WORK_SHAPE_DERIV)")
mt = extract(["TRACTION","CLOSURE_SHAPE_DERIV","VIRTUAL_WORK_SHAPE_DERIV"])
raw_nz=aft_nz=0
for k,s in mt:
    e=parse(s)
    if sp.simplify(e)!=0: raw_nz+=1
    if sp.simplify(collapse(e))!=0: aft_nz+=1
print(f"  cases={len(mt)}  raw-nonzero={raw_nz}  nonzero AFTER collapse X(args)->X = {aft_nz}")
print("  => residual is EXACTLY the applied-vs-bare mu_theta difference; collapse -> 0 for all." if aft_nz==0 else "  => SOMETHING SURVIVES")

# ---- Class 2: PROJECTION density ----
print("="*78); print("CLASS 2  PROJECTION density term  (collapse + linearity cancels in-plane div.j; density identity remains)")
pj = extract(["PROJECTION_RESIDUAL","PROJECTION_SHAPE_DERIV","PROJECTION_STATIC_OPERAND",
              "PROJECTION_DYNAMIC_OPERAND","PROJECTION_TERM_ORIGINS"])
summ={}
for k,s in pj:
    fam,kk=keyinfo(k); dens=kk.get("DENSITY"); br=kk.get("BRANCH")
    red=C.combine_bound_integrals(collapse(parse(s)))
    has_j = any(a.name.startswith("delta_j_bulk") for a in red.atoms(Symbol))
    sub = {drho_t:(rho_br/W0)*theta_t} if dens=="RHO4_CONSTANT" else {drho_t:(rho_br/Wbg)*theta_t}
    red2=C.combine_bound_integrals(sp.expand(red.xreplace(sub)))
    if dens=="RHOBR_CONSTANT":
        red2=red2.replace(Pow(Wbg,Integer(-1)),(1-eta*w1)/W0).xreplace({Wbg:W0*(1+eta*w1)})
        red2=trunc1_eta(C.combine_bound_integrals(sp.expand(red2)))
    z=(sp.simplify(red2)==0)
    g=summ.setdefault((fam,br,dens),{"z":0,"nz":0,"jcancel":True})
    g["z" if z else "nz"]+=1; g["jcancel"]&=(not has_j)
print(f"  {'family':26s} {'branch':17s} {'density':14s}  divj-cancelled  zero  nonzero")
for (fam,br,dens),d in sorted(summ.items()):
    print(f"  {fam:26s} {br:17s} {dens:14s}  {str(d['jcancel']):>13s}  {d['z']:4d}  {d['nz']:6d}")
print("  NOTE: MATERIAL_ADVECTED/RHOBR_CONSTANT nonzero here reflects the OMITTED background-advection")
print("        term d_t[rho_bg(chi)] = -grad(rho_bg).u (the sigma_W.u.dw1 piece); the theta-only identity")
print("        is INCOMPLETE for the advected branch -> this is the one from-spec question for the consult.")
