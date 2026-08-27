"""Finding-relevant subset: the HOLD controls (REP_INVARIANCE, UNIFORM_LIMIT) must show NO genuine nonzero
cross-engine residual (route-invariance / regression hold across engines). Classify with structural zero-test
+ collapse+linearity; a genuine post-reduction nonzero = a candidate finding."""
import sys,os,re; sys.path.insert(0,"research/pde_ledger_v3/scripts")
import sympy as sp
from sympy import Symbol,Function,Integer,Rational,Mul,Add,Pow,I,Float,Dummy
from sympy.core.function import AppliedUndef
import S11c_a_cross_engine_comparator as C
OUT=os.path.expanduser("~/.s11_build/comparator_run.out")
NS=dict(Symbol=Symbol,Function=Function,Integer=Integer,Rational=Rational,Mul=Mul,Add=Add,Pow=Pow,I=I,Float=Float,Dummy=Dummy,Tuple=sp.Tuple,oo=sp.oo,pi=sp.pi)
HEADS=("delta_j_bulk","mu_theta","delta_rho","theta","currentXPerturbation","currentWPerturbation","currentWWave")
def collapse(e): return e.xreplace({a:Symbol(a.func.__name__) for a in e.atoms(AppliedUndef) if any(a.func.__name__.startswith(h) for h in HEADS)})
def allzero(v):
    if isinstance(v,tuple): return all(allzero(x) for x in v)
    return sp.simplify(v)==0
def redzero(v):
    if isinstance(v,tuple): return all(redzero(x) for x in v)
    try: return sp.simplify(C.combine_bound_integrals(collapse(v)))==0
    except Exception: return False
def classify(raw):
    raw=raw.strip()
    if raw.startswith("TextAtom"): return "UNJOINED"
    if raw.startswith("Mismatch"): return "STRUCTURAL"
    try: v=eval(raw,{"__builtins__":{}},NS)
    except Exception: return "PARSE_ERR"
    if allzero(v): return "ZERO"
    return "ZERO_REPR" if redzero(v) else "NONZERO"
def cases(fam):
    grab=False
    pat=re.compile(r'^CASE family='+re.escape(fam)+r' ')
    for l in open(OUT):
        if pat.match(l): key=l.rstrip(); grab=True; continue
        if grab and l.startswith("A_minus_B = "): yield key,l[12:].rstrip(); grab=False
for fam in ["REP_INVARIANCE_RESIDUAL","UNIFORM_LIMIT_RESIDUAL"]:
    ctr={}; nz=[]
    for k,raw in cases(fam):
        c=classify(raw); ctr[c]=ctr.get(c,0)+1
        if c=="NONZERO": nz.append(k)
    print(f"{fam}:  "+"  ".join(f"{k}={v}" for k,v in sorted(ctr.items())))
    for k in nz[:8]: print("    GENUINE NONZERO:", k.split('key=')[1])
