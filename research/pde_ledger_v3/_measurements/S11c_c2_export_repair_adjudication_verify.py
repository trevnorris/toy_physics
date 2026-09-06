"""Orchestrator rule-13 check: is publication_compact faithful, and does the
semantic guard BITE a real value change? Settles Grok(bites)/Agent(inconclusive).
Same-process (no serialization/cross-process dummy artifact)."""
import sys
sys.path.insert(0, '/var/projects/toy_physics/research/pde_ledger_v3/scripts')
import sympy as sp
from sympy.core.symbol import Str
import S11c_c2_selfenergy_fold_sympy_audit as M

# --- replicate astra's algebraic-leaf comparison (publish semantic() leaf branch) ---
def poles(v): return {p for p in v.atoms(sp.Pow) if p.exp.is_negative}
def leaf_diff(original, decoded):
    atoms = (poles(original)|poles(decoded)|
             original.atoms(sp.Function,sp.Derivative)|decoded.atoms(sp.Function,sp.Derivative))
    protected = {v: sp.Dummy() for v in atoms}
    unprotect = {v:k for k,v in protected.items()}
    def norm(value):
        repl={}
        for integral in value.atoms(sp.Integral):
            integrand = sp.expand(norm(integral.function).xreplace(protected)).xreplace(unprotect)
            repl[integral] = sp.S.Zero if integrand==0 else sp.Integral(integrand,*integral.limits)
        return value.xreplace(repl)
    expr = norm(decoded)-norm(original)
    integrals = {v: sp.Dummy() for v in expr.atoms(sp.Integral)}
    return sp.expand(expr.xreplace(integrals))

def leaves(obj, path='', out=None):
    if out is None: out=[]
    if isinstance(obj, sp.MatrixBase):
        for i in range(obj.rows):
            for j in range(obj.cols): leaves(obj[i,j], f'{path}[{i},{j}]', out)
    elif isinstance(obj,(sp.Tuple,tuple,list)):
        for i,x in enumerate(obj): leaves(x, f'{path}[{i}]', out)
    elif isinstance(obj, dict):
        for k,v in obj.items(): leaves(v, f'{path}/{k}', out)
    elif isinstance(obj, sp.Expr):
        out.append((path,obj))
    return out

# --- setup (as run() does) ---
from pathlib import Path
ROOT=M.ROOT
fold,audit=M.load_model(str(ROOT/'scripts/S11c_b_exports.py'),str(ROOT/'scripts/S11c_c1_exports.py'))
la=M.assert_lookups_equal_manifest(M.bind_inputs,fold,M.IMPORT_KEYS)
inputs=la['result']
model=M.build_case(inputs, M.ANCHORINGS[0], M.DENSITIES[0])
op=model['closed']
ls=leaves(op)
print('closed slab operator: %d algebraic leaves'%len(ls))

faithful=0; faithful_fail=[]
bite2=0; bite2_fail=[]
bite_marker=0; bite_marker_fail=[]
MARK=sp.Symbol('ORCH_CORRUPTION_MARKER')
for path,L in ls:
    C=M.publication_compact(L)
    # 1. faithfulness: compact must equal original
    d=leaf_diff(L,C)
    if d==0: faithful+=1
    else: faithful_fail.append((path,str(d)[:80]))
    if L==0: continue  # can't corrupt a zero leaf meaningfully
    # 2. bite: doubling the compact value must be caught
    d2=leaf_diff(L, sp.Integer(2)*C)
    if d2!=0: bite2+=1
    else: bite2_fail.append(path)
    # 3. bite: adding a distinct marker term must be caught (real value change, pole-preserving)
    dm=leaf_diff(L, C+MARK)
    if dm!=0: bite_marker+=1
    else: bite_marker_fail.append(path)

nonzero=[p for p,L in ls if L!=0]
print('FAITHFUL (compact==emitted):      %d/%d  fails=%s'%(faithful,len(ls),faithful_fail[:3]))
print('GUARD BITES x2 corruption:        %d/%d  slipped=%s'%(bite2,len(nonzero),bite2_fail[:3]))
print('GUARD BITES +marker corruption:   %d/%d  slipped=%s'%(bite_marker,len(nonzero),bite_marker_fail[:3]))
# 4. bite inside an Integral integrand (harder case)
for path,L in ls:
    if L.atoms(sp.Integral):
        C=M.publication_compact(L)
        one_int=next(iter(C.atoms(sp.Integral)))
        corrupt=C.xreplace({one_int: sp.Integral(one_int.function+MARK, *one_int.limits)})
        di=leaf_diff(L, corrupt)
        print('INTEGRAND corruption on %s: bites=%s'%(path, di!=0))
        break
