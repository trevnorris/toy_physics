import sys; sys.path.insert(0,".")
import S11c_a_cross_engine_comparator as C
import sympy as sp

py = C.load_py(C.DEFAULT_PY); wl = C.load_wl(C.DEFAULT_WL)
cs = C.extract_family("FACE_SHIFT", py.get("FACE_SHIFT"), wl.get("FACE_SHIFT"))
print("PY cases:", len(cs.py), " WL cases:", len(cs.wl))
def ax(k): return dict(k)
# does WL now carry a DENSITY axis?
wl_axes = set()
for w in cs.wl[:5]: wl_axes |= set(ax(w.key).keys())
print("WL key axes now:", sorted(wl_axes))
# sample density operand both engines (LAB_HELD, RHO4_CONSTANT, MINUS, DELTA_W)
def find(lst):
    for c in lst:
        d=ax(c.key)
        if d.get("BRANCH")=="LAB_HELD" and d.get("DENSITY")=="RHO4_CONSTANT" and d.get("FACE")=="MINUS" and d.get("DOF")=="DELTA_W" and d.get("FIELD")=="BULK_DENSITY":
            return c
    return None
pc, wc = find(cs.py), find(cs.wl)
print("\nPY density operand:", pc.value if pc else None)
print("WL density operand:", wc.value if wc else None)
print("WL still contains rhoBulkBackground?:", "rhoBulkBackground" in str(wc.value) if wc else "n/a")

# residual across all density cases under the perturbation bridge
def resid(pc, wc, face):
    A, B = pc.value, wc.value
    subs={}
    for f in B.atoms(sp.Function):
        if f.func.__name__=="rhoBulkPerturbation":
            subs[f]=sp.Symbol(f"delta_rho_4D_face_{face.lower()}")
    for s in B.atoms(sp.Symbol):
        if s.name=="rhoBulkPerturbationXJETXdw": subs[s]=sp.Symbol(f"delta_rho_4D_face_{face.lower()}_dw")
        if s.name=="e_W": subs[s]=sp.Symbol("eta_bg")*sp.Symbol("w1_profile")
    Bs = B.xreplace(subs)
    try: return sp.simplify(sp.expand(Bs - A))
    except Exception as e: return f"ERR:{e}"

from collections import defaultdict
joined=0; nonzero=[]
pyd={ (ax(c.key).get("BRANCH"),ax(c.key).get("DENSITY"),ax(c.key).get("FACE"),ax(c.key).get("DOF")):c for c in cs.py if ax(c.key).get("FIELD")=="BULK_DENSITY"}
wld={ (ax(c.key).get("BRANCH"),ax(c.key).get("DENSITY"),ax(c.key).get("FACE"),ax(c.key).get("DOF")):c for c in cs.wl if ax(c.key).get("FIELD")=="BULK_DENSITY"}
for k in pyd:
    if k in wld:
        joined+=1
        r = resid(pyd[k], wld[k], k[2])
        if r != 0: nonzero.append((k, r))
print(f"\nDENSITY cases joined by full key: {joined}/{len(pyd)}")
print("nonzero density residuals:", len(nonzero))
for k,r in nonzero[:6]: print("  ", k, "->", r)
