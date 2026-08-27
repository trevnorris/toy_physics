import sys; sys.path.insert(0,".")
import S11c_a_cross_engine_comparator as C
import sympy as sp

py = C.load_py(C.DEFAULT_PY); wl = C.load_wl(C.DEFAULT_WL)
cases = C.extract_family("FACE_SHIFT", py.get("FACE_SHIFT"), wl.get("FACE_SHIFT"))
def ax(k): return dict(k)

def get(lst, want, has_density):
    for c in lst:
        d=ax(c.key)
        ok = d.get("BRANCH")==want[0] and d.get("FACE")==want[1] and d.get("DOF")==want[2] and d.get("FIELD")=="BULK_DENSITY"
        if has_density: ok = ok and d.get("DENSITY")=="RHO4_CONSTANT"
        if ok: return c
    return None

results=[]
for face,fref in (("MINUS","-W_0/2"),("PLUS","W_0/2")):
  for dof in ("DELTA_W","ZETA_C"):
    for branch in ("LAB_HELD","MATERIAL_ADVECTED"):
        pc = get(cases.py,(branch,face,dof),True)
        wc = get(cases.wl,(branch,face,dof),False)
        if pc is None or wc is None:
            results.append((branch,face,dof,"MISSING")); continue
        A = pc.value   # PY
        B = wc.value   # WL (transliterated)
        # symbols
        syms = {s.name: s for s in B.free_symbols}
        def S(n): return syms.get(n, sp.Symbol(n))
        # ---- field-name bridge (WL -> PY), perturbation side only ----
        subs = {}
        # perturbation face value: rhoBulkPerturbation(x1,x2,x3,(fref,time)) -> delta_rho_4D_face_{face}
        for f in B.atoms(sp.Function):
            nm = f.func.__name__
            if nm=="rhoBulkPerturbation":
                subs[f] = sp.Symbol(f"delta_rho_4D_face_{face.lower()}")
        # jets and background derivs (symbols)
        pert_dw = S("rhoBulkPerturbationXJETXdw")
        subs[pert_dw] = sp.Symbol(f"delta_rho_4D_face_{face.lower()}_dw")
        # e_W -> eta_bg*w1_profile   (spec s45: e_W = deltaWidth/W_0, PY deltaWidth=W_0*eta_bg*w1_profile)
        subs[S("e_W")] = sp.Symbol("eta_bg")*sp.Symbol("w1_profile")
        # GROUND the density background per §2b: rho^0 depends on in-plane anchor, not w => d_w rho^0 = 0
        subs[S("rhoBulkBackgroundXJETXdw")] = 0
        subs[S("rhoBulkBackgroundXJETXdwXdw")] = 0
        B_grounded = B.xreplace(subs)
        resid = sp.simplify(sp.expand(B_grounded - A))
        # also the UN-grounded residual (what the WL-only term is)
        subs_keep = dict(subs); subs_keep.pop(S("rhoBulkBackgroundXJETXdw")); subs_keep.pop(S("rhoBulkBackgroundXJETXdwXdw"))
        resid_raw = sp.expand(B.xreplace(subs_keep) - A)
        results.append((branch,face,dof, "grounded_resid="+str(resid), "wl_only_term="+str(resid_raw)))

for r in results:
    print(r)
