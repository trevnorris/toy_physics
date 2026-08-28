import sys; sys.path.insert(0,".")
import S11c_a_cross_engine_comparator as C
import sympy as sp
py = C.load_py(C.DEFAULT_PY); wl = C.load_wl(C.DEFAULT_WL)
cs = C.extract_family("FACE_SHIFT", py.get("FACE_SHIFT"), wl.get("FACE_SHIFT"))
def ax(k): return dict(k)
W = C.BOUND_BINDER
# current fields
FIELDS = ["CURRENT_X1","CURRENT_X2","CURRENT_X3","NORMAL_CURRENT"]
wlmap = {"CURRENT_X1":"currentXPerturbation1","CURRENT_X2":"currentXPerturbation2",
         "CURRENT_X3":"currentXPerturbation3","NORMAL_CURRENT":"currentWPerturbation"}
pymap = {"CURRENT_X1":"delta_j_bulk_1","CURRENT_X2":"delta_j_bulk_2",
         "CURRENT_X3":"delta_j_bulk_3","NORMAL_CURRENT":"delta_j_bulk_4"}
def get(lst, br, den, face, dof, field):
    for c in lst:
        d=ax(c.key)
        if (d.get("BRANCH")==br and d.get("DENSITY")==den and d.get("FACE")==face
            and d.get("DOF")==dof and d.get("FIELD")==field): return c
    return None
nonzero=0; checked=0
for field in FIELDS:
  for face,fref in (("MINUS",-sp.Rational(1,2)),("PLUS",sp.Rational(1,2))):
    for dof in ("DELTA_W","ZETA_C"):
      for br in ("LAB_HELD","MATERIAL_ADVECTED"):
        pc=get(cs.py,br,"RHO4_CONSTANT",face,dof,field); wc=get(cs.wl,br,"RHO4_CONSTANT",face,dof,field)
        if pc is None or wc is None: continue
        A,B=pc.value,wc.value
        subs={}
        # WL currentXPerturbation(x1,x2,x3,(fref*W0,time))_jet -> PY delta_j_bulk_i ; keep w on jets
        for f in B.atoms(sp.Function):
            nm=f.func.__name__
            if nm==wlmap[field]: subs[f]=sp.Symbol(pymap[field])
        for s in B.atoms(sp.Symbol):
            if s.name==wlmap[field]+"XJETXdw": subs[s]=sp.Symbol(f"d_w_{pymap[field].replace('delta_j_bulk','delta_j_bulk')}_{face.lower()}") 
            if s.name=="e_W": subs[s]=sp.Symbol("eta_bg")*sp.Symbol("w1_profile")
        Bs=B.xreplace(subs)
        try: r=sp.simplify(sp.expand(Bs-A))
        except Exception as e: r=f"ERR"
        checked+=1
        if r!=0: nonzero+=1
        if r!=0 and nonzero<=4: print(f"  FLAG {field} {br} {face} {dof}: A={A}\n     B(bridged)={Bs}\n     resid={r}")
print(f"\ncurrent-field cases checked={checked}, nonzero-after-bridge={nonzero}")
