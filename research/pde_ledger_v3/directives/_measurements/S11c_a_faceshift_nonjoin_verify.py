import sys; sys.path.insert(0, ".")
import S11c_a_cross_engine_comparator as C
import sympy as sp

py = C.load_py(C.DEFAULT_PY); wl = C.load_wl(C.DEFAULT_WL)
cases = C.extract_family("FACE_SHIFT", py.get("FACE_SHIFT"), wl.get("FACE_SHIFT"))
print("PY cases:", len(cases.py), " WL cases:", len(cases.wl))
print("extraction_notes:", cases.extraction_notes[:5])

def axes(k):
    # Key is a tuple of (axis, value) pairs or similar; print raw
    return dict(k) if isinstance(k, (tuple,list)) and all(isinstance(x,tuple) and len(x)==2 for x in k) else k

# show a couple PY keys and WL keys
print("\n--- sample PY keys ---")
for pc in cases.py[:3]: print(" ", pc.key)
print("--- sample WL keys ---")
for wc in cases.wl[:3]: print(" ", wc.key)

# find PY density cases (BULK_DENSITY) grouped by (branch,face,dof), collect per representative
def keyget(k, axis):
    d = axes(k)
    return d.get(axis) if isinstance(d, dict) else None

py_density = [pc for pc in cases.py if keyget(pc.key,"FIELD")=="BULK_DENSITY"]
print("\nPY BULK_DENSITY cases:", len(py_density))
# group by (branch,face,dof), compare representatives
from collections import defaultdict
groups = defaultdict(dict)
for pc in py_density:
    d = axes(pc.key)
    g = (d.get("BRANCH"), d.get("FACE"), d.get("DOF"))
    groups[g][d.get("DENSITY")] = pc.value
identical = True
for g, reps in groups.items():
    vals = list(reps.values())
    if len(vals) >= 2:
        same = all(sp.simplify(v - vals[0])==0 for v in vals[1:])
        identical = identical and same
        if not same:
            print("  DIFFER at", g, list(reps.keys()))
print("PY reps identical across DENSITY axis for every (branch,face,dof)?", identical)

# print one literal PY density operand and the matching WL density operand
g0 = ("LAB_HELD","MINUS","DELTA_W")
pc = None
for c in py_density:
    d=axes(c.key)
    if (d.get("BRANCH"),d.get("FACE"),d.get("DOF"))==g0:
        pc=c; break
print("\n=== PY density operand (LAB_HELD,MINUS,DELTA_W,BULK_DENSITY) ===")
print("key:", pc.key)
print("value:", pc.value)

wl_density = [wc for wc in cases.wl if keyget(wc.key,"FIELD")=="BULK_DENSITY"]
print("\nWL BULK_DENSITY cases:", len(wl_density))
wc=None
for c in wl_density:
    d=axes(c.key)
    if d.get("BRANCH")=="LAB_HELD" and d.get("FACE")=="MINUS" and d.get("DOF")=="DELTA_W":
        wc=c;break
if wc:
    print("=== WL density operand (LAB_HELD,MINUS,DELTA_W,BULK_DENSITY) ===")
    print("key:", wc.key)
    print("value:", wc.value)
