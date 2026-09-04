#!/usr/bin/env python3
"""F9c routed-pair severity scan over the frozen base S11c_b_exports.py.

Classifies every (prefixed, bare) F9c pair as identity-resolvable (distinct srepr)
vs genuinely ambiguous (same srepr under different keys), and reports how many are
referenced by the critical-path roots c1/c2 bind. Evidence for design §D3.
Run: python3 _measurements/f9c_pair_scan.py
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scripts"))
import sympy as sp
from sympy import srepr
from sympy.core.function import AppliedUndef
import S11c_b_exports as m

L = m.LEDGER
PREFIXES = ("s11_", "s11b_", "s11c_a_", "s11c_b_", "s11c_c_", "s11c_c1_")


def pairs():
    out = []
    for k in L:
        for p in PREFIXES:
            if k.startswith(p) and k[len(p):] in L:
                out.append((k, k[len(p):])); break
    return out


def sr(k):
    v = L[k]["value"]
    try: return srepr(v)
    except Exception: return repr(v)


ps = pairs()
same = [p for p in ps if sr(p[0]) == sr(p[1])]
diff = [p for p in ps if sr(p[0]) != sr(p[1])]
ROOTS = ["mu_theta_operator", "slab_operator", "coupling_kernel", "face_normal",
         "conormal_deriv", "face_measure_shape_deriv", "face_velocity", "relative_flux",
         "kinematic_balance", "traction", "face_shift", "closure_shape_deriv"]
refnames = set()
for r in ROOTS:
    v = L[r]["value"]
    if isinstance(v, sp.Basic):
        refnames |= {s.name for s in v.free_symbols}
        refnames |= {f.func.__name__ for f in v.atoms(AppliedUndef)}
crit = [(k, b) for (k, b) in ps if b in refnames]
crit_same = [(k, b) for (k, b) in crit if sr(k) == sr(b)]

print(f"total F9c routed pairs (prefixed + bare both present): {len(ps)}")
print(f"  identity-resolvable (distinct srepr):        {len(diff)}")
print(f"  genuinely ambiguous (same srepr):            {len(same)}")
print(f"critical-path roots reference {len(crit)} pair(s); genuinely-ambiguous on the critical path: {len(crit_same)}")
for k, b in crit:
    tag = "GENUINE-AMBIG" if sr(k) == sr(b) else "resolvable"
    print(f"  [{tag}] {b}: bare={L[b]['step']} vs {k}={L[k]['step']}")
