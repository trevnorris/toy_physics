#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

def banner(t: str):
    print("\n" + "="*88)
    print(t)
    print("="*88)

x = sp.symbols('x', positive=True, real=True)

banner("STAGE 127 — GEOMETRIC MOUTH-PENETRATION FAMILIES")

R = sp.sqrt(4107 - 100*sp.pi**2)
gminus = sp.simplify((2*R - 37*sp.sqrt(3))/(20*sp.pi))

# Uniform slab family
g_slab = sp.simplify(2*sp.sin(sp.pi*x/2)/(sp.pi*x))
print("g_slab(x) =", g_slab)

# truncated exponential family
g_exp = sp.simplify(2*(2 + sp.pi*x*sp.exp(-1/x))/((4 + sp.pi**2*x**2)*(1 - sp.exp(-1/x))))
print("g_exp(x) =", g_exp)

# Numerical roots
x0_slab = sp.nsolve(sp.Eq(g_slab, gminus), 0.8, tol=1e-30, maxsteps=200, prec=80)
x0_exp = sp.nsolve(sp.Eq(g_exp, gminus), 0.66, tol=1e-30, maxsteps=200, prec=80)

print("\ng_-^F1 =", sp.N(gminus, 30))
print("x_*^slab =", x0_slab)
print("x_*^exp  =", x0_exp)

res_slab = sp.N(g_slab.subs(x, x0_slab) - gminus, 30)
res_exp = sp.N(g_exp.subs(x, x0_exp) - gminus, 30)
print("residual slab =", res_slab)
print("residual exp  =", res_exp)

if abs(float(res_slab)) > 1e-20 or abs(float(res_exp)) > 1e-20:
    raise AssertionError("Penetration-family compensation solve failed.")

print("\nConclusion: moderate positive source penetration reaches the unique lower compensated branch.")
