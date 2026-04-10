
#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr) -> None:
    expr_s = sp.simplify(sp.factor(sp.together(sp.expand(expr))))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")

"""
5pn_stage110_penetration_families.py

Stage 110 audit: geometric mouth-penetration families.
"""

banner("STAGE 110 — GEOMETRIC MOUTH-PENETRATION FAMILIES")

x = sp.symbols("x", positive=True, real=True)
r_F1 = sp.sqrt(sp.Rational(12, 1) / sp.pi**2 * (sp.Rational(37, 20))**2 - 1)
g_minus = sp.simplify(r_F1 - sp.Rational(1, 2) * sp.sqrt(1 + r_F1**2))

g_slab = sp.simplify(2 * sp.sin(sp.pi * x / 2) / (sp.pi * x))
g_exp = sp.simplify(2 * (2 + sp.pi * x * sp.exp(-1 / x)) / ((4 + sp.pi**2 * x**2) * (1 - sp.exp(-1 / x))))

print("g_slab(x) =")
sp.pprint(g_slab)
print("g_exp(x) =")
sp.pprint(g_exp)

x_slab = sp.nsolve(sp.N(g_slab - g_minus, 50), 0.8)
x_exp = sp.nsolve(sp.N(g_exp - g_minus, 50), 0.66)

print("x_*^slab =", x_slab)
print("x_*^exp  =", x_exp)

banner("STAGE 110 FINAL LEDGER")
print("Simple positive penetration families reach the exact lower compensated branch at")
print("  x_*^slab ≈", x_slab)
print("  x_*^exp  ≈", x_exp)
print("So the canonical Family-1 mouth branch is reached at moderate penetration depth,")
print("not by any exotic sign-changing source.")
