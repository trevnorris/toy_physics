
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
5pn_stage113_exact_mouth_bias_map.py

Stage 113 audit: exact mouth-bias map and Family-1 compensation point.
"""

banner("STAGE 113 — EXACT MOUTH-BIAS MAP")

Pi = sp.symbols("Pi", positive=True, real=True)
g_Pi = sp.simplify(2 * Pi * (2 * Pi * sp.exp(Pi) + sp.pi) / ((4 * Pi**2 + sp.pi**2) * (sp.exp(Pi) - 1)))
print("g_Pi =")
sp.pprint(g_Pi)

print("limit Pi->0+ =", sp.simplify(sp.limit(g_Pi, Pi, 0, dir='+')))
print("limit Pi->oo =", sp.simplify(sp.limit(g_Pi, Pi, sp.oo)))

r_F1 = sp.sqrt(sp.Rational(12, 1) / sp.pi**2 * (sp.Rational(37, 20))**2 - 1)
g_minus = sp.simplify(r_F1 - sp.Rational(1, 2) * sp.sqrt(1 + r_F1**2))
Pi_star = sp.nsolve(sp.N(g_Pi - g_minus, 60), 1.5)
print("Pi_* =", Pi_star)

banner("STAGE 113 FINAL LEDGER")
print("The explicit exponential boundary-layer family has exact mouth bias")
print("  g_Pi = 2 Pi (2 Pi e^Pi + pi) / ((4 Pi^2 + pi^2)(e^Pi - 1))")
print("with range 2/pi -> 1 as Pi goes 0+ -> infinity.")
print("The unique Family-1 compensation point solving g_Pi = g_-^F1 is")
print("  Pi_* ≈", Pi_star)
