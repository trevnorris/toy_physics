
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
5pn_stage109_positive_source_families.py

Stage 109 audit: explicit positive source families and the Family-1 compensation point.
"""

banner("STAGE 109 — EXPLICIT POSITIVE SOURCE FAMILIES")

g_minus = sp.sqrt(sp.Rational(12, 1) / sp.pi**2 * (sp.Rational(37, 20))**2 - 1) - sp.Rational(1, 2) * sp.sqrt(1 + (sp.Rational(12, 1) / sp.pi**2 * (sp.Rational(37, 20))**2 - 1))

g_match = sp.pi / 4
print("mathfrak_g_match = pi/4 =", sp.N(g_match, 20))
print("mathfrak_g_-^F1 =", sp.N(g_minus, 20))
print("Delta g_match =", sp.N(g_match - g_minus, 20))

traction_ratio = sp.simplify(g_match / g_minus)
print("T_m^{(-)} / T_m^{match} = (pi/4)/g_- =", sp.N(traction_ratio, 20))

xi = sp.symbols("xi", real=True)
g_xi = sp.simplify((1 - xi) * sp.pi / 4 + xi * 2 / sp.pi)
print("mathfrak_g_xi =", g_xi)
xi_star = sp.simplify(sp.solve(sp.Eq(g_xi, g_minus), xi)[0])
print("xi_* =", sp.simplify(xi_star))
print("xi_* (numeric) =", sp.N(xi_star, 20))

banner("STAGE 109 FINAL LEDGER")
print("The self-matched derivative source has mathfrak_g = pi/4, already close to the")
print("Family-1 lower compensated branch.")
print("A simple positive convex family")
print("  sigma_xi = (1-xi) k cos(kz) + xi/L")
print("hits the exact lower branch at")
print("  xi_* ≈", sp.N(xi_star, 15))
print("so only a modest positive broadening is needed.")
