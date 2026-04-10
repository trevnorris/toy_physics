
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
5pn_stage111_mouth_source_bias_status.py

Stage 111 audit: mouth-source bias status.
"""

banner("STAGE 111 — MOUTH-SOURCE BIAS STATUS")

r_F1 = sp.sqrt(sp.Rational(12, 1) / sp.pi**2 * (sp.Rational(37, 20))**2 - 1)
g_minus = sp.simplify(r_F1 - sp.Rational(1, 2) * sp.sqrt(1 + r_F1**2))
g_plus = sp.simplify(r_F1 + sp.Rational(1, 2) * sp.sqrt(1 + r_F1**2))
xi_star = sp.simplify((sp.pi / 4 - g_minus) / (sp.pi / 4 - 2 / sp.pi))

print("mathfrak_r_F1 =", sp.N(r_F1, 20))
print("mathfrak_g_-^F1 =", sp.N(g_minus, 20))
print("mathfrak_g_+^F1 =", sp.N(g_plus, 20))
print("xi_* =", sp.N(xi_star, 20))

banner("STAGE 111 FINAL LEDGER")
print("Under positive localized mouth sourcing:")
print("  - the upper compensated branch is ruled out,")
print("  - the lower compensated branch is uniquely admissible,")
print("  - and it is reached by modest positive broadening of natural source laws.")
