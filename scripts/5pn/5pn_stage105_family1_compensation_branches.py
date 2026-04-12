
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
5pn_stage105_family1_compensation_branches.py

Stage 105 audit: explicit Family-1 compensated coupling branches.
"""

banner("STAGE 105 — FAMILY-1 COMPENSATION BRANCHES")

r_F1 = sp.sqrt(sp.Rational(12, 1) / sp.pi**2 * (sp.Rational(37, 20))**2 - 1)
g_minus = sp.simplify(r_F1 - sp.Rational(1, 2) * sp.sqrt(1 + r_F1**2))
g_plus = sp.simplify(r_F1 + sp.Rational(1, 2) * sp.sqrt(1 + r_F1**2))

print("mathfrak_r_F1 =", sp.N(r_F1, 20))
print("mathfrak_g_-^F1 =", sp.N(g_minus, 20))
print("mathfrak_g_+^F1 =", sp.N(g_plus, 20))

# Useful normalized traction factors from the exact traction law T_m ~ 1/mathfrak_g
tau_minus = sp.simplify(1 / g_minus)
tau_plus = sp.simplify(1 / g_plus)
print("normalized traction factor on lower branch  ~ 1/g_- =", sp.N(tau_minus, 20))
print("normalized traction factor on upper branch  ~ 1/g_+ =", sp.N(tau_plus, 20))

banner("STAGE 105 FINAL LEDGER")
print("On the explicit Family-1 geometric branch the two compensated coupling values are")
print("  mathfrak_g_-^F1 =", sp.N(g_minus, 15))
print("  mathfrak_g_+^F1 =", sp.N(g_plus, 15))
print("These are the two candidate normalized mouth-coupling branches that later positivity")
print("and mouth-source theorems must test.")
