
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
5pn_stage107_family1_pre_source_status.py

Stage 107 audit: pre-source-law status on the Family-1 branch.
"""

banner("STAGE 107 — FAMILY-1 PRE-SOURCE STATUS")

r_F1 = sp.sqrt(sp.Rational(12, 1) / sp.pi**2 * (sp.Rational(37, 20))**2 - 1)
g_minus = sp.simplify(r_F1 - sp.Rational(1, 2) * sp.sqrt(1 + r_F1**2))
g_plus = sp.simplify(r_F1 + sp.Rational(1, 2) * sp.sqrt(1 + r_F1**2))

print("mathfrak_r_F1 =", sp.N(r_F1, 20))
print("lower candidate  =", sp.N(g_minus, 20))
print("upper candidate  =", sp.N(g_plus, 20))
print("point-source reference =", 1)
print("self-matched derivative =", sp.N(sp.pi / 4, 20))

banner("STAGE 107 FINAL LEDGER")
print("Before any explicit mouth-source law is imposed, the Family-1 branch has")
print("one admissible-looking candidate")
print("  mathfrak_g_-^F1 ≈", sp.N(g_minus, 15))
print("and one obviously over-large candidate")
print("  mathfrak_g_+^F1 ≈", sp.N(g_plus, 15))
print("The next theorem gate is therefore not more outlet algebra, but an honest")
print("positive mouth-source law that can select between them.")
