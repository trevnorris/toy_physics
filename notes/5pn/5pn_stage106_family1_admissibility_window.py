
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
5pn_stage106_family1_admissibility_window.py

Stage 106 audit: basic admissibility window of the Family-1 compensated branches.
"""

banner("STAGE 106 — FAMILY-1 ADMISSIBILITY WINDOW")

r_F1 = sp.sqrt(sp.Rational(12, 1) / sp.pi**2 * (sp.Rational(37, 20))**2 - 1)
g_minus = sp.simplify(r_F1 - sp.Rational(1, 2) * sp.sqrt(1 + r_F1**2))
g_plus = sp.simplify(r_F1 + sp.Rational(1, 2) * sp.sqrt(1 + r_F1**2))

g_point = sp.Integer(1)
g_match = sp.pi / 4
g_uniform = 2 / sp.pi

print("2/pi      =", sp.N(g_uniform, 20))
print("g_-^F1    =", sp.N(g_minus, 20))
print("pi/4      =", sp.N(g_match, 20))
print("point src =", g_point)
print("g_+^F1    =", sp.N(g_plus, 20))

banner("STAGE 106 FINAL LEDGER")
print("The explicit Family-1 branch already sits in the window")
print("  2/pi < mathfrak_g_-^F1 < pi/4 < 1 < mathfrak_g_+^F1.")
print("So the lower compensated branch is the only one compatible with any later")
print("positive-source theorem, and it already lies close to the natural self-matched")
print("derivative profile value pi/4.")
