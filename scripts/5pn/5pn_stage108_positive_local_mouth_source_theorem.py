
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
5pn_stage108_positive_local_mouth_source_theorem.py

Stage 108 audit: positive local mouth-source theorem.
"""

banner("STAGE 108 — POSITIVE LOCAL MOUTH-SOURCE THEOREM")

L, z = sp.symbols("L z", positive=True, real=True)

g_sigma = sp.Symbol("mathfrak_g[sigma]", real=True)
print("For any positive normalized sigma(z) on [0,L],")
print("  mathfrak_g[sigma] = Integral[sigma(z) cos(pi z/(2L)), {z,0,L}]")
print("and since 0 <= cos(pi z/(2L)) <= 1, one has")
print("  0 <= mathfrak_g[sigma] <= 1.")

r_F1 = sp.sqrt(sp.Rational(12, 1) / sp.pi**2 * (sp.Rational(37, 20))**2 - 1)
g_minus = sp.simplify(r_F1 - sp.Rational(1, 2) * sp.sqrt(1 + r_F1**2))
g_plus = sp.simplify(r_F1 + sp.Rational(1, 2) * sp.sqrt(1 + r_F1**2))

print("mathfrak_g_-^F1 =", sp.N(g_minus, 20))
print("mathfrak_g_+^F1 =", sp.N(g_plus, 20))

banner("STAGE 108 FINAL LEDGER")
print("Because every positive normalized mouth source obeys 0 <= mathfrak_g <= 1,")
print("the explicit Family-1 values imply:")
print("  mathfrak_g_+^F1 > 1  -> upper branch impossible,")
print("  0 < mathfrak_g_-^F1 < 1 -> lower branch is the unique physically admissible canonical candidate.")
