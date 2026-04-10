
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
5pn_stage102_parent_balance_family.py

Stage 102 audit: one-parameter parent compensation family.
"""

banner("STAGE 102 — ONE-PARAMETER PARENT COMPENSATION FAMILY")

K_s, K_q, lam = sp.symbols("K_s K_q lambda", positive=True, real=True)
g_s, g_q = sp.symbols("g_s g_q", positive=True, real=True)
a = sp.symbols("a", positive=True, real=True)

r = sp.symbols("mathfrak_r", real=True)
g = sp.symbols("mathfrak_g", real=True)

r_def = sp.simplify(lam / sp.sqrt(K_s * K_q))
g_def = sp.simplify(g_q * sp.sqrt(K_s) / (g_s * sp.sqrt(K_q)))
print("mathfrak_r =", r_def)
print("mathfrak_g =", g_def)

balance = sp.expand(g_s**2 * (K_s * K_q + lam**2) - 4 * (K_s * g_q - lam * g_s)**2)
balance_rg = sp.simplify(balance / (g_s**2 * K_s * K_q))
expect_zero("normalized balance law", balance_rg - (1 + r_def**2 - 4 * (g_def - r_def)**2))

g_branches = sp.solve(sp.Eq(1 + r**2, 4 * (g - r)**2), g)
print("mathfrak_g branches =", g_branches)

L_W = sp.symbols("L_W", positive=True, real=True)
r_c = sp.symbols("r_c", nonnegative=True, real=True)
L_W_formula = sp.simplify(sp.pi * a * sp.sqrt((1 + r_c) / 3) / 2)
L_W_r = sp.simplify(L_W_formula.subs(r_c, r**2))
print("L_W(r) =")
sp.pprint(L_W_r)

banner("STAGE 102 FINAL LEDGER")
print("The exact core-balance surface reduces to")
print("  1 + mathfrak_r^2 = 4 (mathfrak_g - mathfrak_r)^2.")
print("So the compensated canonical outlet is controlled by only two normalized parent ratios:")
print("  mathfrak_r = lambda/sqrt(K_s K_q),")
print("  mathfrak_g = g_q sqrt(K_s)/(g_s sqrt(K_q)).")
print("Its two exact branches are")
for idx, branch in enumerate(g_branches, start=1):
    print(f"  branch {idx}: mathfrak_g =", branch)
print("and the D/N mixed-tube length becomes")
print("  L_W = (pi a / 2) sqrt((1 + mathfrak_r^2)/3).")
