
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
5pn_stage115_mouth_boundary_layer_status.py

Stage 115 audit: mouth boundary-layer status after explicit source-law extraction.
"""

banner("STAGE 115 — MOUTH BOUNDARY-LAYER STATUS")

Pi = sp.symbols("Pi", positive=True, real=True)
g_Pi = sp.simplify(2 * Pi * (2 * Pi * sp.exp(Pi) + sp.pi) / ((4 * Pi**2 + sp.pi**2) * (sp.exp(Pi) - 1)))
r_F1 = sp.sqrt(sp.Rational(12, 1) / sp.pi**2 * (sp.Rational(37, 20))**2 - 1)
g_minus = sp.simplify(r_F1 - sp.Rational(1, 2) * sp.sqrt(1 + r_F1**2))
Pi_star = sp.nsolve(sp.N(g_Pi - g_minus, 60), 1.5)

print("Pi_* =", Pi_star)
print("g(Pi_*) =", sp.N(g_Pi.subs(Pi, Pi_star), 20))
print("g_-^F1  =", sp.N(g_minus, 20))

banner("STAGE 115 FINAL LEDGER")
print("The mouth-source problem is no longer profile selection or branch sign.")
print("It has collapsed to one explicit parent threshold:")
print("  Pi_m = Pi_* ≈", Pi_star)
print("or equivalently")
print("  T_m - q_* A0' = Pi_* Theta_sigma / L.")
