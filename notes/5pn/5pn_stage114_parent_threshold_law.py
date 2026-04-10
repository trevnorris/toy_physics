
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
5pn_stage114_parent_threshold_law.py

Stage 114 audit: direct parent threshold for the canonical mouth point.
"""

banner("STAGE 114 — DIRECT PARENT THRESHOLD")

Pi = sp.symbols("Pi", positive=True, real=True)
Theta_sigma, L = sp.symbols("Theta_sigma L", positive=True, real=True)
T_m, qstar, A0p = sp.symbols("T_m q_* A0prime", real=True)

g_Pi = sp.simplify(2 * Pi * (2 * Pi * sp.exp(Pi) + sp.pi) / ((4 * Pi**2 + sp.pi**2) * (sp.exp(Pi) - 1)))
r_F1 = sp.sqrt(sp.Rational(12, 1) / sp.pi**2 * (sp.Rational(37, 20))**2 - 1)
g_minus = sp.simplify(r_F1 - sp.Rational(1, 2) * sp.sqrt(1 + r_F1**2))
Pi_star = sp.nsolve(sp.N(g_Pi - g_minus, 60), 1.5)

gprime_star = sp.N(sp.diff(g_Pi, Pi).subs(Pi, Pi_star), 30)
print("Pi_* =", Pi_star)
print("g'_*(Pi_*) =", gprime_star)

V1_star = sp.simplify(Pi_star * Theta_sigma / L)
print("Threshold V_1 =", V1_star)
print("Equivalent threshold:")
print("  T_m - q_* A0' = Pi_* Theta_sigma / L")

banner("STAGE 114 FINAL LEDGER")
print("The canonical lower compensated branch is selected when")
print("  Pi_m = (L/Theta_sigma) (T_m - q_* A0') = Pi_*")
print("with")
print("  Pi_* ≈", Pi_star)
print("and local sensitivity")
print("  g'_*(Pi_*) ≈", gprime_star)
