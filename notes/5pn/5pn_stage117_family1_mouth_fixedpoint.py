
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
5pn_stage117_family1_mouth_fixedpoint.py

Stage 117 audit: Family-1 shell + first mixed D/N tube reduction.
"""

banner("STAGE 117 — FAMILY-1 MOUTH FIXED-POINT")

Pi = sp.symbols("Pi", positive=True, real=True)
M_s, M_q = sp.symbols("M_s M_q", real=True)
kappa_q = sp.pi / 2

S_q = sp.simplify(
    Pi * (kappa_q * sp.tanh(kappa_q) + Pi * (sp.exp(-Pi) * sp.sech(kappa_q) - 1))
    / ((1 - sp.exp(-Pi)) * (kappa_q**2 - Pi**2))
)
print("S_q(Pi) =")
sp.pprint(S_q)

r_F1 = sp.sqrt(sp.Rational(12, 1) / sp.pi**2 * (sp.Rational(37, 20))**2 - 1)
g_minus = sp.simplify(r_F1 - sp.Rational(1, 2) * sp.sqrt(1 + r_F1**2))
g_Pi = sp.simplify(2 * Pi * (2 * Pi * sp.exp(Pi) + sp.pi) / ((4 * Pi**2 + sp.pi**2) * (sp.exp(Pi) - 1)))
Pi_star = sp.nsolve(sp.N(g_Pi - g_minus, 60), 1.5)
S_q_star = sp.N(S_q.subs(Pi, Pi_star), 30)

print("Pi_* =", Pi_star)
print("S_q(Pi_*) =", S_q_star)

gain_line = sp.simplify(Pi_star - M_q * S_q_star)
print("M_s(Pi_*) =", gain_line)

banner("STAGE 117 FINAL LEDGER")
print("On the first explicit Family-1 mouth-layer branch one has")
print("  Pi = M_s + M_q S_q(Pi)")
print("with")
print("  kappa_s = 0,   kappa_q = pi/2.")
print("At the canonical Family-1 point")
print("  Pi_* ≈", Pi_star)
print("  S_q(Pi_*) ≈", S_q_star)
print("so the exact gain line is")
print("  M_s = Pi_* - M_q S_q(Pi_*).")
