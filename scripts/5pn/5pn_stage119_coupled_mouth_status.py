
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
5pn_stage119_coupled_mouth_status.py

Stage 119 audit: coupled mouth-layer status after the explicit solve.
"""

banner("STAGE 119 — COUPLED MOUTH STATUS")

Pi = sp.symbols("Pi", positive=True, real=True)
S_q = sp.simplify(
    Pi * ((sp.pi / 2) * sp.tanh(sp.pi / 2) + Pi * (sp.exp(-Pi) * sp.sech(sp.pi / 2) - 1))
    / ((1 - sp.exp(-Pi)) * (sp.pi**2 / 4 - Pi**2))
)
r_F1 = sp.sqrt(sp.Rational(12, 1) / sp.pi**2 * (sp.Rational(37, 20))**2 - 1)
g_minus = sp.simplify(r_F1 - sp.Rational(1, 2) * sp.sqrt(1 + r_F1**2))
g_Pi = sp.simplify(2 * Pi * (2 * Pi * sp.exp(Pi) + sp.pi) / ((4 * Pi**2 + sp.pi**2) * (sp.exp(Pi) - 1)))
Pi_star = sp.nsolve(sp.N(g_Pi - g_minus, 60), 1.5)
S_q_star = sp.N(S_q.subs(Pi, Pi_star), 30)
Sigma_star = sp.N(Pi_star / (4 - S_q_star), 30)
Ms_star = sp.N(4 * Sigma_star, 30)
Mq_star = sp.N(-Sigma_star, 30)
mixed_feedback = sp.N(Mq_star * S_q_star, 30)

print("Pi_* =", Pi_star)
print("S_q(Pi_*) =", S_q_star)
print("M_s^* =", Ms_star)
print("M_q^* =", Mq_star)
print("M_q^* S_q(Pi_*) =", mixed_feedback)

banner("STAGE 119 FINAL LEDGER")
print("The coupled mouth-layer problem is no longer a free gain pair.")
print("On the outlet-consistent Family-1 branch it is fixed by one moderate amplitude:")
print("  Sigma_m^* ≈", Sigma_star)
print("with a nontrivial mixed-lane feedback")
print("  M_q^* S_q(Pi_*) ≈", mixed_feedback)
print("that shifts the shell demand while keeping the total bias at Pi_*.")
