
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
5pn_stage118_outlet_consistent_mouth_closure.py

Stage 118 audit: outlet-consistent mouth closure.
"""

banner("STAGE 118 — OUTLET-CONSISTENT MOUTH CLOSURE")

Pi = sp.symbols("Pi", positive=True, real=True)
Sigma_m = sp.symbols("Sigma_m", positive=True, real=True)

S_q = sp.simplify(
    Pi * ((sp.pi / 2) * sp.tanh(sp.pi / 2) + Pi * (sp.exp(-Pi) * sp.sech(sp.pi / 2) - 1))
    / ((1 - sp.exp(-Pi)) * (sp.pi**2 / 4 - Pi**2))
)

fixed_point = sp.simplify(Sigma_m * (4 - S_q))
print("Pi = Sigma_m [4 - S_q(Pi)]")
sp.pprint(fixed_point)

r_F1 = sp.sqrt(sp.Rational(12, 1) / sp.pi**2 * (sp.Rational(37, 20))**2 - 1)
g_minus = sp.simplify(r_F1 - sp.Rational(1, 2) * sp.sqrt(1 + r_F1**2))
g_Pi = sp.simplify(2 * Pi * (2 * Pi * sp.exp(Pi) + sp.pi) / ((4 * Pi**2 + sp.pi**2) * (sp.exp(Pi) - 1)))
Pi_star = sp.nsolve(sp.N(g_Pi - g_minus, 60), 1.5)
S_q_star = sp.N(S_q.subs(Pi, Pi_star), 30)
Sigma_star = sp.N(Pi_star / (4 - S_q_star), 30)
M_s_star = sp.N(4 * Sigma_star, 30)
M_q_star = sp.N(-Sigma_star, 30)

print("Pi_* =", Pi_star)
print("Sigma_m^* =", Sigma_star)
print("M_s^* =", M_s_star)
print("M_q^* =", M_q_star)

banner("STAGE 118 FINAL LEDGER")
print("Imposing the outlet-consistent shell:mixed ratio 4:-1 collapses the mouth problem to")
print("  Pi = Sigma_m [4 - S_q(Pi)].")
print("At the canonical Family-1 point this selects")
print("  Sigma_m^* ≈", Sigma_star)
print("  M_s^* ≈", M_s_star)
print("  M_q^* ≈", M_q_star)
print("so the canonical branch is realized by a moderate one-parameter mouth gain.")
