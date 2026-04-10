
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
5pn_stage120_explicit_core_to_mouth_gain_map.py

Stage 120 audit: explicit core-to-mouth gain map.
"""

banner("STAGE 120 — EXPLICIT CORE-TO-MOUTH GAIN MAP")

L, Theta_sigma = sp.symbols("L Theta_sigma", positive=True, real=True)
K_s, K_q, lam = sp.symbols("K_s K_q lambda", positive=True, real=True)
g_s, g_q = sp.symbols("g_s g_q", real=True)

rho_c = sp.simplify(g_s**2 / K_s)
sigma_c = sp.simplify((K_s * g_q - lam * g_s)**2 / (K_s * (K_s * K_q + lam**2)))

M_s = sp.simplify(L * rho_c / Theta_sigma)
M_q = sp.simplify(-L * sigma_c / Theta_sigma)

print("rho_c =", rho_c)
print("sigma_c =", sigma_c)
print("M_s =")
sp.pprint(M_s)
print("M_q =")
sp.pprint(M_q)

Pi = sp.symbols("Pi", positive=True, real=True)
S_q = sp.simplify(
    Pi * ((sp.pi / 2) * sp.tanh(sp.pi / 2) + Pi * (sp.exp(-Pi) * sp.sech(sp.pi / 2) - 1))
    / ((1 - sp.exp(-Pi)) * (sp.pi**2 / 4 - Pi**2))
)
print("Family-1 fixed-point law:")
print("  Pi = M_s + M_q S_q(Pi)")

banner("STAGE 120 FINAL LEDGER")
print("The actual coupled mouth-layer gains are derived directly from the explicit core outlet:")
print("  M_s = L g_s^2 / (K_s Theta_sigma)")
print("  M_q = -L (K_s g_q - lambda g_s)^2 / (K_s (K_s K_q + lambda^2) Theta_sigma)")
print("So the mouth fixed-point ambiguity has collapsed from an abstract gain pair")
print("to one explicit set of parent core quantities.")
