
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
5pn_stage116_coupled_mouth_fixedpoint.py

Stage 116 audit: full coupled mouth-layer fixed-point law.
"""

banner("STAGE 116 — FULL COUPLED MOUTH-LAYER FIXED-POINT LAW")

Pi, kappa, x, G = sp.symbols("Pi kappa x G", positive=True, real=True)

Sigma = sp.simplify(Pi * sp.exp(-Pi * x) / (1 - sp.exp(-Pi)))
C = sp.simplify(G * Pi / ((1 - sp.exp(-Pi)) * (kappa**2 - Pi**2)))
A = sp.simplify(C * (kappa * sp.sinh(kappa) + Pi * sp.exp(-Pi)) / (kappa * sp.cosh(kappa)))
u = sp.simplify(A * sp.sinh(kappa * x) - C * sp.cosh(kappa * x) + C * sp.exp(-Pi * x))

expect_zero("ODE solve", -sp.diff(u, x, 2) + kappa**2 * u - G * Sigma)
expect_zero("Dirichlet mouth", u.subs(x, 0))
expect_zero("Neumann bottom", sp.diff(u, x).subs(x, 1))

S = sp.simplify(sp.diff(u, x).subs(x, 0) / G)
S_note = sp.simplify(
    Pi * (kappa * sp.tanh(kappa) + Pi * (sp.exp(-Pi) * sp.sech(kappa) - 1))
    / ((1 - sp.exp(-Pi)) * (kappa**2 - Pi**2))
)
expect_zero("mouth-derivative kernel", S - S_note)
print("S(Pi,kappa) =")
sp.pprint(S_note)

print("S(Pi,0) =", sp.simplify(sp.limit(S_note, kappa, 0)))

banner("STAGE 116 FINAL LEDGER")
print("The exact scalar D/N mouth-response kernel is")
print("  S(Pi,kappa) = Pi [kappa tanh(kappa) + Pi(e^{-Pi} sech(kappa) - 1)] / ((1-e^{-Pi})(kappa^2-Pi^2)).")
print("So the fully coupled two-channel mouth bias obeys")
print("  Pi = sum_alpha M_alpha S(Pi, kappa_alpha).")
print("This is the first explicit fixed-point replacement for the earlier effective slope datum.")
