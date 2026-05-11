
#!/usr/bin/env python3
"""
moving_throat_pde_stage133_coupled_mouth_fixedpoint_sympy_audit.py

SymPy audit for the explicit coupled mouth-layer fixed-point law.

Checks
------
1. Solve the scalar D/N problem with exponential source on x in [0,1].
2. Verify the exact mouth-derivative kernel S(Pi,kappa).
3. Verify the static-shell limit S(Pi,0)=1.
4. Record the general two-channel fixed-point law after diagonalization.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 116 — SCALAR D/N RESPONSE KERNEL")

x, Pi, kappa, G = sp.symbols("x Pi kappa G", positive=True, real=True)
Sigma = Pi * sp.exp(-Pi*x) / (1 - sp.exp(-Pi))

C = G * Pi / ((1 - sp.exp(-Pi)) * (kappa**2 - Pi**2))
A = sp.simplify(C * (kappa * sp.sinh(kappa) + Pi * sp.exp(-Pi)) / (kappa * sp.cosh(kappa)))

u = sp.simplify(A * sp.sinh(kappa*x) - C * sp.cosh(kappa*x) + C * sp.exp(-Pi*x))

res = sp.simplify(-sp.diff(u, x, 2) + kappa**2 * u - G * Sigma)
bc0 = sp.simplify(u.subs(x, 0))
bc1 = sp.simplify(sp.diff(u, x).subs(x, 1))

S = sp.simplify(sp.diff(u, x).subs(x, 0) / G)
S_target = sp.simplify(
    Pi * (kappa * sp.tanh(kappa) + Pi * (sp.exp(-Pi) / sp.cosh(kappa) - 1))
    / ((1 - sp.exp(-Pi)) * (kappa**2 - Pi**2))
)

print("u(x) =")
sp.pprint(u)
expect_zero("ODE residual", res)
expect_zero("u(0)", bc0)
expect_zero("u'(1)", bc1)
expect_zero("mouth derivative kernel", S - S_target)

kk = sp.symbols("kk", positive=True, real=True)
S0 = sp.simplify(sp.limit(S_target.subs(kappa, kk), kk, 0))
print("S(Pi,0) =", S0)
if sp.simplify(S0 - 1) != 0:
    raise AssertionError("Static-shell limit failed.")

banner("STAGE 116 — GENERAL TWO-CHANNEL FIXED-POINT LAW")
print("Pi = Mplus*S(Pi,kappa_plus) + Mminus*S(Pi,kappa_minus)")
print("with")
print("  S(Pi,kappa) =")
sp.pprint(S_target)
