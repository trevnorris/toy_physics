#!/usr/bin/env python3
"""
moving_throat_pde_stage150_full_profile_residual_sympy_audit.py

Lightweight SymPy audit for the full-profile Family-1 mouth potential and its curvature residual.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.together(sp.expand(expr))
    expr = sp.simplify(expr)
    print(f"{name} =", expr)
    if expr != 0:
        raise AssertionError(f"{name} is not zero")

x, Pi, Sigma = sp.symbols("x Pi Sigma", positive=True, real=True)
k = sp.pi/2

banner("FULL-PROFILE FAMILY-1 MOUTH RESIDUAL")

# Static shell profile (exact)
Ts = (1 - sp.exp(-Pi*x)) / (Pi * (1 - sp.exp(-Pi))) - x * sp.exp(-Pi)/(1 - sp.exp(-Pi))

# Mixed D/N channel profile (exact, but keep algebra lightweight)
Cq = Pi / ((1 - sp.exp(-Pi)) * (k**2 - Pi**2))
Aq = Cq * (k * sp.sinh(k) + Pi * sp.exp(-Pi)) / (k * sp.cosh(k))
Tq = Aq * sp.sinh(k*x) - Cq * sp.cosh(k*x) + Cq * sp.exp(-Pi*x)

# Hand-derived closed form: T_q'(0) = Aq*k - Cq*Pi
# (differentiate Tq = Aq*sinh(k*x) - Cq*cosh(k*x) + Cq*exp(-Pi*x) at x=0).
# Build the slope compactly from FREE coefficient symbols so the PRINTED form is provably
# the real slope; then .subs the concrete Aq, Cq definitions for the load-bearing checks.
Aq_s, Cq_s = sp.symbols("Aq Cq")
Sq_symbolic = Aq_s*k - Cq_s*Pi
Sq = Sq_symbolic.subs({Aq_s: Aq, Cq_s: Cq})

print("S_q(Pi) =", Sq_symbolic)
print("S_q(Pi) [expanded] =")
sp.pprint(Sq)

expect_zero("T_s(0)", Ts.subs(x, 0))
expect_zero("T_q(0)", Tq.subs(x, 0))
expect_zero("T_s'(0)-1", sp.diff(Ts, x).subs(x, 0) - 1)
expect_zero("T_q'(0)-S_q", sp.diff(Tq, x).subs(x, 0) - Sq)

# Use the exact compensated residual form.
R = Sigma * (4*Ts - Tq - (4 - Sq)*x)
expect_zero("R(0)", R.subs(x, 0))
expect_zero("R'(0)", sp.diff(R, x).subs(x, 0))

# Direct second derivative at the mouth.
R2 = sp.simplify(sp.diff(R, x, 2).subs(x, 0))
target_R2 = sp.simplify(-3 * Sigma * Pi / (1 - sp.exp(-Pi)))
print("R''(0) =")
sp.pprint(R2)
expect_zero("R''(0) - target", R2 - target_R2)

print("\nTheorem:")
print("  R(0)=0, R'(0)=0, and R''(0) = -3 Sigma Pi/(1-exp(-Pi)) < 0.")
print("  So the full compensated mouth potential is tangent-matched but sublinear at the mouth.")
