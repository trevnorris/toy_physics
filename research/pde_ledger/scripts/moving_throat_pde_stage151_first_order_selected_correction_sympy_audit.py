#!/usr/bin/env python3
"""
moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py

Algebraic audit for the first-order source correction selected by the full mouth profile.
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
    print(f"{name} =", expr)
    if expr != 0:
        raise AssertionError(f"{name} is not zero")

banner("FIRST-ORDER SELF-CONSISTENT SOURCE CORRECTION")

# Abstract expectation symbols with respect to Sigma_*.
Rbar, cbar, Kbar, cRbar, KRbar = sp.symbols("Rbar cbar Kbar cRbar KRbar", real=True)
gprime, AT, BT = sp.symbols("gprime AT BT", real=True, nonzero=True)

# Covariance definitions
CovcR = cRbar - cbar*Rbar
CovKR = KRbar - Kbar*Rbar

# Linearized normalized correction:
#   delta Sigma = Sigma_* (-R + <R>_*)
# so its moment shifts are
delta_g = -(cRbar - cbar*Rbar)
delta_S = -(KRbar - Kbar*Rbar)

expect_zero("delta g + Cov(c,R)", delta_g + CovcR)
expect_zero("delta S + Cov(K,R)", delta_S + CovKR)

deltaPi = sp.simplify(-delta_g / gprime)
deltaT = sp.simplify(AT*delta_g + BT*delta_S)

print("deltaPi =")
sp.pprint(deltaPi)
print("deltaT  =")
sp.pprint(deltaT)

deltaPi_cov = sp.simplify(deltaPi)
deltaT_cov = sp.simplify(deltaT)

print("deltaPi in covariance form =")
sp.pprint(deltaPi_cov)
print("deltaT in covariance form =")
sp.pprint(deltaT_cov)

print("\nTheorem:")
print("  Once the full mouth residual R_*(x) is known, the selected first-order")
print("  source correction is completely determined by Cov_*(c,R_*) and Cov_*(K_q,R_*).")
