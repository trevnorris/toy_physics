#!/usr/bin/env python3
"""
moving_throat_pde_stage154_coevolving_core_mouth_sympy_audit.py

SymPy-backed audit for the exact co-evolving Family-1 core–mouth map.

Checks:
1. Exact Family-1 core ratio law R(g).
2. Exact compensation equivalence g = g_* <=> R = 1/4.
3. Exact defect-transport expansion of R around the lower compensated branch.
4. Linearized slope identity dPi in terms of dSigma0, dS, dR.
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

banner("STAGE 137 — EXACT CO-EVOLVING CORE–MOUTH MAP")

g, r, dg = sp.symbols("g r dg", real=True)
R = (g - r)**2 / (1 + r**2)
print("R(g) =", sp.simplify(R))

g_star = r - sp.sqrt(1 + r**2) / 2
expect_zero("R(g_star) - 1/4", R.subs(g, g_star) - sp.Rational(1, 4))

R_shift = sp.expand(R.subs(g, g_star + dg))
R_shift_expected = sp.expand(sp.Rational(1, 4) - dg / sp.sqrt(1 + r**2) + dg**2 / (1 + r**2))
expect_zero("exact shifted R formula", R_shift - R_shift_expected)

banner("Linearized slope identity")
Sigma0, dSigma0, Sstar, dS = sp.symbols("Sigma0 dSigma0 Sstar dS", real=True)
Rstar, dR = sp.symbols("Rstar dR", real=True)
Pi = (Sigma0 + dSigma0) * (1 - (Rstar + dR) * (Sstar + dS))
Pi_lin = sp.expand(sp.series(Pi, dSigma0, 0, 2).removeO())
# manually drop quadratic dR*dS, dSigma*dR, dSigma*dS terms
Pi_lin = sp.expand(Pi).subs({
    dSigma0*dR: 0,
    dSigma0*dS: 0,
    dR*dS: 0,
    dSigma0*dR*dS: 0,
})
Pi0 = Sigma0 * (1 - Rstar*Sstar)
dPi = sp.expand(Pi_lin - Pi0)
dPi_expected = sp.expand((1 - Rstar*Sstar)*dSigma0 - Sigma0*(Rstar*dS + Sstar*dR))
expect_zero("dPi identity", dPi - dPi_expected)

print("\nCarry-forward formulas:")
print("  R(g) = (g-r)^2/(1+r^2)")
print("  g = g_*  <=>  R = 1/4 on the lower branch")
print("  delta R = -delta g/sqrt(1+r^2) + delta g^2/(1+r^2)")
print("  delta Pi = (1-R_*S_*) delta Sigma0 - Sigma0_* (R_* delta S + S_* delta R)")
