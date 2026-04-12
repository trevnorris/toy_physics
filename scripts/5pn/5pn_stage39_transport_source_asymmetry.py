#!/usr/bin/env python3
"""
5pn_stage39_transport_source_asymmetry.py

SymPy audit for Moving-Throat PDE Stage 39.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    simplified = sp.simplify(sp.expand(expr))
    print(f"{name} = {simplified}")
    if simplified != 0:
        raise AssertionError(f"{name} is not zero")


banner("STAGE 39 — TRANSPORT ORIGIN OF THE LOWEST-LANE SOURCE ASYMMETRY")

# Symbols
Pe, x = sp.symbols("Pe x", positive=True, real=True)
pi = sp.pi

# Stationary zero-flux transport branch on x in [0,1]
p_Pe = sp.simplify(Pe * sp.exp(Pe * x) / (sp.exp(Pe) - 1))
print("p_Pe(x) =", p_Pe)
expect_zero("normalized source density", sp.simplify(sp.integrate(p_Pe, (x, 0, 1)) - 1))

# D/N lowest mode on [0,1]
chi0 = sp.sqrt(2) * sp.sin(pi * x / 2)
I_W = sp.simplify(sp.integrate(chi0, (x, 0, 1)))
E_chi = sp.simplify(sp.integrate(p_Pe * chi0, (x, 0, 1)))
Omega_Pe = sp.simplify(E_chi / I_W)

print("\nchi_0(x) =", chi0)
print("I_W      =", I_W)
print("E_Pe[chi_0] =", E_chi)
print("Omega_Pe =", Omega_Pe)

expect_zero(
    "Omega_Pe - closed form",
    Omega_Pe
    - pi * Pe * (2 * Pe * sp.exp(Pe) + pi)
    / ((4 * Pe**2 + pi**2) * (sp.exp(Pe) - 1)),
)

# Exact score-function / covariance identities
E_x = sp.simplify(sp.integrate(x * p_Pe, (x, 0, 1)))
dp_identity = sp.simplify(sp.diff(p_Pe, Pe) - (x - E_x) * p_Pe)
expect_zero("dp_Pe/dPe - (x-E[x]) p_Pe", dp_identity)

E_xchi = sp.simplify(sp.integrate(x * p_Pe * chi0, (x, 0, 1)))
cov_identity = sp.simplify(sp.diff(E_chi, Pe) - (E_xchi - E_x * E_chi))
expect_zero("dE[chi]/dPe - Cov(chi,x)", cov_identity)

expect_zero("Omega_Pe(0) - 1", sp.simplify(sp.limit(Omega_Pe, Pe, 0) - 1))
expect_zero("Omega_Pe(infty) - pi/2", sp.simplify(sp.limit(Omega_Pe, Pe, sp.oo) - pi / 2))
print("small-Pe series =", sp.series(Omega_Pe, Pe, 0, 3).removeO())
expect_zero(
    "large-Pe asymptotic coefficient",
    sp.simplify(sp.limit((pi / 2 - Omega_Pe) * Pe**2, Pe, sp.oo) - pi**3 / 8),
)

banner("STAGE 39 THEOREM LEDGER")
print("The Stage-36 exponential family is exactly the stationary zero-flux branch of")
print("  d_t sigma + d_x(-D_sigma d_x sigma + v_sigma sigma) = 0.")
print("Its asymmetry parameter is the physical Peclet number")
print("  Pe = v_sigma L / D_sigma.")
print()
print("On the constructive branch Pe >= 0, the overlap boost")
print("  Omega_Pe = pi Pe (2 Pe e^Pe + pi) / [(4 Pe^2 + pi^2)(e^Pe-1)]")
print("runs monotonically from 1 to pi/2.")
