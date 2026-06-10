#!/usr/bin/env python3
"""
Moving-Throat PDE — Stage 56 SymPy audit.

Checks:
1. stationary zero-flux drift-diffusion branch gives the normalized exponential source family,
2. exact D/N lowest-mode overlap reproduces the Stage-053 boost formula with alpha = Pe,
3. derivative of the overlap boost equals the exact covariance expression,
4. small- and large-Pe asymptotics are correct.
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


banner("STAGE 56 — TRANSPORT ORIGIN OF THE SOURCE-SHAPE ASYMMETRY")

s, L = sp.symbols("s L", positive=True, real=True)
Dsig, vsig = sp.symbols("D_sigma v_sigma", positive=True, real=True)
Pe = sp.symbols("Pe", positive=True, real=True)

# Zero-flux stationary transport law.
sigma = sp.exp(vsig * s / Dsig)
J = -Dsig * sp.diff(sigma, s) + vsig * sigma
expect_zero("zero-flux transport residual", J)

# Normalized branch with total source strength L.
sigma_Pe = sp.simplify(Pe * sp.exp(Pe * s / L) / (sp.exp(Pe) - 1))
norm = sp.simplify(sp.integrate(sigma_Pe, (s, 0, L)))
expect_zero("normalization int sigma_Pe ds - L", norm - L)

# D/N lowest mode and overlap.
chi0 = sp.sqrt(2 / L) * sp.sin(sp.pi * s / (2 * L))
I_W = sp.simplify(sp.integrate(chi0, (s, 0, L)))
print("I_W =", I_W)

I_Pe = sp.simplify(sp.integrate(sigma_Pe * chi0, (s, 0, L)))
Omega_Pe = sp.simplify(I_Pe / I_W)
print("I_Pe =", I_Pe)
print("Omega_Pe =", Omega_Pe)

Omega_expected = sp.simplify(
    sp.pi * Pe * (2 * Pe * sp.exp(Pe) + sp.pi)
    / ((4 * Pe**2 + sp.pi**2) * (sp.exp(Pe) - 1))
)
expect_zero("Omega_Pe - expected formula", Omega_Pe - Omega_expected)

# Endpoint limits.
Omega0 = sp.simplify(sp.limit(Omega_Pe, Pe, 0))
OmegaInf = sp.simplify(sp.limit(Omega_Pe, Pe, sp.oo))
print("Omega_Pe(0) =", Omega0)
print("lim Pe->+inf Omega_Pe =", OmegaInf)
if Omega0 != 1:
    raise AssertionError("Twin baseline was not recovered.")
if sp.simplify(OmegaInf - sp.pi / 2) != 0:
    raise AssertionError("Upper finite-throat overlap limit is wrong.")

# Covariance identity dOmega/dPe = Cov(chi0,s)/I_W.
p_Pe = sp.simplify(sigma_Pe / L)
Echi = sp.simplify(sp.integrate(p_Pe * chi0, (s, 0, L)))
Es = sp.simplify(sp.integrate(p_Pe * s, (s, 0, L)))
Echis = sp.simplify(sp.integrate(p_Pe * chi0 * s, (s, 0, L)))
Cov = sp.simplify(Echis - Echi * Es)
expect_zero("dOmega/dPe - Cov/I_W", sp.diff(Omega_Pe, Pe) - Cov / I_W)

# Small- and large-bias asymptotics.
Omega_small = sp.series(Omega_Pe, Pe, 0, 3).removeO()
print("Omega_Pe small-Pe series =", sp.expand(Omega_small))
small_expected = sp.expand(1 + ((4 - sp.pi) / (2 * sp.pi)) * Pe)
expect_zero(
    "small-Pe linear coefficient",
    sp.expand(Omega_small).coeff(Pe, 1) - sp.expand(small_expected).coeff(Pe, 1),
)

Omega_large = sp.series(Omega_Pe, Pe, sp.oo, 3).removeO()
print("Omega_Pe large-Pe series =", sp.expand(Omega_large))
expect_zero(
    "large-Pe asymptotic through O(Pe^-2)",
    sp.expand(Omega_large - (sp.pi / 2 - sp.pi**3 / (8 * Pe**2))),
)

print("\nStage 56 audit passed.")
