#!/usr/bin/env python3
"""
moving_throat_pde_stage8_minimal_isotropic_normalization_sympy_audit.py

SymPy audit for Stage 8 of the moving-throat PDE program.

Scope
-----
This script verifies the minimal isotropic single-support/single-port closure that
follows after the Stage-7 angular theorem:

  • exact zero-frequency coefficients B0, Z0, N0,
  • exact closed formula for P0 = N0 / D0,
  • exact target equation on the minimal isotropic branch,
  • stability-domain positivity structure,
  • and exact monotonic derivatives with respect to bare stiffness K and
    support-softening X = C^2/varpi^2.

This is still a reduced-sector theorem, but it turns the remaining normalization
problem into one explicit scalar equation in radial/axial overlap amplitudes.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    simplified = sp.simplify(sp.expand(expr))
    print(f"{name} = {simplified}")
    if simplified != 0:
        raise AssertionError(f"{name} is not zero")


# ---------------------------------------------------------------------------
# Symbols
# ---------------------------------------------------------------------------

K, varpi = sp.symbols("K varpi", positive=True, real=True)
C = sp.symbols("C", real=True)
OmegaU, OmegaW = sp.symbols("OmegaU OmegaW", positive=True, real=True)
GU, GW, R = sp.symbols("GU GW R", real=True)
G, c, c_s, a = sp.symbols("G c c_s a", positive=True, real=True)
mhat = sp.symbols("mhat", positive=True, real=True)

X = sp.symbols("X", nonnegative=True, real=True)


# ---------------------------------------------------------------------------
# I. Minimal isotropic zero-frequency coefficients
# ---------------------------------------------------------------------------

def zero_frequency_coefficients() -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    banner("SECTION I — MINIMAL ISOTROPIC ZERO-FREQUENCY COEFFICIENTS")

    Delta = sp.simplify(OmegaU**2 * OmegaW**2 - R**2)
    Q = sp.simplify(GU**2 * OmegaW**2 + 2 * GU * GW * R + GW**2 * OmegaU**2)
    P = sp.simplify(OmegaU**2 * GW + R * GU)

    B0 = sp.simplify(C**2 / varpi**2)
    Z0 = sp.simplify(Q / Delta)
    N0 = sp.simplify(P**2 / Delta**2)
    D0 = sp.simplify(K - B0 - Z0)

    print("Delta =", Delta)
    print("Q     =", Q)
    print("P     =", P)
    print("B0    =", B0)
    print("Z0    =", Z0)
    print("N0    =", N0)
    print("D0    =", D0)

    return Delta, Q, P, B0, Z0, N0, D0


# ---------------------------------------------------------------------------
# II. Exact normalization formula
# ---------------------------------------------------------------------------

def normalization_formula() -> None:
    banner("SECTION II — EXACT MINIMAL ISOTROPIC NORMALIZATION FORMULA")

    Delta, Q, P, B0, Z0, N0, D0 = zero_frequency_coefficients()
    P0 = sp.simplify(N0 / D0)
    P0_compact = sp.simplify(P**2 / (Delta * (K * Delta - Delta * C**2 / varpi**2 - Q)))

    subbanner("II.1 — P0 in raw and compact form")
    print("P0 raw    =", P0)
    print("P0 compact=", P0_compact)
    expect_zero("P0 - P0_compact", P0 - P0_compact)

    target = sp.simplify(54 * G * c_s**5 / (5 * a**5 * c**5))
    equation_residual = sp.simplify(mhat**2 * P0_compact - target)

    subbanner("II.2 — Exact target equation")
    print("Target residual =")
    sp.pprint(equation_residual)


# ---------------------------------------------------------------------------
# III. Stability and positivity structure
# ---------------------------------------------------------------------------

def stability_and_positivity() -> None:
    banner("SECTION III — STABILITY AND POSITIVITY STRUCTURE")

    Delta, Q, P, B0, Z0, N0, D0 = zero_frequency_coefficients()
    compact_denom = sp.simplify(Delta * D0)
    expect_zero("Delta*D0 - (K*Delta - Delta*C^2/varpi^2 - Q)", compact_denom - (K * Delta - Delta * C**2 / varpi**2 - Q))
    expect_zero("N0 - P^2/Delta^2", N0 - P**2 / Delta**2)

    print("If Delta > 0 and D0 > 0, then P0 > 0 whenever P != 0.")


# ---------------------------------------------------------------------------
# IV. Exact monotonic derivatives on the stable branch
# ---------------------------------------------------------------------------

def monotonic_derivatives() -> None:
    banner("SECTION IV — EXACT MONOTONIC DERIVATIVES")

    Delta, Q, P, B0, Z0, N0, D0 = zero_frequency_coefficients()
    P0 = sp.simplify(N0 / (K - X - Q / Delta))

    dP0_dK = sp.simplify(sp.diff(P0, K))
    dP0_dX = sp.simplify(sp.diff(P0, X))

    subbanner("IV.1 — Derivatives with respect to K and X = C^2/varpi^2")
    print("dP0/dK =", dP0_dK)
    print("dP0/dX =", dP0_dX)

    expect_zero("dP0/dK + N0/(K - X - Q/Delta)^2", dP0_dK + N0 / (K - X - Q / Delta)**2)
    expect_zero("dP0/dX - N0/(K - X - Q/Delta)^2", dP0_dX - N0 / (K - X - Q / Delta)**2)
    expect_zero("dP0/dX + dP0/dK", dP0_dX + dP0_dK)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    normalization_formula()
    stability_and_positivity()
    monotonic_derivatives()
