#!/usr/bin/env python3
"""
moving_throat_pde_stage9_concrete_axial_overlaps_sympy_audit.py

SymPy audit for Stage 9 of the moving-throat PDE program.

Scope
-----
This script verifies the first concrete finite-throat axial branch used after the
Stage-8 minimal isotropic normalization closure:

  • exact N/N zero mode and D/N half-wave normalization,
  • exact overlap law kappa_n between the constant mode and the D/N ladder,
  • the minimal-branch value kappa = 2*sqrt(2)/pi,
  • exact substitution of the concrete overlaps into the Stage-8 quantities
    (C, G_U, G_W, R, Delta, Q, P, B0, Z0, N0, D0),
  • and exact solution of the normalization equation for the required wall
    stiffness K_req on the branch.

This remains a reduced-sector theorem, but it converts the overlap problem from
formal integrals into one explicit branch-level algebraic normalization test.
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

s, L = sp.symbols("s L", positive=True, real=True)
n = sp.symbols("n", integer=True, nonnegative=True)

lambda_B, lambda_U, lambda_W, lambda_R = sp.symbols(
    "lambda_B lambda_U lambda_W lambda_R", real=True
)
varpi = sp.symbols("varpi", positive=True, real=True)
Omega_U, Omega_W = sp.symbols("Omega_U Omega_W", positive=True, real=True)
K = sp.symbols("K", real=True)
G, c, c_s, a = sp.symbols("G c c_s a", positive=True, real=True)
mhat = sp.symbols("mhat", positive=True, real=True)
K_eta, T_Omega = sp.symbols("K_eta T_Omega", real=True)


# ---------------------------------------------------------------------------
# I. Concrete finite-throat modes
# ---------------------------------------------------------------------------

def concrete_modes() -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    banner("SECTION I — CONCRETE FINITE-THROAT MODES")

    u0 = 1 / sp.sqrt(L)
    f_n = sp.sqrt(2 / L) * sp.sin((n + sp.Rational(1, 2)) * sp.pi * s / L)
    f0 = sp.simplify(f_n.subs(n, 0))

    print("u0(s) =", u0)
    print("f_n(s) =", f_n)
    print("f0(s) =", f0)

    expect_zero("int u0^2 - 1", sp.integrate(u0**2, (s, 0, L)) - 1)
    expect_zero("int f0^2 - 1", sp.integrate(f0**2, (s, 0, L)) - 1)

    return u0, f_n, f0


# ---------------------------------------------------------------------------
# II. Exact overlap law kappa_n
# ---------------------------------------------------------------------------

def overlap_law() -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    banner("SECTION II — EXACT OVERLAP LAW")

    u0, f_n, f0 = concrete_modes()
    kappa_n = sp.simplify(sp.integrate(u0 * f_n, (s, 0, L)))
    kappa_n_expected = sp.sqrt(2) / ((n + sp.Rational(1, 2)) * sp.pi)
    kappa = sp.simplify(kappa_n.subs(n, 0))
    kappa_expected = sp.simplify(2 * sp.sqrt(2) / sp.pi)

    subbanner("II.1 — General D/N overlap with the constant zero mode")
    print("kappa_n =", kappa_n)
    print("kappa_n_expected =", kappa_n_expected)
    expect_zero("kappa_n - expected", kappa_n - kappa_n_expected)

    subbanner("II.2 — Lowest-branch overlap")
    print("kappa =", kappa)
    expect_zero("kappa - 2*sqrt(2)/pi", kappa - kappa_expected)
    print("kappa numeric =", sp.N(kappa, 15))

    return u0, f0, kappa


# ---------------------------------------------------------------------------
# III. Concrete branch substitution into Stage-8 quantities
# ---------------------------------------------------------------------------

def branch_substitution() -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    banner("SECTION III — CONCRETE BRANCH SUBSTITUTION INTO STAGE-8 QUANTITIES")

    u0, f0, kappa = overlap_law()

    # Concrete overlaps on the chosen branch:
    I_eta_phi = sp.simplify(sp.integrate(u0 * f0, (s, 0, L)))
    I_eta_u = sp.simplify(sp.integrate(u0 * u0, (s, 0, L)))
    I_eta_w = sp.simplify(sp.integrate(u0 * f0, (s, 0, L)))
    I_u_w = sp.simplify(sp.integrate(u0 * f0, (s, 0, L)))

    subbanner("III.1 — Explicit overlap integrals on the branch")
    print("I_(eta,phi) =", I_eta_phi)
    print("I_(eta,u)   =", I_eta_u)
    print("I_(eta,w)   =", I_eta_w)
    print("I_(u,w)     =", I_u_w)

    expect_zero("I_(eta,phi) - kappa", I_eta_phi - kappa)
    expect_zero("I_(eta,u) - 1", I_eta_u - 1)
    expect_zero("I_(eta,w) - kappa", I_eta_w - kappa)
    expect_zero("I_(u,w) - kappa", I_u_w - kappa)

    C = sp.simplify(lambda_B * I_eta_phi)
    GU = sp.simplify(lambda_U * I_eta_u)
    GW = sp.simplify(lambda_W * I_eta_w)
    R = sp.simplify(lambda_R * I_u_w)

    Delta = sp.simplify(Omega_U**2 * Omega_W**2 - R**2)
    Q = sp.simplify(GU**2 * Omega_W**2 + 2 * GU * GW * R + GW**2 * Omega_U**2)
    P = sp.simplify(Omega_U**2 * GW + R * GU)

    B0 = sp.simplify(C**2 / varpi**2)
    Z0 = sp.simplify(Q / Delta)
    N0 = sp.simplify(P**2 / Delta**2)
    D0 = sp.simplify(K - B0 - Z0)
    P0 = sp.simplify(N0 / D0)

    subbanner("III.2 — Reduced coefficients")
    print("C     =", C)
    print("G_U   =", GU)
    print("G_W   =", GW)
    print("R     =", R)
    print("Delta =", Delta)
    print("Q     =", Q)
    print("P     =", P)
    print("B0    =", B0)
    print("Z0    =", Z0)
    print("N0    =", N0)
    print("D0    =", D0)
    print("P0    =", P0)

    # Check against the hand-derived compact formulas.
    Delta_expected = sp.simplify(Omega_U**2 * Omega_W**2 - kappa**2 * lambda_R**2)
    Q_expected = sp.simplify(
        lambda_U**2 * Omega_W**2
        + 2 * kappa**2 * lambda_U * lambda_W * lambda_R
        + kappa**2 * lambda_W**2 * Omega_U**2
    )
    P_expected = sp.simplify(kappa * (Omega_U**2 * lambda_W + lambda_R * lambda_U))
    B0_expected = sp.simplify(kappa**2 * lambda_B**2 / varpi**2)
    P0_expected = sp.simplify(P_expected**2 / (Delta_expected * (K * Delta_expected - Delta_expected * B0_expected - Q_expected)))

    expect_zero("Delta - Delta_expected", Delta - Delta_expected)
    expect_zero("Q - Q_expected", Q - Q_expected)
    expect_zero("P - P_expected", P - P_expected)
    expect_zero("B0 - B0_expected", B0 - B0_expected)
    expect_zero("P0 - P0_expected", P0 - P0_expected)

    return kappa, Delta, Q, P, B0, Z0, N0, D0


# ---------------------------------------------------------------------------
# IV. Exact branch-level normalization test and required wall stiffness
# ---------------------------------------------------------------------------

def normalization_test() -> None:
    banner("SECTION IV — BRANCH-LEVEL NORMALIZATION TEST")

    kappa, Delta, Q, P, B0, Z0, N0, D0 = branch_substitution()
    target = sp.simplify(54 * G * c_s**5 / (5 * a**5 * c**5))
    residual = sp.simplify(mhat**2 * N0 / D0 - target)

    subbanner("IV.1 — Target residual")
    print("Target residual =")
    sp.pprint(residual)

    # Solve exactly for the required wall stiffness.
    K_req = sp.solve(sp.Eq(residual, 0), K)[0]
    K_req_expected = sp.simplify(
        B0 + Q / Delta + mhat**2 * P**2 / (target * Delta**2)
    )

    subbanner("IV.2 — Exact required wall stiffness")
    print("K_req =")
    sp.pprint(sp.simplify(K_req))
    expect_zero("K_req - expected", K_req - K_req_expected)

    # Constant wall profile => no axial-gradient contribution.
    K_geom = sp.simplify(K_eta + 6 * T_Omega)
    print("On the constant wall branch, bare quadrupole wall stiffness is K =", K_geom)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    normalization_test()
