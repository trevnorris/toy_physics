#!/usr/bin/env python3
"""
moving_throat_pde_stage10_nonconstant_axial_family_sympy_audit.py

SymPy audit for Stage 27 of the moving-throat PDE program.

Scope
-----
This script verifies the first nonconstant finite-throat wall/brane family built
on the exact N/N and D/N axial bases:

  • exact N/N orthonormality for u0 and u1,
  • exact D/N half-wave normalization for f0,
  • exact overlap constants kappa0 and kappa1,
  • exact one-parameter overlap law kappa(theta),
  • exact blind-angle and max-coupling branch values,
  • exact wall stiffness expectation K_geo(theta),
  • exact substitution into the Stage-8/9 minimal isotropic module,
  • exact recovery of Stage 9 at theta=0,
  • and the exact blind-angle no-go for the outgoing quadrupole normalization.

This remains a reduced-sector theorem, but it shows that the Stage-9 branch is
not an artifact of the constant wall profile: it is the theta=0 member of an
exact finite-throat nonconstant profile family.
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
    simplified = sp.simplify(sp.expand_trig(sp.expand(expr)))
    print(f"{name} = {simplified}")
    if simplified != 0:
        raise AssertionError(f"{name} is not zero")


# ---------------------------------------------------------------------------
# Symbols
# ---------------------------------------------------------------------------

s, L = sp.symbols("s L", positive=True, real=True)
theta = sp.symbols("theta", real=True)

lambda_B, lambda_U, lambda_W, lambda_R = sp.symbols(
    "lambda_B lambda_U lambda_W lambda_R", real=True
)
varpi = sp.symbols("varpi", positive=True, real=True)
Omega_U, Omega_W = sp.symbols("Omega_U Omega_W", positive=True, real=True)
K_eta, T_Omega, T_w = sp.symbols("K_eta T_Omega T_w", real=True)
G, c, c_s, a = sp.symbols("G c c_s a", positive=True, real=True)
mhat = sp.symbols("mhat", positive=True, real=True)


# ---------------------------------------------------------------------------
# I. Exact finite-throat basis
# ---------------------------------------------------------------------------

def finite_throat_basis() -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    banner("SECTION I — EXACT FINITE-THROAT N/N AND D/N BASIS")

    u0 = 1 / sp.sqrt(L)
    u1 = sp.sqrt(2 / L) * sp.cos(sp.pi * s / L)
    f0 = sp.sqrt(2 / L) * sp.sin(sp.pi * s / (2 * L))
    chi = sp.cos(theta) * u0 + sp.sin(theta) * u1

    print("u0(s)  =", u0)
    print("u1(s)  =", u1)
    print("f0(s)  =", f0)
    print("chi(s) =", chi)

    expect_zero("int u0^2 - 1", sp.integrate(u0**2, (s, 0, L)) - 1)
    expect_zero("int u1^2 - 1", sp.integrate(u1**2, (s, 0, L)) - 1)
    expect_zero("int u0*u1", sp.integrate(u0 * u1, (s, 0, L)))
    expect_zero("int f0^2 - 1", sp.integrate(f0**2, (s, 0, L)) - 1)

    return u0, u1, f0, chi


# ---------------------------------------------------------------------------
# II. Exact overlap law
# ---------------------------------------------------------------------------

def overlap_law() -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    banner("SECTION II — EXACT OVERLAP LAW")

    u0, u1, f0, chi = finite_throat_basis()
    kappa0 = sp.simplify(sp.integrate(u0 * f0, (s, 0, L)))
    kappa1 = sp.simplify(sp.integrate(u1 * f0, (s, 0, L)))
    kappa = sp.simplify(sp.integrate(chi * f0, (s, 0, L)))
    rho = sp.simplify(sp.sqrt(kappa0**2 + kappa1**2))

    subbanner("II.1 — Base overlap constants")
    print("kappa0 =", kappa0)
    print("kappa1 =", kappa1)
    print("kappa(theta) =", kappa)
    print("rho =", rho)

    expect_zero("kappa0 - 2*sqrt(2)/pi", kappa0 - 2 * sp.sqrt(2) / sp.pi)
    expect_zero("kappa1 + 4/(3*pi)", kappa1 + 4 / (3 * sp.pi))
    expect_zero(
        "kappa(theta) - [kappa0 cos(theta) + kappa1 sin(theta)]",
        kappa - (kappa0 * sp.cos(theta) + kappa1 * sp.sin(theta)),
    )
    expect_zero("rho - 2*sqrt(22)/(3*pi)", rho - 2 * sp.sqrt(22) / (3 * sp.pi))

    subbanner("II.2 — Blind angle and max-coupling branch")
    # Use exact trigonometric values rather than atan to keep everything algebraic.
    blind_subs = {
        sp.cos(theta): sp.sqrt(2) / sp.sqrt(11),
        sp.sin(theta): sp.Integer(3) / sp.sqrt(11),
    }
    max_subs = {
        sp.cos(theta): sp.Integer(3) / sp.sqrt(11),
        sp.sin(theta): -sp.sqrt(2) / sp.sqrt(11),
    }

    kappa_blind = sp.simplify(kappa.subs(blind_subs))
    kappa_max = sp.simplify(kappa.subs(max_subs))
    sin2_max = sp.simplify((sp.sin(theta) ** 2).subs(max_subs))

    print("kappa(blind) =", kappa_blind)
    print("kappa(max)   =", kappa_max)
    print("sin^2(theta_max) =", sin2_max)

    expect_zero("kappa(blind)", kappa_blind)
    expect_zero("kappa(max) - rho", kappa_max - rho)
    expect_zero("sin^2(theta_max) - 2/11", sin2_max - sp.Rational(2, 11))

    return u0, u1, f0, chi, kappa


# ---------------------------------------------------------------------------
# III. Wall stiffness expectation on the coherent family
# ---------------------------------------------------------------------------

def wall_stiffness() -> tuple[sp.Expr, sp.Expr]:
    banner("SECTION III — EXACT WALL STIFFNESS EXPECTATION")

    u0, u1, f0, chi, kappa = overlap_law()
    G_eta = -T_w * sp.diff(chi, s, 2) + (K_eta + 6 * T_Omega) * chi
    K_geo = sp.simplify(sp.integrate(chi * G_eta, (s, 0, L)))
    K_geo_expected = sp.simplify(K_eta + 6 * T_Omega + T_w * sp.pi**2 * sp.sin(theta) ** 2 / L**2)

    print("K_geo(theta) =", K_geo)
    expect_zero("K_geo - expected", K_geo - K_geo_expected)

    # Exact value on the max-coupling branch.
    max_subs = {
        sp.cos(theta): sp.Integer(3) / sp.sqrt(11),
        sp.sin(theta): -sp.sqrt(2) / sp.sqrt(11),
    }
    K_geo_max = sp.simplify(K_geo.subs(max_subs))
    print("K_geo(theta_max) =", K_geo_max)
    expect_zero(
        "K_geo(theta_max) - [K_eta + 6 T_Omega + 2 T_w pi^2/(11 L^2)]",
        K_geo_max - (K_eta + 6 * T_Omega + 2 * T_w * sp.pi**2 / (11 * L**2)),
    )

    return kappa, K_geo


# ---------------------------------------------------------------------------
# IV. Substitution into the Stage-8/9 minimal isotropic module
# ---------------------------------------------------------------------------

def branch_substitution() -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    banner("SECTION IV — PROFILE-DRESSED STAGE-8/9 BRANCH")

    _, _, _, _, kappa = overlap_law()
    K_geo = wall_stiffness()[1]

    C = sp.simplify(lambda_B * kappa)
    G_U = sp.simplify(lambda_U)
    G_W = sp.simplify(lambda_W * kappa)
    R = sp.simplify(lambda_R * kappa)

    Delta = sp.simplify(Omega_U**2 * Omega_W**2 - R**2)
    Q = sp.simplify(G_U**2 * Omega_W**2 + 2 * G_U * G_W * R + G_W**2 * Omega_U**2)
    P = sp.simplify(Omega_U**2 * G_W + R * G_U)

    B0 = sp.simplify(C**2 / varpi**2)
    Z0 = sp.simplify(Q / Delta)
    N0 = sp.simplify(P**2 / Delta**2)
    D0 = sp.simplify(K_geo - B0 - Z0)

    subbanner("IV.1 — Reduced branch quantities")
    print("C     =", C)
    print("G_U   =", G_U)
    print("G_W   =", G_W)
    print("R     =", R)
    print("Delta =", Delta)
    print("Q     =", Q)
    print("P     =", P)
    print("B0    =", B0)
    print("Z0    =", Z0)
    print("N0    =", N0)
    print("D0    =", D0)

    # Exact compact formulas.
    Delta_expected = sp.simplify(Omega_U**2 * Omega_W**2 - lambda_R**2 * kappa**2)
    Q_expected = sp.simplify(lambda_U**2 * Omega_W**2 + 2 * lambda_U * lambda_W * lambda_R * kappa**2 + lambda_W**2 * Omega_U**2 * kappa**2)
    P_expected = sp.simplify(kappa * (Omega_U**2 * lambda_W + lambda_R * lambda_U))
    B0_expected = sp.simplify(lambda_B**2 * kappa**2 / varpi**2)

    expect_zero("Delta - expected", Delta - Delta_expected)
    expect_zero("Q - expected", Q - Q_expected)
    expect_zero("P - expected", P - P_expected)
    expect_zero("B0 - expected", B0 - B0_expected)

    subbanner("IV.2 — Exact recovery of Stage 9 at theta = 0")
    theta0 = {sp.cos(theta): 1, sp.sin(theta): 0}
    kappa0 = sp.simplify(2 * sp.sqrt(2) / sp.pi)

    expect_zero("kappa(theta=0) - kappa0", kappa.subs(theta0) - kappa0)
    expect_zero(
        "K_geo(theta=0) - (K_eta + 6 T_Omega)",
        K_geo.subs(theta0) - (K_eta + 6 * T_Omega),
    )
    expect_zero(
        "Delta(theta=0) - [Omega_U^2 Omega_W^2 - lambda_R^2 kappa0^2]",
        Delta.subs(theta0) - (Omega_U**2 * Omega_W**2 - lambda_R**2 * kappa0**2),
    )

    return kappa, K_geo, Delta, Q, P, N0


# ---------------------------------------------------------------------------
# V. Exact normalization condition and blind-angle no-go
# ---------------------------------------------------------------------------

def normalization_and_no_go() -> None:
    banner("SECTION V — NORMALIZATION EQUATION AND BLIND-ANGLE NO-GO")

    kappa, K_geo, Delta, Q, P, N0 = branch_substitution()
    target = sp.simplify(54 * G * c_s**5 / (5 * a**5 * c**5))
    B0 = sp.simplify(lambda_B**2 * kappa**2 / varpi**2)
    K_req = sp.simplify(B0 + Q / Delta + mhat**2 * P**2 / (target * Delta**2))

    print("K_req(theta) =")
    print(K_req)

    subbanner("V.1 — Blind-angle no-go")
    blind_subs = {
        sp.cos(theta): sp.sqrt(2) / sp.sqrt(11),
        sp.sin(theta): sp.Integer(3) / sp.sqrt(11),
    }
    P_blind = sp.simplify(P.subs(blind_subs))
    N0_blind = sp.simplify(N0.subs(blind_subs))
    lhs_blind = sp.simplify((mhat**2 * N0 / (K_geo - B0 - Q / Delta)).subs(blind_subs))

    print("P(blind)   =", P_blind)
    print("N0(blind)  =", N0_blind)
    print("lhs(blind) =", lhs_blind)

    expect_zero("P(blind)", P_blind)
    expect_zero("N0(blind)", N0_blind)
    expect_zero("lhs(blind)", lhs_blind)

    print("Since the target 54*G*c_s^5/(5*a^5*c^5) is strictly positive, the blind-angle")
    print("branch is an exact no-go for the outgoing quadrupole normalization bridge.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    normalization_and_no_go()
