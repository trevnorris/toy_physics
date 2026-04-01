#!/usr/bin/env python3
"""
moving_throat_pde_stage11_loaded_profile_selection_sympy_audit.py

SymPy audit for Stage 11 of the moving-throat PDE program.

Scope
-----
This script verifies the first loaded axial wall-profile selection model in the
{u0,u1} N/N basis:

  • exact bare wall stiffness matrix,
  • exact rank-1 attractive loading from the D/N overlap vector,
  • exact trace/determinant/eigenvalue formulas,
  • exact profile-angle equation tan(2 theta),
  • weak-loading and strong-loading limits,
  • exact identification of the strong-loading angle with the Stage-10
    max-coupling branch,
  • and the exact softening threshold alpha_crit.

This is the first reduced spectral problem in which the Stage-10 profile angle
becomes an output rather than a free parameter.
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
# Symbols and exact overlap constants from Stage 10
# ---------------------------------------------------------------------------

alpha = sp.symbols("alpha", positive=True, real=True)
theta = sp.symbols("theta", real=True)
K_eta, T_Omega = sp.symbols("K_eta T_Omega", real=True)
T_w, L = sp.symbols("T_w L", positive=True, real=True)

kappa0 = sp.simplify(2 * sp.sqrt(2) / sp.pi)
kappa1 = sp.simplify(-4 / (3 * sp.pi))
DeltaK = sp.simplify(T_w * sp.pi**2 / L**2)
K0 = sp.simplify(K_eta + 6 * T_Omega)
K1 = sp.simplify(K0 + DeltaK)


# ---------------------------------------------------------------------------
# I. Bare and loaded wall operators
# ---------------------------------------------------------------------------

def loaded_operator() -> tuple[sp.Matrix, sp.Expr, sp.Expr]:
    banner("SECTION I — BARE AND LOADED WALL OPERATORS")

    v = sp.Matrix([kappa0, kappa1])
    K_bare = sp.Matrix([[K0, 0], [0, K1]])
    K_eff = sp.simplify(K_bare - alpha * (v * v.T))

    print("v =")
    sp.pprint(v)
    print("K_bare =")
    sp.pprint(K_bare)
    print("K_eff =")
    sp.pprint(K_eff)

    # Exact sign data.
    subbanner("I.1 — Exact sign ledger")
    print("kappa0^2 - kappa1^2 =", sp.simplify(kappa0**2 - kappa1**2))
    print("2 kappa0 kappa1    =", sp.simplify(2 * kappa0 * kappa1))

    return K_eff, K0, K1


# ---------------------------------------------------------------------------
# II. Trace, determinant, eigenvalues
# ---------------------------------------------------------------------------

def eigenvalue_ledger() -> tuple[sp.Matrix, sp.Expr, sp.Expr, sp.Expr]:
    banner("SECTION II — TRACE, DETERMINANT, AND EXACT EIGENVALUES")

    K_eff, K0_expr, K1_expr = loaded_operator()

    tr_eff = sp.simplify(sp.trace(K_eff))
    det_eff = sp.simplify(K_eff.det())
    tr_expected = sp.simplify(K0 + K1 - alpha * (kappa0**2 + kappa1**2))
    det_expected = sp.simplify(K0 * K1 - alpha * (K1 * kappa0**2 + K0 * kappa1**2))

    print("trace(K_eff) =", tr_eff)
    print("det(K_eff)   =", det_eff)
    expect_zero("trace - expected", tr_eff - tr_expected)
    expect_zero("det - expected", det_eff - det_expected)

    disc = sp.simplify((DeltaK + alpha * (kappa0**2 - kappa1**2))**2 + 4 * alpha**2 * kappa0**2 * kappa1**2)
    lam_minus = sp.simplify((tr_expected - sp.sqrt(disc)) / 2)
    lam_plus = sp.simplify((tr_expected + sp.sqrt(disc)) / 2)

    # Check characteristic polynomial.
    x = sp.symbols("x", real=True)
    char = sp.expand((x - lam_minus) * (x - lam_plus) - (x**2 - tr_eff * x + det_eff))
    expect_zero("characteristic-factorization check", char)

    print("lambda_- =")
    sp.pprint(lam_minus)
    print("lambda_+ =")
    sp.pprint(lam_plus)

    return K_eff, lam_minus, lam_plus, det_eff


# ---------------------------------------------------------------------------
# III. Exact profile-angle equation
# ---------------------------------------------------------------------------

def angle_equation() -> tuple[sp.Expr, sp.Expr]:
    banner("SECTION III — EXACT PROFILE-ANGLE EQUATION")

    K_eff, lam_minus, lam_plus, det_eff = eigenvalue_ledger()

    q = sp.Matrix([sp.cos(theta), sp.sin(theta)])
    E = sp.simplify((q.T * K_eff * q)[0] / 2)
    dE = sp.simplify(sp.expand_trig(sp.diff(E, theta)))

    rhs = sp.simplify(2 * alpha * kappa0 * kappa1 / (DeltaK + alpha * (kappa0**2 - kappa1**2)))
    # Stationarity condition written as numerator of tan(2 theta) relation.
    stationarity_expected = sp.simplify((DeltaK + alpha * (kappa0**2 - kappa1**2)) * sp.sin(2 * theta) - 2 * alpha * kappa0 * kappa1 * sp.cos(2 * theta))

    expect_zero("dE/dtheta - stationarity_expected/2", dE - stationarity_expected / 2)

    print("tan(2 theta) =", rhs)
    negative_tan = sp.simplify(-rhs)
    negative_tan_target = sp.simplify(
        2 * alpha * (-kappa0 * kappa1) / (DeltaK + alpha * (kappa0**2 - kappa1**2))
    )
    expect_zero("-tan(2 theta) - manifestly positive form", negative_tan - negative_tan_target)
    print("-tan(2 theta) =", negative_tan)

    # Weak-loading expansion: theta = alpha*k0*k1/DeltaK + O(alpha^2).
    weak = sp.series(rhs / 2, alpha, 0, 2).removeO()
    print("weak-loading theta leading term =", weak)
    expect_zero("weak-loading coefficient - k0*k1/DeltaK", weak / alpha - (kappa0 * kappa1 / DeltaK))

    return rhs, det_eff


# ---------------------------------------------------------------------------
# IV. Strong-loading limit and max-coupling branch
# ---------------------------------------------------------------------------

def strong_loading_limit() -> tuple[sp.Expr, sp.Expr]:
    banner("SECTION IV — STRONG-LOADING LIMIT AND MAX-COUPLING BRANCH")

    rhs, det_eff = angle_equation()
    rhs_inf = sp.simplify(sp.limit(rhs, alpha, sp.oo))
    tmax = sp.simplify(kappa1 / kappa0)
    tan2_tmax = sp.simplify(2 * tmax / (1 - tmax**2))

    print("lim_{alpha->oo} tan(2 theta) =", rhs_inf)
    print("tan(theta_max)               =", tmax)
    print("tan(2 theta_max)            =", tan2_tmax)

    expect_zero("strong-loading limit - tan(2 theta_max)", rhs_inf - tan2_tmax)
    expect_zero("tan(theta_max) + sqrt(2)/3", tmax + sp.sqrt(2) / 3)

    print("Thus the strong-loading eigenvector aligns with v/|v|, i.e. the Stage-10")
    print("max-coupling branch rather than the blind-angle branch.")

    return rhs_inf, det_eff


# ---------------------------------------------------------------------------
# V. Softening threshold
# ---------------------------------------------------------------------------

def softening_threshold() -> None:
    banner("SECTION V — EXACT SOFTENING THRESHOLD")

    rhs_inf, det_eff = strong_loading_limit()
    alpha_crit = sp.solve(sp.Eq(det_eff, 0), alpha)[0]
    alpha_crit_expected = sp.simplify(K0 * K1 / (K1 * kappa0**2 + K0 * kappa1**2))

    print("alpha_crit =", alpha_crit)
    expect_zero("alpha_crit - expected", alpha_crit - alpha_crit_expected)

    # Check determinant sign switch.
    eps = sp.symbols("eps", positive=True, real=True)
    det_below = sp.simplify(det_eff.subs(alpha, alpha_crit * (1 - eps)))
    det_at = sp.simplify(det_eff.subs(alpha, alpha_crit))
    print("det(alpha_crit) =", det_at)
    expect_zero("det(alpha_crit)", det_at)
    print("det(alpha_crit*(1-eps)) =")
    sp.pprint(sp.factor(det_below))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    softening_threshold()
