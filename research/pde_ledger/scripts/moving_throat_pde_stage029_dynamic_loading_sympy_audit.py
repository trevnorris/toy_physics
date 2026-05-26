#!/usr/bin/env python3
"""
moving_throat_pde_stage029_dynamic_loading_sympy_audit.py

SymPy audit for Stage 029 of the moving-throat PDE program.

Scope
-----
This script verifies the first exact dynamic-loading lift of the Stage-11 profile
selection model. It checks:

  • the full Schur-complement elimination of the coupled wall/U/W/phi block,
  • the exact decomposition of the wall self-energy into
        Xi(omega) I_2 + alpha(omega) v v^T,
  • the conservative static limits Xi_0 and alpha_0,
  • the refined Stage-11 angle law with the isotropic shift included,
  • the first-order outgoing expansion of alpha(omega),
  • and the selected-mode odd coefficient projected onto the conservative lower
    wall eigenvector.

This is the first point in the program where the Stage-11 loading parameter is
computed exactly from the coupled wall/BdG/Maxwell/mixed operator instead of
being inserted by hand.
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


def expect_zero(name: str, expr: sp.Expr | sp.Matrix) -> None:
    if isinstance(expr, sp.MatrixBase):
        simplified = expr.applyfunc(lambda z: sp.simplify(sp.expand(z)))
        print(f"{name} =")
        sp.pprint(simplified)
        if any(entry != 0 for entry in simplified):
            raise AssertionError(f"{name} is not zero")
    else:
        simplified = sp.simplify(sp.expand(expr))
        print(f"{name} = {simplified}")
        if simplified != 0:
            raise AssertionError(f"{name} is not zero")


# ---------------------------------------------------------------------------
# Symbols and exact overlap constants
# ---------------------------------------------------------------------------

omega = sp.symbols("omega", real=True)
K0, M0, M1 = sp.symbols("K0 M0 M1", positive=True, real=True)
DeltaK_ax = sp.symbols("DeltaK_ax", positive=True, real=True)
K1 = K0 + DeltaK_ax  # K1 > K0 by physical assumption (axial stiffness ordering)
varpi, Omega_U, Omega_W = sp.symbols("varpi Omega_U Omega_W", positive=True, real=True)
lambda_B, lambda_U, lambda_W, lambda_R = sp.symbols(
    "lambda_B lambda_U lambda_W lambda_R", real=True
)
Pi = sp.symbols("Pi")
Gamma_port = sp.symbols("Gamma_port", positive=True, real=True)

kappa0 = sp.simplify(2 * sp.sqrt(2) / sp.pi)
kappa1 = sp.simplify(-4 / (3 * sp.pi))
sigma = sp.simplify(kappa0**2 + kappa1**2)
xi = sp.simplify(kappa0**2 - kappa1**2)
eta = sp.simplify(kappa0 * kappa1)

v = sp.Matrix([kappa0, kappa1])
I2 = sp.eye(2)

Aphi = sp.simplify(varpi**2 - omega**2)
AU = sp.simplify(Omega_U**2 - omega**2)
AW = sp.simplify(Omega_W**2 - omega**2 - Pi)
Delta_UW = sp.simplify(AU * AW - lambda_R**2 * sigma)

Xi = sp.simplify(lambda_U**2 / AU)
alpha = sp.simplify(lambda_B**2 / Aphi + (AU * lambda_W + lambda_R * lambda_U) ** 2 / (AU * Delta_UW))


# ---------------------------------------------------------------------------
# I. Exact Schur-complement decomposition
# ---------------------------------------------------------------------------

def schur_decomposition() -> None:
    banner("SECTION I — EXACT SCHUR-COMPLEMENT DECOMPOSITION")

    # Internal fields ordered as (u0_int, u1_int, W, phi).
    Mint = sp.Matrix(
        [
            [AU, 0, -lambda_R * kappa0, 0],
            [0, AU, -lambda_R * kappa1, 0],
            [-lambda_R * kappa0, -lambda_R * kappa1, AW, 0],
            [0, 0, 0, Aphi],
        ]
    )

    # Coupling matrix from wall q=(q0,q1)^T to internal fields.
    C = sp.Matrix(
        [
            [lambda_U, 0],
            [0, lambda_U],
            [lambda_W * kappa0, lambda_W * kappa1],
            [lambda_B * kappa0, lambda_B * kappa1],
        ]
    )

    Sigma = sp.simplify(C.T * Mint.inv() * C)
    Sigma_expected = sp.simplify(Xi * I2 + alpha * (v * v.T))

    print("Sigma_wall(omega) =")
    sp.pprint(Sigma)
    print("Sigma_expected =")
    sp.pprint(Sigma_expected)
    expect_zero("Sigma - (Xi I + alpha vv^T)", Sigma - Sigma_expected)

    print("sigma =", sigma)
    print("xi    =", xi)
    print("eta   =", eta)


# ---------------------------------------------------------------------------
# II. Conservative static data and refined Stage-11 profile law
# ---------------------------------------------------------------------------

def conservative_profile_selection() -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    banner("SECTION II — CONSERVATIVE STATIC DATA AND PROFILE SELECTION")

    Xi0 = sp.simplify(Xi.subs({omega: 0, Pi: 0}))
    Delta0 = sp.simplify(Delta_UW.subs({omega: 0, Pi: 0}))
    alpha0 = sp.simplify(alpha.subs({omega: 0, Pi: 0}))

    K0t = sp.simplify(K0 - Xi0)
    K1t = sp.simplify(K1 - Xi0)
    DeltaK = sp.simplify(K1t - K0t)

    print("Xi_0    =", Xi0)
    print("Delta_0 =", Delta0)
    print("alpha_0 =", alpha0)
    # Sanity check: K0t = K0 - Xi0, K1t = K1 - Xi0, K1 = K0 + DeltaK_ax, so
    # (K1t - K0t) - DeltaK_ax = 0 is algebraically forced by the construction.
    # Kept here as a typo guard, not as a physics check.
    expect_zero("DeltaK_tilde - DeltaK_bare (Xi_0 cancellation)", DeltaK - DeltaK_ax)

    Keff0 = sp.Matrix(
        [
            [K0t - alpha0 * kappa0**2, -alpha0 * kappa0 * kappa1],
            [-alpha0 * kappa0 * kappa1, K1t - alpha0 * kappa1**2],
        ]
    )

    tr_eff = sp.simplify(sp.trace(Keff0))
    det_eff = sp.simplify(Keff0.det())
    disc = sp.simplify((DeltaK + alpha0 * xi) ** 2 + 4 * alpha0**2 * eta**2)

    lam_minus = sp.simplify((tr_eff - sp.sqrt(disc)) / 2)
    lam_plus = sp.simplify((tr_eff + sp.sqrt(disc)) / 2)

    print("lambda_- =")
    sp.pprint(lam_minus)
    print("lambda_+ =")
    sp.pprint(lam_plus)

    theta = sp.symbols("theta", real=True)
    q = sp.Matrix([sp.cos(theta), sp.sin(theta)])
    E = sp.simplify((q.T * Keff0 * q)[0] / 2)
    dE = sp.simplify(sp.expand_trig(sp.diff(E, theta)))

    stationarity = sp.simplify(
        (DeltaK + alpha0 * xi) * sp.sin(2 * theta) - 2 * alpha0 * eta * sp.cos(2 * theta)
    )
    expect_zero("dE/dtheta - stationarity/2", dE - stationarity / 2)

    tan2theta = sp.simplify(2 * alpha0 * eta / (DeltaK + alpha0 * xi))
    print("tan(2 theta_-) =", tan2theta)

    return Xi0, alpha0, lam_minus


# ---------------------------------------------------------------------------
# III. First-order outgoing dressing of alpha(omega)
# ---------------------------------------------------------------------------

def outgoing_loading() -> tuple[sp.Expr, sp.Expr]:
    banner("SECTION III — OUTGOING DRESSING OF THE DIRECTIONAL LOAD")

    # First-order expansion in the outgoing port Pi.
    alpha_cons = sp.simplify(alpha.subs(Pi, 0))
    beta = sp.simplify(sp.diff(alpha, Pi).subs(Pi, 0))
    Delta_cons = sp.simplify(AU * (Omega_W**2 - omega**2) - lambda_R**2 * sigma)
    beta_clean = sp.simplify((AU * lambda_W + lambda_R * lambda_U) ** 2 / Delta_cons**2)

    print("alpha_cons(omega) =")
    sp.pprint(alpha_cons)
    print("beta(omega) =")
    sp.pprint(beta)
    expect_zero("beta - clean transfer factor", beta - beta_clean)

    eps = sp.symbols("eps")
    alpha_series = sp.series(alpha.subs(Pi, eps), eps, 0, 2).removeO()
    expect_zero("alpha - (alpha_cons + beta*Pi) at O(Pi)", sp.expand(alpha_series.subs(eps, Pi) - (alpha_cons + beta_clean * Pi)))

    # Now impose the compact outgoing l=2 port law.
    alpha_out = sp.simplify(alpha.subs(Pi, sp.I * Gamma_port * omega**5))
    beta5 = sp.simplify(beta_clean.subs(omega, 0) * Gamma_port)

    # Series check only through O(omega^5).
    alpha_out_series = sp.series(alpha_out, omega, 0, 6).removeO()
    print("alpha_out(omega) through O(omega^5) =")
    sp.pprint(alpha_out_series)
    print("beta_5 =", beta5)
    beta5_extracted = sp.simplify(sp.expand(alpha_out_series).coeff(omega, 5) / sp.I)
    expect_zero("extracted beta_5 - expected beta_5", beta5_extracted - beta5)

    return beta5, alpha_out_series


def alpha0_symbolic() -> sp.Expr:
    """Conservative static directional load alpha_0."""
    return sp.simplify(alpha.subs({omega: 0, Pi: 0}))


# ---------------------------------------------------------------------------
# IV. Selected-mode odd coefficient
# ---------------------------------------------------------------------------

def selected_mode_projection() -> None:
    banner("SECTION IV — SELECTED-MODE ODD QUADRUPOLE COEFFICIENT")

    Xi0, alpha0, lam_minus = conservative_profile_selection()
    beta5, alpha_out_series = outgoing_loading()

    K0t = sp.simplify(K0 - Xi0)
    K1t = sp.simplify(K1 - Xi0)
    DeltaK = sp.simplify(K1t - K0t)

    theta = sp.symbols("theta", real=True)
    q = sp.Matrix([sp.cos(theta), sp.sin(theta)])
    kappa_theta_sq = sp.simplify((q.T * (v * v.T) * q)[0])
    print("kappa(theta)^2 =")
    sp.pprint(kappa_theta_sq)

    # Project the odd rank-1 piece onto a unit vector q.
    odd_projection = sp.simplify(-sp.I * beta5 * kappa_theta_sq * omega**5)
    print("Projected odd operator on unit q =")
    sp.pprint(odd_projection)

    # Hellmann-Feynman identity for the selected lower mode.
    al = sp.symbols('alpha_load', real=True)
    disc = sp.simplify((DeltaK + al * xi) ** 2 + 4 * al**2 * eta**2)
    tr_eff = sp.simplify((K0t + K1t) - al * sigma)
    lam_minus_template = sp.simplify((tr_eff - sp.sqrt(disc)) / 2)
    kappa_sel_template = sp.simplify(-sp.diff(lam_minus_template, al))
    kappa_sel_sq = sp.simplify(kappa_sel_template.subs(al, alpha0))
    print("kappa_sel^2 = - d lambda_- / d alpha |_(alpha_0) =")
    sp.pprint(kappa_sel_sq)

    # Direct eigenvector projection of v onto the lower eigenvector of K_eff(al).
    # Independent of the Hellmann-Feynman derivation above; the two must agree.
    K_eff_al = sp.Matrix(
        [
            [K0t - al * kappa0**2, -al * kappa0 * kappa1],
            [-al * kappa0 * kappa1, K1t - al * kappa1**2],
        ]
    )
    null = (K_eff_al - lam_minus_template * sp.eye(2)).nullspace()
    assert null, "lower-eigenvector nullspace is empty"
    vec_lo = null[0]
    norm_sq = sp.simplify((vec_lo.T * vec_lo)[0])
    kappa_sel_sq_direct_template = sp.simplify(((vec_lo.T * v)[0]) ** 2 / norm_sq)
    expect_zero(
        "kappa_sel^2 closed-form vs eigenvector projection",
        sp.simplify(kappa_sel_template - kappa_sel_sq_direct_template),
    )

    # Verify limiting values on the template before substitution.
    expect_zero("weak-loading kappa_sel^2 -> kappa0^2", sp.simplify(kappa_sel_template.subs(al, 0) - kappa0**2))
    expect_zero("strong-loading kappa_sel^2 -> sigma", sp.simplify(sp.limit(kappa_sel_template, al, sp.oo) - sigma))

    # Anchor the paper's eq selected-odd against the script-computed pieces.
    delta_D_paper = (
        -sp.I * Gamma_port
        * (Omega_U**2 * lambda_W + lambda_R * lambda_U) ** 2
        / (Omega_U**2 * Omega_W**2 - lambda_R**2 * sigma) ** 2
        * kappa_sel_sq
        * omega**5
    )
    delta_D_script = -sp.I * beta5 * kappa_sel_sq * omega**5
    expect_zero(
        "delta D_-^odd (script) - delta D_-^odd (paper formula)",
        delta_D_script - delta_D_paper,
    )

    print("Therefore the selected lower-mode odd coefficient is")
    print("  delta D_-^(odd)(omega) = - i beta_5 kappa_sel^2 omega^5 + O(omega^7)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    schur_decomposition()
    selected_mode_projection()
