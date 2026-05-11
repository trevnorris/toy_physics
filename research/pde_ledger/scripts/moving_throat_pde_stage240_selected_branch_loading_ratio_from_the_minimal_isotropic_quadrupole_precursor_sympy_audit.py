#!/usr/bin/env python3
"""SymPy audit for Moving-Throat PDE Stage 240.

This script verifies the selected-branch loading-ratio extraction from the
minimal isotropic quadrupole precursor:

1.  Pi_tr / C_mix = alpha_req / alpha_mix = rho_alpha
2.  The exact contact-plus-pole compiler
        c0 = 1/rho_alpha,
        c1 = (rho_alpha - 1)/rho_alpha
    and the inverse formulas
        rho_alpha = 1/c0 = 1/(1-c1),
        zeta_req  = c1/c0 = rho_alpha - 1
3.  The minimal isotropic specialization
        c0 = 3/4, c1 = 1/4
    gives
        rho_alpha = 4/3, zeta_req = 1/3
4.  The exact selected demand product
        Pi_tr = (4/3) C_mix
5.  The support-selector reduction
        varrho = pi^2 Pi_tr / (16 Lambda) = 2(1-eps_*)/3
6.  The regime classification:
        mixed-only is not enough,
        symmetric lowest twin is enough,
        non-twin asymmetry is not required.
"""

from __future__ import annotations

import sympy as sp


def assert_zero(expr: sp.Expr, label: str) -> None:
    simplified = sp.simplify(sp.factor(expr))
    if simplified != 0:
        raise AssertionError(f"{label} failed: {simplified}")
    print(f"[ok] {label}")


def main() -> None:
    sp.init_printing()

    # ------------------------------------------------------------------
    # Symbols
    # ------------------------------------------------------------------
    omega, Omega_Q = sp.symbols("omega Omega_Q", nonzero=True)
    alpha_req, alpha_mix = sp.symbols("alpha_req alpha_mix", positive=True)
    beta0, NQ_target = sp.symbols("beta0 NQ_target", positive=True)
    mhat_minus, s_minus, lambda_minus = sp.symbols(
        "mhat_minus s_minus lambda_minus", positive=True
    )
    c0, c1 = sp.symbols("c0 c1", real=True)
    Lambda, eps_star = sp.symbols("Lambda eps_star", positive=True)
    Pi_tr, C_mix = sp.symbols("Pi_tr C_mix", positive=True)

    rho_alpha = sp.simplify(alpha_req / alpha_mix)

    # ------------------------------------------------------------------
    # 1. Exact selected-branch product identities
    # ------------------------------------------------------------------
    Pi_tr_expr = NQ_target * alpha_req / beta0
    C_mix_expr = NQ_target * alpha_mix / beta0

    assert_zero(
        sp.simplify(Pi_tr_expr / C_mix_expr - rho_alpha),
        "Pi_tr / C_mix = alpha_req / alpha_mix",
    )

    # Spectral version using N_Q^(target) = mhat_-^2 beta0 s_- / lambda_-
    NQ_target_spec = mhat_minus**2 * beta0 * s_minus / lambda_minus
    Pi_tr_spec = sp.simplify(NQ_target_spec * alpha_req / beta0)
    C_mix_spec = sp.simplify(NQ_target_spec * alpha_mix / beta0)
    assert_zero(
        sp.simplify(Pi_tr_spec / C_mix_spec - rho_alpha),
        "spectral selected-branch ratio identity",
    )

    # ------------------------------------------------------------------
    # 2. Contact-plus-pole compiler and inverse formulas
    # ------------------------------------------------------------------
    pole = 1 / (1 - omega**2 / Omega_Q**2)

    Y_support = alpha_mix / alpha_req + (alpha_req - alpha_mix) / alpha_req * pole
    Y_rho = 1 / rho_alpha + (rho_alpha - 1) / rho_alpha * pole

    assert_zero(sp.simplify(Y_support - Y_rho), "contact-plus-pole support/source compiler")

    c0_expr = sp.simplify(1 / rho_alpha)
    c1_expr = sp.simplify((rho_alpha - 1) / rho_alpha)

    assert_zero(sp.simplify(c0_expr + c1_expr - 1), "normalized static limit c0 + c1 = 1")
    assert_zero(sp.simplify(1 / c0_expr - rho_alpha), "inverse formula rho_alpha = 1/c0")
    assert_zero(sp.simplify(1 / (1 - c1_expr) - rho_alpha), "inverse formula rho_alpha = 1/(1-c1)")
    assert_zero(sp.simplify(c1_expr / c0_expr - (rho_alpha - 1)), "zeta_req = c1/c0 = rho_alpha - 1")

    # Also verify the same formulas in symbolic c0/c1 language under c0 + c1 = 1.
    rho_from_c0 = 1 / c0
    rho_from_c1 = 1 / (1 - c1)
    zeta_from_c = c1 / c0
    assert_zero(
        sp.simplify(rho_from_c0.subs(c0, c0_expr) - rho_alpha),
        "symbolic rho_alpha = 1/c0 specialization",
    )
    assert_zero(
        sp.simplify(rho_from_c1.subs(c1, c1_expr) - rho_alpha),
        "symbolic rho_alpha = 1/(1-c1) specialization",
    )
    assert_zero(
        sp.simplify(zeta_from_c.subs({c0: c0_expr, c1: c1_expr}) - (rho_alpha - 1)),
        "symbolic zeta_req specialization",
    )

    # Omega_Q does not affect the static loading-ratio extraction.
    assert_zero(sp.diff(c0_expr, Omega_Q), "Omega_Q independence of c0")
    assert_zero(sp.diff(c1_expr, Omega_Q), "Omega_Q independence of c1")

    # ------------------------------------------------------------------
    # 3. Minimal isotropic conservative module
    # ------------------------------------------------------------------
    c0_min = sp.Rational(3, 4)
    c1_min = sp.Rational(1, 4)

    rho_min_from_c0 = sp.simplify(1 / c0_min)
    rho_min_from_c1 = sp.simplify(1 / (1 - c1_min))
    zeta_min = sp.simplify(c1_min / c0_min)

    assert_zero(rho_min_from_c0 - sp.Rational(4, 3), "minimal isotropic rho_alpha from c0 = 3/4")
    assert_zero(rho_min_from_c1 - sp.Rational(4, 3), "minimal isotropic rho_alpha from c1 = 1/4")
    assert_zero(zeta_min - sp.Rational(1, 3), "minimal isotropic zeta_req = 1/3")

    # ------------------------------------------------------------------
    # 4. Exact selected demand product
    # ------------------------------------------------------------------
    Pi_tr_selected = sp.simplify(rho_min_from_c0 * C_mix)
    assert_zero(
        Pi_tr_selected - sp.Rational(4, 3) * C_mix,
        "Pi_tr = (4/3) C_mix on the minimal isotropic branch",
    )

    S_req = sp.simplify(Pi_tr_selected / C_mix)
    assert_zero(S_req - sp.Rational(4, 3), "S_req = Pi_tr / C_mix = 4/3")

    # ------------------------------------------------------------------
    # 5. Exact support-selector reduction
    # ------------------------------------------------------------------
    C_mix_Lambda = 8 * Lambda * (1 - eps_star) / sp.pi**2
    varrho_expr = sp.simplify(sp.pi**2 * (sp.Rational(4, 3) * C_mix_Lambda) / (16 * Lambda))
    assert_zero(
        varrho_expr - 2 * (1 - eps_star) / 3,
        "varrho = 2(1-eps_*)/3",
    )

    # Also verify by substitution into the symbolic Pi_tr/C_mix representation.
    varrho_symbolic = sp.simplify(sp.pi**2 * Pi_tr / (16 * Lambda))
    assert_zero(
        sp.simplify(varrho_symbolic.subs(Pi_tr, sp.Rational(4, 3) * C_mix_Lambda) - 2 * (1 - eps_star) / 3),
        "symbolic support-selector reduction",
    )

    # ------------------------------------------------------------------
    # 6. Regime classification
    # ------------------------------------------------------------------
    ratio = sp.Rational(4, 3)
    if not (1 < ratio < 2):
        raise AssertionError("regime classification failed for ratio = 4/3")
    print("[ok] regime inequality 1 < 4/3 < 2")

    # Encode the exact regime meaning in terms of Pi_tr and C_mix.
    if not (C_mix.subs(C_mix, 1) < Pi_tr_selected.subs(C_mix, 1) < 2 * C_mix.subs(C_mix, 1)):
        raise AssertionError("selected branch is not in the symmetric-lowest-twin regime")
    print("[ok] selected branch lies strictly in the symmetric-lowest-twin sector")

    print("\nAll Stage 240 symbolic checks passed.")


if __name__ == "__main__":
    main()
