#!/usr/bin/env python3
"""SymPy audit for Moving-Throat PDE Stage 224.

This script verifies the exact primitive ranking on the selected twin-support
branch:

1.  Selected-branch reduction
        eps_* = 1 - 3 varrho / 2,
        sigma = 4/(3 varrho) - 2
2.  The selected curve lies strictly inside the symmetric lowest-twin window.
3.  Primitive coherent weight identities:
        w_chi = 2 w_Z,
        w_Z > w_Lambda,
        w_Z > w_W,
        w_W > |w_U|
4.  Exact selected-branch thresholds:
        varrho_WLambda = 2(1+beta^2)/[3(2+beta^2)],
        varrho_ULambda = 2(1+beta^2)/[3(1+beta+beta^2)]
5.  Exact sign-transfer laws:
        w_Lambda - w_W     ~ (varrho_WLambda - varrho),
        w_Lambda - |w_U|   ~ (varrho_ULambda - varrho)
6.  Exact ordering:
        0 < varrho_WLambda < varrho_ULambda < 2/3
7.  Numerical threshold windows implied by 0 < beta < 2/11.
"""

from __future__ import annotations

import sympy as sp


def assert_zero(expr: sp.Expr, label: str) -> None:
    simplified = sp.simplify(sp.factor(expr))
    if simplified != 0:
        raise AssertionError(f"{label} failed: {simplified}")
    print(f"[ok] {label}")


def assert_positive_sample(expr: sp.Expr, subs_map: dict[sp.Symbol, sp.Expr], label: str) -> None:
    value = sp.N(expr.subs(subs_map))
    if not (value > 0):
        raise AssertionError(f"{label} is not positive at sample point: {value}")
    print(f"[ok] {label} > 0 at sample point")


def main() -> None:
    sp.init_printing()

    # ------------------------------------------------------------------
    # Symbols
    # ------------------------------------------------------------------
    beta, varrho = sp.symbols("beta varrho", positive=True)
    eps_star, sigma = sp.symbols("eps_star sigma", positive=True)

    # ------------------------------------------------------------------
    # 1. Selected-branch reduction
    # ------------------------------------------------------------------
    eps_sel = sp.simplify(1 - sp.Rational(3, 2) * varrho)
    sigma_sel = sp.simplify(2 * eps_sel / (1 - eps_sel))

    assert_zero(eps_sel - (1 - sp.Rational(3, 2) * varrho), "eps_* = 1 - 3 varrho / 2")
    assert_zero(sigma_sel - (4 / (3 * varrho) - 2), "sigma = 4/(3 varrho) - 2")

    # ------------------------------------------------------------------
    # 2. Twin-window inclusion
    # ------------------------------------------------------------------
    mixed_upper = 1 / varrho - 2
    twin_upper = 2 / varrho - 2

    assert_zero(
        sp.simplify(sigma_sel - mixed_upper - 1 / (3 * varrho)),
        "sigma_sel - (1/varrho - 2) = 1/(3 varrho)",
    )
    assert_zero(
        sp.simplify(twin_upper - sigma_sel - 2 / (3 * varrho)),
        "(2/varrho - 2) - sigma_sel = 2/(3 varrho)",
    )

    # ------------------------------------------------------------------
    # 3. Primitive coherent weights
    # ------------------------------------------------------------------
    N = 5 * (1 - eps_star) ** 2 + 6 * eps_star**2 * (1 + beta**2)

    w_Lambda = eps_star**2 * (1 + beta**2) / N
    w_Z = (1 - 2 * eps_star + (2 + beta**2) * eps_star**2) / N
    w_chi = 2 * (1 - 2 * eps_star + (2 + beta**2) * eps_star**2) / N
    w_W = eps_star * (1 - eps_star) / N
    w_Umag = beta * eps_star * (1 - eps_star) / N

    assert_zero(w_chi - 2 * w_Z, "w_chi = 2 w_Z")
    assert_zero(
        sp.simplify(w_Z - w_Lambda - (1 - eps_star) ** 2 / N),
        "w_Z - w_Lambda = (1-eps_*)^2 / N",
    )
    assert_zero(
        sp.simplify(
            w_Z
            - w_W
            - (beta**2 * eps_star**2 + 3 * (eps_star - sp.Rational(1, 2)) ** 2 + sp.Rational(1, 4)) / N
        ),
        "w_Z - w_W exact positive form",
    )
    assert_zero(
        sp.simplify(w_W - w_Umag - eps_star * (1 - eps_star) * (1 - beta) / N),
        "w_W - |w_U| exact positive form",
    )

    # Pull back to the selected twin-support curve.
    wL_sel = sp.simplify(w_Lambda.subs(eps_star, eps_sel))
    wZ_sel = sp.simplify(w_Z.subs(eps_star, eps_sel))
    wchi_sel = sp.simplify(w_chi.subs(eps_star, eps_sel))
    wW_sel = sp.simplify(w_W.subs(eps_star, eps_sel))
    wU_sel = sp.simplify(w_Umag.subs(eps_star, eps_sel))

    # ------------------------------------------------------------------
    # 4. Exact selected-branch thresholds
    # ------------------------------------------------------------------
    varrho_WLambda = sp.simplify(2 * (1 + beta**2) / (3 * (2 + beta**2)))
    varrho_ULambda = sp.simplify(2 * (1 + beta**2) / (3 * (1 + beta + beta**2)))

    assert_zero(
        sp.simplify(eps_sel - 1 / (2 + beta**2) - sp.Rational(3, 2) * (varrho_WLambda - varrho)),
        "selected-branch q_W = q_Lambda threshold",
    )
    assert_zero(
        sp.simplify(
            eps_sel
            - beta / (1 + beta + beta**2)
            - sp.Rational(3, 2) * (varrho_ULambda - varrho)
        ),
        "selected-branch |q_U| = q_Lambda threshold",
    )

    # ------------------------------------------------------------------
    # 5. Factorized sign-transfer laws
    # ------------------------------------------------------------------
    D = 18 * beta**2 * varrho**2 - 24 * beta**2 * varrho + 8 * beta**2 + 33 * varrho**2 - 24 * varrho + 8

    assert_zero(
        sp.simplify(
            wL_sel - wW_sel
            - (2 - 3 * varrho) * (2 + beta**2) * (varrho_WLambda - varrho) / D
        ),
        "factorized sign law for w_Lambda - w_W",
    )
    assert_zero(
        sp.simplify(
            wL_sel - wU_sel
            - (2 - 3 * varrho) * (1 + beta + beta**2) * (varrho_ULambda - varrho) / D
        ),
        "factorized sign law for w_Lambda - |w_U|",
    )

    # ------------------------------------------------------------------
    # 6. Exact ordering of the thresholds
    # ------------------------------------------------------------------
    assert_zero(
        sp.simplify(
            varrho_ULambda
            - varrho_WLambda
            - 2 * (1 + beta**2) * (1 - beta) / (3 * (1 + beta + beta**2) * (2 + beta**2))
        ),
        "varrho_ULambda - varrho_WLambda exact positive form",
    )
    assert_zero(
        sp.simplify(
            sp.Rational(2, 3)
            - varrho_ULambda
            - 2 * beta / (3 * (1 + beta + beta**2))
        ),
        "2/3 - varrho_ULambda exact positive form",
    )

    # ------------------------------------------------------------------
    # 7. Numerical windows from 0 < beta < 2/11
    # ------------------------------------------------------------------
    beta_max = sp.Rational(2, 11)

    assert_zero(
        varrho_WLambda.subs(beta, 0) - sp.Rational(1, 3),
        "varrho_WLambda(beta=0) = 1/3",
    )
    assert_zero(
        varrho_WLambda.subs(beta, beta_max) - sp.Rational(125, 369),
        "varrho_WLambda(beta=2/11) = 125/369",
    )
    assert_zero(
        varrho_ULambda.subs(beta, 0) - sp.Rational(2, 3),
        "varrho_ULambda(beta=0) = 2/3",
    )
    assert_zero(
        varrho_ULambda.subs(beta, beta_max) - sp.Rational(250, 441),
        "varrho_ULambda(beta=2/11) = 250/441",
    )

    assert_zero(
        sp.simplify(sp.diff(varrho_WLambda, beta) - 4 * beta / (3 * (beta**2 + 2) ** 2)),
        "d varrho_WLambda / d beta exact form",
    )
    assert_zero(
        sp.simplify(
            sp.diff(varrho_ULambda, beta) + 2 * (1 - beta**2) / (3 * (1 + beta + beta**2) ** 2)
        ),
        "d varrho_ULambda / d beta exact form",
    )

    # ------------------------------------------------------------------
    # 8. Representative region checks
    # ------------------------------------------------------------------
    beta_sample = sp.Rational(1, 10)
    rhoW_sample = sp.simplify(varrho_WLambda.subs(beta, beta_sample))
    rhoU_sample = sp.simplify(varrho_ULambda.subs(beta, beta_sample))

    rho_region_1 = sp.simplify(rhoW_sample / 2)
    rho_region_2 = sp.simplify((rhoW_sample + rhoU_sample) / 2)
    rho_region_3 = sp.simplify((rhoU_sample + sp.Rational(2, 3)) / 2)

    # Region I: q_chi > q_Z > q_Lambda > q_W > |q_U|
    subs1 = {beta: beta_sample, varrho: rho_region_1}
    assert_positive_sample(wchi_sel - wZ_sel, subs1, "Region I: w_chi - w_Z")
    assert_positive_sample(wZ_sel - wL_sel, subs1, "Region I: w_Z - w_Lambda")
    assert_positive_sample(wL_sel - wW_sel, subs1, "Region I: w_Lambda - w_W")
    assert_positive_sample(wW_sel - wU_sel, subs1, "Region I: w_W - |w_U|")

    # Region II: q_chi > q_Z > q_W > q_Lambda > |q_U|
    subs2 = {beta: beta_sample, varrho: rho_region_2}
    assert_positive_sample(wchi_sel - wZ_sel, subs2, "Region II: w_chi - w_Z")
    assert_positive_sample(wZ_sel - wW_sel, subs2, "Region II: w_Z - w_W")
    assert_positive_sample(wW_sel - wL_sel, subs2, "Region II: w_W - w_Lambda")
    assert_positive_sample(wL_sel - wU_sel, subs2, "Region II: w_Lambda - |w_U|")

    # Region III: q_chi > q_Z > q_W > |q_U| > q_Lambda
    subs3 = {beta: beta_sample, varrho: rho_region_3}
    assert_positive_sample(wchi_sel - wZ_sel, subs3, "Region III: w_chi - w_Z")
    assert_positive_sample(wZ_sel - wW_sel, subs3, "Region III: w_Z - w_W")
    assert_positive_sample(wW_sel - wU_sel, subs3, "Region III: w_W - |w_U|")
    assert_positive_sample(wU_sel - wL_sel, subs3, "Region III: |w_U| - w_Lambda")

    print("\nAll Stage-224 symbolic checks passed.")


if __name__ == "__main__":
    main()
