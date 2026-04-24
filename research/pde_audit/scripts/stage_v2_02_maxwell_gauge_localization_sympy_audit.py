#!/usr/bin/env python3
"""
Stage V2-02 — Maxwell gauge-localization audit for the moving-throat PDE program.

This script checks the algebra behind the localized 4+1 Maxwell sector with a
possibly weighted Lorenz gauge-fixing functional

    L_gf = - H(w)/(2 xi mu0) (partial_M A^M)^2.

The current papers use H(w)=1 while the Maxwell kinetic term uses Z(w).  The
script verifies:

1. the EOM and divergence identity for general H(w), including H=1 and H=Z;
2. the zero-mode reduction of the kinetic and gauge-fixing coefficients;
3. the divergence of the unweighted zero-mode gauge-fixing norm on a noncompact
   w-line;
4. the finite weighted/locally supported alternatives and their effective 4D
   gauge parameter;
5. gauge invariance of the mixed-sector observables E_w and C_a.

The output is intended as a deterministic derivation ledger, not a numerical fit.
"""

from __future__ import annotations

import sympy as sp


def header(title: str) -> None:
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


def check_zero(name: str, expr: sp.Expr) -> None:
    simplified = sp.simplify(expr)
    print(f"{name}: {simplified}")
    assert simplified == 0, f"FAILED: {name} -> {simplified}"


def main() -> None:
    # ------------------------------------------------------------------
    # 1. Gaussian localization integrals and zero-mode coefficients.
    # ------------------------------------------------------------------
    header("1. Localization integrals and zero-mode gauge-fixing norm")

    w = sp.symbols("w", real=True)
    lam = sp.symbols("lambda", positive=True, finite=True)
    R = sp.symbols("R", positive=True, finite=True)
    xi = sp.symbols("xi", positive=True, finite=True)
    mu0 = sp.symbols("mu0", positive=True, finite=True)

    Z = sp.exp(-w**2 / lam**2)
    Z_int = sp.integrate(Z, (w, -sp.oo, sp.oo))
    H_unweighted_int_R = sp.integrate(1, (w, -R, R))
    H_weighted_int = Z_int
    ratio_unweighted_to_kinetic = sp.simplify(H_unweighted_int_R / Z_int)
    ratio_weighted_to_kinetic = sp.simplify(H_weighted_int / Z_int)
    xi4_regulated_unweighted = sp.simplify(xi * Z_int / H_unweighted_int_R)
    xi4_weighted = sp.simplify(xi * Z_int / H_weighted_int)

    print(f"Z(w) = {Z}")
    print(f"Z_int = ∫ Z dw = {Z_int}")
    print(f"H_int[unweighted; -R,R] = {H_unweighted_int_R}")
    print(f"H_int[weighted H=Z] = {H_weighted_int}")
    print(f"H_unweighted/Z_int = {ratio_unweighted_to_kinetic}")
    print(f"limit_{'{'}R->oo{'}'} H_unweighted/Z_int = {sp.limit(ratio_unweighted_to_kinetic, R, sp.oo)}")
    print(f"H_weighted/Z_int = {ratio_weighted_to_kinetic}")
    print()
    print("For a finite gauge-fixing weight H_int, matching")
    print("    -H_int/(2 xi mu0) (∂·a)^2 = -Z_int/(2 xi_4 mu0) (∂·a)^2")
    print("gives xi_4 = xi Z_int/H_int.")
    print(f"xi_4[unweighted box -R,R] = {xi4_regulated_unweighted}")
    print(f"limit_{'{'}R->oo{'}'} xi_4[unweighted] = {sp.limit(xi4_regulated_unweighted, R, sp.oo)}")
    print(f"xi_4[weighted H=Z] = {xi4_weighted}")

    assert Z_int == sp.sqrt(sp.pi) * lam
    assert sp.limit(ratio_unweighted_to_kinetic, R, sp.oo) == sp.oo
    assert xi4_weighted == xi

    # ------------------------------------------------------------------
    # 2. Divergence of the Maxwell kinetic operator: antisymmetry check.
    # ------------------------------------------------------------------
    header("2. Antisymmetry check: divergence of ∂_M(Z F^{MN})")

    t, x = sp.symbols("t x", real=True)
    F_xw = sp.Function("F_xw")(t, x, w)
    Zfun = sp.Function("Z")(w)

    # In the two-coordinate subblock (x,w), take F^{xw}=F and F^{wx}=-F.
    # The double divergence is ∂_x∂_w(Z F^{wx}) + ∂_w∂_x(Z F^{xw}).
    double_div = sp.diff(sp.diff(Zfun * (-F_xw), w), x) + sp.diff(sp.diff(Zfun * F_xw, x), w)
    check_zero("∂_N∂_M(Z F^{MN}) in an antisymmetric 2D subblock", double_div)

    # ------------------------------------------------------------------
    # 3. General gauge-fixing weight H(w): EOM and divergence identity.
    # ------------------------------------------------------------------
    header("3. General gauge-fixing weight H(w): divergence identity")

    H = sp.Function("H")(w)
    B = sp.Function("B")(t, x, w)  # B = ∂·A
    box_B = -sp.diff(B, t, 2) + sp.diff(B, x, 2) + sp.diff(B, w, 2)
    box_HB = -sp.diff(H * B, t, 2) + sp.diff(H * B, x, 2) + sp.diff(H * B, w, 2)
    expanded_box_HB = sp.expand(box_HB)
    expected_box_HB = sp.expand(H * box_B + 2 * sp.diff(H, w) * sp.diff(B, w) + sp.diff(H, w, 2) * B)
    check_zero("□(H B) - [H □B + 2 H' B_w + H'' B]", expanded_box_HB - expected_box_HB)

    print("General weighted gauge-fixing EOM term:")
    print("    +(1/xi) ∂^N[H(w) B],  where B = ∂·A")
    print("General divergence consistency condition:")
    print("    (1/xi) □[H(w) B] = mu0 ∂_N J^N")
    print("For H=1, this reduces to (1/xi) □B = mu0 ∂·J.")
    print("For H=Z, the equation is (1/xi) □[Z B] = mu0 ∂·J.")

    # Unweighted and weighted specializations.
    H1 = sp.Integer(1)
    check_zero("H=1 specialization extra H' terms", (2 * sp.diff(H1, w) * sp.diff(B, w) + sp.diff(H1, w, 2) * B))
    HZ = Z
    extra_Z_terms = sp.simplify(2 * sp.diff(HZ, w) * sp.diff(B, w) + sp.diff(HZ, w, 2) * B)
    print(f"For Gaussian H=Z, extra w-localization terms in □(ZB) are: {extra_Z_terms}")

    # ------------------------------------------------------------------
    # 4. Zero-mode action audit.
    # ------------------------------------------------------------------
    header("4. Zero-mode action audit")

    b, f2 = sp.symbols("b f2")  # b=∂_mu a^mu, f2=f_{mu nu}f^{mu nu}
    S_kin_zero = -Z_int * f2 / (4 * mu0)
    S_gf_unweighted_R = -H_unweighted_int_R * b**2 / (2 * xi * mu0)
    S_gf_weighted = -Z_int * b**2 / (2 * xi * mu0)
    print(f"Zero-mode kinetic density coefficient: {S_kin_zero}")
    print(f"Zero-mode unweighted gauge-fixing coefficient on [-R,R]: {S_gf_unweighted_R}")
    print(f"Zero-mode weighted H=Z gauge-fixing coefficient: {S_gf_weighted}")
    print("Noncompact verdict:")
    print("    H=1 gives a divergent gauge-fixing action for any zero mode with b=∂_mu a^mu != 0.")
    print("    H=1 is finite only if Lorenz gauge b=0 is imposed before reduction, or if the w-line is regulated.")
    print("    H=Z gives a finite localized gauge-fixed zero-mode action with the same Z_int normalization as the Maxwell kinetic term.")

    # Verify general effective xi matching.
    H_int = sp.symbols("H_int", positive=True)
    xi4 = xi * Z_int / H_int
    match_expr = H_int / xi - Z_int / xi4
    check_zero("general 4D gauge parameter match H_int/xi = Z_int/xi_4", match_expr)

    # ------------------------------------------------------------------
    # 5. Gauge invariance of mixed-sector observables.
    # ------------------------------------------------------------------
    header("5. Gauge invariance of mixed-sector observables")

    chi = sp.Function("chi")(t, x, w)
    A0 = sp.Function("A0")(t, x, w)
    Aa = sp.Function("Aa")(t, x, w)
    Aw = sp.Function("Aw")(t, x, w)

    A0_p = A0 - sp.diff(chi, t)
    Aa_p = Aa + sp.diff(chi, x)
    Aw_p = Aw + sp.diff(chi, w)

    Ew = -sp.diff(Aw, t) - sp.diff(A0, w)
    Ew_p = -sp.diff(Aw_p, t) - sp.diff(A0_p, w)
    Ca = sp.diff(Aw, x) - sp.diff(Aa, w)
    Ca_p = sp.diff(Aw_p, x) - sp.diff(Aa_p, w)

    check_zero("delta E_w", Ew_p - Ew)
    check_zero("delta C_a", Ca_p - Ca)

    # ------------------------------------------------------------------
    # 6. Final pass/fail verdicts.
    # ------------------------------------------------------------------
    header("6. Stage V2-02 verdict")

    verdicts = [
        ("Bulk EOM with H=1 as written", "PASS as a formal gauge-fixed bulk equation"),
        ("Divergence consistency for H=1", "PASS if the total current is conserved and B obeys □B=0 / Lorenz gauge is imposed"),
        ("Noncompact zero-mode gauge-fixed action with H=1", "FAIL unless ∂_mu a^mu=0 is imposed before reduction or the w-line is regulated"),
        ("Clean gauge-fixed brane reduction with H=1", "CONDITIONAL: impose Lorenz gauge before reduction, then choose any 3+1 gauge fixing afterward"),
        ("Weighted localized gauge fixing H=Z", "PASS for finite zero-mode normalization; modifies the gauge-condition propagation to □(ZB)=0"),
        ("Mixed-sector gauge invariants E_w, C_a", "PASS; they are exact gauge-invariant observables"),
    ]
    for item, verdict in verdicts:
        print(f"{item}: {verdict}")

    print("\nRecommended Volume-2 patch:")
    print("    Either (A) keep the current H=1 term only as a bulk gauge-fixing device and")
    print("        state explicitly that Lorenz gauge is imposed before the zero-mode reduction,")
    print("        with 3+1 gauge fixing chosen only after reduction; or")
    print("    (B) replace it in the parent ledger by a localized gauge-fixing functional")
    print("        H(w)=Z(w) or another finite profile H_loc(w), so the gauge-fixed")
    print("        zero-mode action is finite without a pre-reduction Lorenz constraint.")


if __name__ == "__main__":
    main()
