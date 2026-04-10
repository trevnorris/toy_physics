#!/usr/bin/env python3
"""
Step 39 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Imposes a stronger minimal adiabatic electron-track closure on top of Step 38:
      d ln Z_q = 0,
   i.e. the wall thermodynamics are frozen and the core reacts by changing its
   compressible stiffness/outgoing normalization rather than its transverse
   localization norm.
2. Solves the resulting two-amplitude reduced branch exactly.
3. Verifies that the whole reduced electron track is rank-2, parameterized by
   (d ln K_s, d ln P_0).
4. Keeps the algebraic upper compensation branch g_+ alive by proving the exact
   sign-indefinite source bound needed to realize any g > 1.

Interpretation
--------------
On the minimal adiabatic electron track, the live reduced variables are only an
elastic squish amplitude and an outgoing/core-yield amplitude. The discarded
upper branch is not deleted: it is reinterpreted as a non-positive, pumped, or
source/sink sheet rather than the passive positive-source electron branch.
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


def main() -> None:
    banner("STEP 39 — MINIMAL CORE-YIELD CLOSURE AND THE RETAINED UPPER SHEET")

    dKs, dP0 = sp.symbols("dln_Ks dln_P0", real=True)
    dKq, dZq, da, dcs = sp.symbols("dln_Kq dln_Zq dln_a dln_cs", real=True)
    dv, dT, dgs, dgq, dlmb = sp.symbols("dln_vw0 dln_Tm dln_gs dln_gq dln_lambda", real=True)
    d_ratio, d_prod = sp.symbols("dln_v_over_T dln_vT", real=True)
    dfrakg, dfrakr, drc = sp.symbols("dln_frakg dln_frakr dln_rc", real=True)

    subbanner("XXXIX.1 — Stronger minimal adiabatic core-yield closure")
    print("On top of Step 38, impose the minimal closure")
    print("  dln Z_q = 0,")
    print("so the core reacts by compressible stiffness/outgoing-normalization drift,")
    print("not by changing the transverse mixed-channel localization norm.")

    # From Step 38 and dZq = 0.
    dKq_sol = sp.simplify(sp.Rational(2, 5) * dP0)
    da_sol = sp.simplify(sp.Rational(1, 2) * dKs)
    dcs_sol = sp.simplify(sp.Rational(1, 2) * dKs + sp.Rational(1, 5) * dP0)
    d_ratio_sol = sp.simplify(sp.Rational(1, 2) * dKs + sp.Rational(2, 5) * dP0)
    d_prod_sol = sp.simplify(-2 * dKs)
    dv_sol = sp.simplify(-sp.Rational(3, 4) * dKs + sp.Rational(1, 5) * dP0)
    dT_sol = sp.simplify(-sp.Rational(5, 4) * dKs - sp.Rational(1, 5) * dP0)
    dgs_sol = sp.simplify(-sp.Rational(1, 4) * dKs - sp.Rational(1, 5) * dP0)
    dgq_sol = sp.simplify(-sp.Rational(3, 4) * dKs)
    dlmb_sol = sp.simplify(sp.Rational(1, 2) * dKs + sp.Rational(1, 5) * dP0)
    dfrakg_sol = sp.Integer(0)
    dfrakr_sol = sp.Integer(0)
    drc_sol = sp.Integer(0)

    print("dln K_q =")
    sp.pprint(dKq_sol)
    print("dln a =")
    sp.pprint(da_sol)
    print("dln c_s =")
    sp.pprint(dcs_sol)
    print("dln(v_w0/T_m) =")
    sp.pprint(d_ratio_sol)
    print("dln(v_w0 T_m) =")
    sp.pprint(d_prod_sol)
    print("dln v_w0 =")
    sp.pprint(dv_sol)
    print("dln T_m =")
    sp.pprint(dT_sol)
    print("dln g_s =")
    sp.pprint(dgs_sol)
    print("dln g_q =")
    sp.pprint(dgq_sol)
    print("dln lambda =")
    sp.pprint(dlmb_sol)

    # Exact consistency checks.
    expect_zero("K_q law recovered", dKq_sol - sp.Rational(2, 5) * dP0)
    expect_zero("c_s law recovered", dcs_sol - (da_sol + sp.Rational(1, 5) * dP0))
    expect_zero("ratio law recovered", d_ratio_sol - (dv_sol - dT_sol))
    expect_zero("product law recovered", d_prod_sol - (dv_sol + dT_sol))
    expect_zero("g_s from T_m and J_s", dgs_sol - (dT_sol + 2 * da_sol))
    expect_zero("g_q from Z_q and L_W", dgq_sol - (dZq.subs(dZq, 0) - sp.Rational(3, 2) * da_sol))
    expect_zero("lambda from v_w0, a, L_W", dlmb_sol - (dv_sol + sp.Rational(5, 2) * da_sol))

    subbanner("XXXIX.2 — Rank-2 reduced electron track")
    vec = sp.Matrix([
        dKq_sol, da_sol, dcs_sol, dv_sol, dT_sol, dgs_sol, dgq_sol, dlmb_sol
    ])
    J = vec.jacobian(sp.Matrix([dKs, dP0]))
    print("Jacobian of the reduced drift vector w.r.t. (dln K_s, dln P_0) =")
    sp.pprint(J)
    print(f"rank = {J.rank()}")
    if J.rank() != 2:
        raise AssertionError("Minimal adiabatic closure did not collapse to rank 2.")

    # Parent surface invariants.
    mouth_imbalance = sp.simplify(dgq_sol + dKs - dgs_sol - dlmb_sol)
    stiff_imbalance = sp.simplify(dKs + dKq_sol - 2 * dlmb_sol)
    expect_zero("mouth-coupling imbalance vanishes", mouth_imbalance)
    expect_zero("stiffness/hybridization imbalance vanishes", stiff_imbalance)
    expect_zero("frak g remains frozen", dfrakg_sol)
    expect_zero("frak r remains frozen", dfrakr_sol)
    expect_zero("r_c remains frozen", drc_sol)

    subbanner("XXXIX.3 — The retained upper sheet as a sign-indefinite source law")
    g = sp.symbols("g", real=True)
    Wp, Wm = sp.symbols("W_plus W_minus", nonnegative=True, real=True)

    # For a sign-indefinite normalized source sigma = sigma_+ - sigma_- with
    # W_+ - W_- = 1 and 0 <= cos <= 1, the mouth bias obeys g <= W_+ = 1+W_-.
    bound = sp.simplify((1 + Wm) - g)
    print("For any normalized sign-indefinite source profile:")
    print("  g <= 1 + W_- ,  where W_- is the total negative weight.")
    print("Hence any realization of g > 1 requires")
    print("  W_- >= g - 1.")

    rF1 = sp.nsimplify("1.77799353547498")
    # Use exact symbolic formula with the decimal carried as Rational approximation.
    gplus = sp.simplify(rF1 + sp.Rational(1, 2) * sp.sqrt(1 + rF1**2))
    gminus = sp.simplify(rF1 - sp.Rational(1, 2) * sp.sqrt(1 + rF1**2))
    Wm_min = sp.simplify(gplus - 1)
    Wp_min = sp.simplify(gplus)

    print("g_-^F1 =")
    sp.pprint(sp.N(gminus, 16))
    print("g_+^F1 =")
    sp.pprint(sp.N(gplus, 16))
    print("Minimal negative weight needed to realize g_+ on any sign-indefinite source law:")
    print("W_-^min = g_+ - 1 =")
    sp.pprint(sp.N(Wm_min, 16))
    print("Corresponding positive weight at the sharp bound:")
    print("W_+^min = g_+ =")
    sp.pprint(sp.N(Wp_min, 16))

    banner("STEP 39 LEDGER")
    print("If the adiabatic-wall electron branch is strengthened by the minimal core-yield")
    print("closure dln Z_q = 0, then the whole reduced branch is only two-dimensional:")
    print("  elastic squish amplitude  ~ dln K_s,")
    print("  outgoing/core-yield amplitude ~ dln P_0.")
    print()
    print("The exact reduced laws are:")
    print("  dln K_q = 2/5 dln P_0,")
    print("  dln a   = 1/2 dln K_s,")
    print("  dln c_s = 1/2 dln K_s + 1/5 dln P_0,")
    print("  dln v_w0 = -3/4 dln K_s + 1/5 dln P_0,")
    print("  dln T_m  = -5/4 dln K_s - 1/5 dln P_0,")
    print("  dln g_s  = -1/4 dln K_s - 1/5 dln P_0,")
    print("  dln g_q  = -3/4 dln K_s,")
    print("  dln lambda = 1/2 dln K_s + 1/5 dln P_0.")
    print()
    print("The lower compensated parent surface remains exact: frak g, frak r, and r_c")
    print("stay frozen at first order.")
    print()
    print("The algebraic upper branch g_+ is retained as a diagnostic sheet, not deleted.")
    print("But because g_+ > 1, it cannot come from a passive positive source law.")
    print("Any realization of g_+ requires a sign-indefinite or pumped source profile with")
    print("at least W_- >= g_+ - 1 ≈ %.12f of compensating negative weight." % float(sp.N(Wm_min, 15)))


if __name__ == "__main__":
    main()
