#!/usr/bin/env python3
"""
Step 38 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Imposes the user's adiabatic-wall ground-state closure on the lower
   compensated Family-1 branch:
      - wall entropy/enthalpy weighting frozen,
      - wall density frozen,
      - therefore d ln Theta_w = 0.
2. Solves the exact lower-branch transport laws under that closure.
3. Pushes the result through the parent-action transport formulas for
   (g_s, g_q, lambda).
4. Proves that the exact lower compensated parent surface stays tangent:
      d ln(g_q K_s / g_s lambda) = 0,
      d ln(K_s K_q / lambda^2) = 0,
      hence delta_perp = 0.

Interpretation
--------------
With the adiabatic-wall ground-state choice, the electron-track branch no longer
fights the lower compensated canonical outlet family. The remaining live bundle
freedom reduces to three reduced observables:
    (K_s, K_q, P_0).
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
    banner("STEP 38 — ADIABATIC-WALL GROUND-STATE CLOSURE")

    dKs, dKq, dP0 = sp.symbols("dln_Ks dln_Kq dln_P0", real=True)
    dTheta = sp.Integer(0)

    drho = sp.symbols("dln_rho_w", real=True)
    da = sp.symbols("dln_a", real=True)
    dcs = sp.symbols("dln_c_s", real=True)
    dZq = sp.symbols("dln_Z_q", real=True)
    dv = sp.symbols("dln_vw0", real=True)
    dT = sp.symbols("dln_Tm", real=True)
    dgs = sp.symbols("dln_gs", real=True)
    dgq = sp.symbols("dln_gq", real=True)
    dlmb = sp.symbols("dln_lambda", real=True)
    drc = sp.symbols("dln_rc", real=True)
    dfrakg = sp.symbols("dln_frakg", real=True)
    dfrakr = sp.symbols("dln_frakr", real=True)

    subbanner("XXXVIII.1 — Adiabatic-wall closure")
    print("Assume the wall entropy/enthalpy weighting and wall density are frozen:")
    print("  dln lambda_mu = 0,  dln rho_w = 0")
    print("Hence the explicit wall datum is frozen as well:")
    print("  dln Theta_w = 0")

    subbanner("XXXVIII.2 — Exact lower-branch inversion under dln Theta_w = 0")
    # Step 37 frozen-wall corollary.
    drho_sol = sp.Integer(0)
    da_sol = sp.simplify(sp.Rational(1, 2) * dKs)
    dcs_sol = sp.simplify(sp.Rational(1, 2) * dKs + sp.Rational(1, 5) * dP0)
    dZq_sol = sp.simplify(dKq - sp.Rational(2, 5) * dP0)

    print("dln rho_w =")
    sp.pprint(drho_sol)
    print("dln a =")
    sp.pprint(da_sol)
    print("dln c_s =")
    sp.pprint(dcs_sol)
    print("dln Z_q =")
    sp.pprint(dZq_sol)

    subbanner("XXXVIII.3 — Exact transport of v_w0 and T_m")
    # From the Step-36/149 sum-difference laws with dTheta = 0.
    d_ratio = sp.simplify(sp.Rational(1, 2) * dKs + sp.Rational(2, 5) * dP0)
    d_prod = sp.simplify(-2 * dKs + dKq - sp.Rational(2, 5) * dP0)
    dv_sol = sp.simplify(sp.Rational(1, 2) * (d_ratio + d_prod))
    dT_sol = sp.simplify(sp.Rational(1, 2) * (d_prod - d_ratio))

    print("dln(v_w0/T_m) =")
    sp.pprint(d_ratio)
    print("dln(v_w0 T_m) =")
    sp.pprint(d_prod)
    print("dln v_w0 =")
    sp.pprint(dv_sol)
    print("dln T_m =")
    sp.pprint(dT_sol)

    expect_zero("ratio identity recovered", d_ratio - (dv_sol - dT_sol))
    expect_zero("product identity recovered", d_prod - (dv_sol + dT_sol))

    subbanner("XXXVIII.4 — Parent-action bundle transport")
    # From the exact bundle transport laws in the notes with dTheta = 0.
    dgs_sol = sp.simplify(-sp.Rational(1, 4) * dKs + sp.Rational(1, 2) * dKq - sp.Rational(2, 5) * dP0)
    dgq_sol = sp.simplify(-sp.Rational(3, 4) * dKs + dKq - sp.Rational(2, 5) * dP0)
    dlmb_sol = sp.simplify(sp.Rational(1, 2) * dKs + sp.Rational(1, 2) * dKq)

    print("dln g_s =")
    sp.pprint(dgs_sol)
    print("dln g_q =")
    sp.pprint(dgq_sol)
    print("dln lambda =")
    sp.pprint(dlmb_sol)

    # Cross-check directly from parent definitions g_s = T_m J_s, J_s ~ a^2 ell,
    # g_q ~ Z_q L_W^{-3/2}, lambda ~ v_w0 a^2 ell L_W^{1/2}, with ell and rho_w fixed.
    dJs_direct = sp.simplify(2 * da_sol)
    dgq_direct = sp.simplify(dZq_sol - sp.Rational(3, 2) * da_sol)
    dgs_direct = sp.simplify(dT_sol + dJs_direct)
    dlmb_direct = sp.simplify(dv_sol + 2 * da_sol + sp.Rational(1, 2) * da_sol)

    expect_zero("g_s direct = bundle law", dgs_direct - dgs_sol)
    expect_zero("g_q direct = bundle law", dgq_direct - dgq_sol)
    expect_zero("lambda direct = bundle law", dlmb_direct - dlmb_sol)

    subbanner("XXXVIII.5 — Exact preservation of the lower compensated parent family")
    mouth_imbalance = sp.simplify(dgq_sol + dKs - dgs_sol - dlmb_sol)
    stiff_imbalance = sp.simplify(dKs + dKq - 2 * dlmb_sol)
    dfrakg_sol = sp.simplify(dgq_sol - dgs_sol + sp.Rational(1, 2) * dKs - sp.Rational(1, 2) * dKq)
    dfrakr_sol = sp.simplify(dlmb_sol - sp.Rational(1, 2) * dKs - sp.Rational(1, 2) * dKq)
    drc_sol = sp.simplify(2 * dlmb_sol - dKs - dKq)

    print("dln(g_q K_s / g_s lambda) =")
    sp.pprint(mouth_imbalance)
    print("dln(K_s K_q / lambda^2) =")
    sp.pprint(stiff_imbalance)
    print("dln frak g =")
    sp.pprint(dfrakg_sol)
    print("dln frak r =")
    sp.pprint(dfrakr_sol)
    print("dln r_c =")
    sp.pprint(drc_sol)

    expect_zero("mouth-coupling imbalance vanishes", mouth_imbalance)
    expect_zero("stiffness/hybridization imbalance vanishes", stiff_imbalance)
    expect_zero("frak g is frozen", dfrakg_sol)
    expect_zero("frak r is frozen", dfrakr_sol)
    expect_zero("r_c is frozen", drc_sol)

    banner("STEP 38 LEDGER")
    print("Under the adiabatic-wall ground-state closure")
    print("  dln lambda_mu = 0,  dln rho_w = 0,  hence dln Theta_w = 0,")
    print("the lower compensated Family-1 branch reduces to")
    print("  dln a   = 1/2 dln K_s,")
    print("  dln c_s = 1/2 dln K_s + 1/5 dln P_0,")
    print("  dln Z_q = dln K_q - 2/5 dln P_0,")
    print("with exact transported auxiliaries")
    print("  dln v_w0 = -3/4 dln K_s + 1/2 dln K_q,")
    print("  dln T_m  = -5/4 dln K_s + 1/2 dln K_q - 2/5 dln P_0,")
    print("  dln g_s  = -1/4 dln K_s + 1/2 dln K_q - 2/5 dln P_0,")
    print("  dln g_q  = -3/4 dln K_s + dln K_q - 2/5 dln P_0,")
    print("  dln lambda = 1/2(dln K_s + dln K_q).")
    print()
    print("Most importantly, the exact parent compensation surface is preserved identically:")
    print("  dln(g_q K_s / g_s lambda) = 0,")
    print("  dln(K_s K_q / lambda^2) = 0,")
    print("  dln frak g = dln frak r = dln r_c = 0.")
    print()
    print("So the adiabatic-wall ground-state electron track does not leave the lower")
    print("compensated canonical branch at first order. The remaining live bundle")
    print("freedom is only the reduced triple (K_s, K_q, P_0).")


if __name__ == "__main__":
    main()
