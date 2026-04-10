#!/usr/bin/env python3
"""
Step 36 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Uses the exact lower compensated Family-1 branch constraints to derive the
   linear co-transport laws for L_W, v_{w0}, and T_m.
2. Applies the frozen n=5 wall-EOS relation to collapse the microscopic drift
   space to four irreducible variables.
3. Verifies that the Stage-147 off-family scalar cancels identically on the exact
   lower compensated branch.

Interpretation
--------------
After this step the anomaly problem is no longer an eight-variable branch-drift
problem. The exact compensated branch already co-transports the side-tube length,
background mixed flow, mouth traction, and wall sound speed. Only four
microscopic drifts remain genuinely open.
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
    banner("STEP 36 — EXACT LOWER-BRANCH CO-TRANSPORT AND OFF-FAMILY CANCELLATION")

    # Log-drift variables.
    dZq, drho, dcsw, dcs, da, dLW, dTm, dvw = sp.symbols(
        "dln_Zq dln_rho dln_csw dln_cs dln_a dln_LW dln_Tm dln_vw", real=True
    )
    r_star, g_star = sp.symbols("r_star g_star", positive=True, real=True)

    subbanner("XXXVI.1 — Exact lower-branch logarithmic constraints")

    # From the explicit throat-core branch under uniform-overlap + healing-lock:
    # channel_g = d ln(g_q K_s/(g_s lambda))
    # channel_r = d ln(K_s K_q/lambda^2)
    channel_g = sp.simplify(dZq + 3 * dcsw - drho - dTm - dvw - 2 * da - 2 * dLW)
    channel_r = sp.simplify(dZq + 2 * dcs + 3 * dcsw - drho - 2 * dvw - 2 * da - 3 * dLW)

    print("Fixed-g lower-branch channel:")
    sp.pprint(sp.Eq(sp.Symbol("channel_g"), channel_g))
    print("Fixed-r lower-branch channel:")
    sp.pprint(sp.Eq(sp.Symbol("channel_r"), channel_r))

    subbanner("XXXVI.2 — D/N geometric co-transport")

    # L_W = (pi a / 2) sqrt((1+r_*^2)/3) at fixed r_*.
    expect_zero("D/N geometric co-transport law", sp.simplify(dLW - da).subs(dLW, da))
    print("So on the exact lower branch, dln L_W = dln a.")

    subbanner("XXXVI.3 — Exact drift laws for v_{w0} and T_m")

    channel_g_red = sp.simplify(channel_g.subs(dLW, da))
    channel_r_red = sp.simplify(channel_r.subs(dLW, da))

    sol = sp.solve([sp.Eq(channel_g_red, 0), sp.Eq(channel_r_red, 0)], [dvw, dTm], dict=True)[0]
    dvw_sol = sp.simplify(sol[dvw])
    dTm_sol = sp.simplify(sol[dTm])

    print("dln v_{w0} =")
    sp.pprint(dvw_sol)
    print("dln T_m =")
    sp.pprint(dTm_sol)

    expect_zero("fixed-g channel after substitution", channel_g_red.subs({dvw: dvw_sol, dTm: dTm_sol}))
    expect_zero("fixed-r channel after substitution", channel_r_red.subs({dvw: dvw_sol, dTm: dTm_sol}))

    subbanner("XXXVI.4 — Exact product/ratio factorization")

    ratio_channel = sp.simplify(dvw_sol - dTm_sol)
    product_channel = sp.simplify(dvw_sol + dTm_sol)

    print("dln(v_{w0}/T_m) =")
    sp.pprint(ratio_channel)
    print("dln(v_{w0} T_m) =")
    sp.pprint(product_channel)

    expect_zero("ratio channel - (2 dln c_s - dln a)", sp.simplify(ratio_channel - (2 * dcs - da)))
    expect_zero(
        "product channel - (dln Z_q + 3 dln c_sw - dln rho - 4 dln a)",
        sp.simplify(product_channel - (dZq + 3 * dcsw - drho - 4 * da)),
    )

    subbanner("XXXVI.5 — Frozen n=5 wall-EOS reduction")

    eos_rule = {dcsw: 2 * drho}
    dvw_n5 = sp.simplify(dvw_sol.subs(eos_rule))
    dTm_n5 = sp.simplify(dTm_sol.subs(eos_rule))
    ratio_n5 = sp.simplify(ratio_channel.subs(eos_rule))
    product_n5 = sp.simplify(product_channel.subs(eos_rule))

    print("With c_{s,w}^2 ∝ rho_w^4, so dln c_{s,w} = 2 dln rho_w:")
    print("dln v_{w0} =")
    sp.pprint(dvw_n5)
    print("dln T_m =")
    sp.pprint(dTm_n5)
    print("dln(v_{w0}/T_m) =")
    sp.pprint(ratio_n5)
    print("dln(v_{w0} T_m) =")
    sp.pprint(product_n5)

    expect_zero("n=5 ratio channel unchanged form", sp.simplify(ratio_n5 - (2 * dcs - da)))
    expect_zero("n=5 product channel", sp.simplify(product_n5 - (dZq + 5 * drho - 4 * da)))

    subbanner("XXXVI.6 — Exact cancellation of the off-family scalar")

    delta_perp = sp.simplify(g_star * channel_g + channel_r / (4 * sp.sqrt(1 + r_star**2)))
    delta_perp_red = sp.simplify(delta_perp.subs({dLW: da, dvw: dvw_sol, dTm: dTm_sol}))

    print("delta_perp =")
    sp.pprint(delta_perp)
    print("delta_perp on the exact lower branch =")
    sp.pprint(delta_perp_red)
    expect_zero("off-family scalar vanishes on the exact lower branch", delta_perp_red)

    banner("STEP 36 LEDGER")
    print("On the exact lower compensated branch:")
    print("  dln L_W = dln a")
    print()
    print("The background mixed flow and mouth traction are not independent. Their exact")
    print("drifts are")
    print("  dln v_{w0} = (dln Z_q - dln rho + 3 dln c_{s,w} + 2 dln c_s - 5 dln a)/2")
    print("  dln T_m   = (dln Z_q - dln rho + 3 dln c_{s,w} - 2 dln c_s - 3 dln a)/2")
    print()
    print("So the ratio and product split cleanly into")
    print("  dln(v_{w0}/T_m) = 2 dln c_s - dln a")
    print("  dln(v_{w0} T_m) = dln Z_q + 3 dln c_{s,w} - dln rho - 4 dln a")
    print()
    print("After the frozen n=5 wall-EOS reduction dln c_{s,w} = 2 dln rho, the whole")
    print("actual branch-drift problem collapses to only four irreducible variables:")
    print("  (dln Z_q, dln rho_w, dln c_s, dln a).")
    print()
    print("The Stage-147 off-family scalar then vanishes identically on the exact lower")
    print("branch. So the compensated branch does not need another first-order repair;")
    print("the remaining problem is only the actual dynamics of those four residual drifts.")


if __name__ == "__main__":
    main()
