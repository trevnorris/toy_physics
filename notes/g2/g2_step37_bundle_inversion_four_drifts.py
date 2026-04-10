#!/usr/bin/env python3
"""
Step 37 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Takes the four exact branch laws for Theta_w, K_s, K_q, and P_0 on the lower
   compensated Family-1 branch.
2. Solves them exactly for the remaining four irreducible microscopic drifts
   (Z_q, rho_w, c_s, a).
3. Rewrites the result in full-bundle language using P_0 = N_0 / D_0.
4. Records the frozen-wall corollary d ln Theta_w = 0.

Interpretation
--------------
After this step the last microscopic freedom is no longer diffuse. The four
residual branch drifts are exact algebraic images of the bundle observables
(Theta_w, K_s, K_q, P_0), or equivalently (Theta_w, K_s, K_q, N_0/D_0).
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
    banner("STEP 37 — EXACT BUNDLE INVERSION OF THE LAST FOUR IRREDUCIBLE DRIFTS")

    dTheta, dKs, dKq, dP0 = sp.symbols("dln_Theta dln_Ks dln_Kq dln_P0", real=True)
    dZq, drho, dcs, da = sp.symbols("dln_Zq dln_rho dln_cs dln_a", real=True)
    dN0, dD0 = sp.symbols("dln_N0 dln_D0", real=True)

    subbanner("XXXVII.1 — The four exact branch laws")

    eq1 = sp.Eq(dTheta, 2 * drho)
    eq2 = sp.Eq(dKs, 2 * da + drho)
    eq3 = sp.Eq(dKq, dZq + 2 * dcs - 2 * da)
    eq4 = sp.Eq(dP0, 5 * (dcs - da))

    print("dln Theta_w = 2 dln rho_w")
    print("dln K_s     = 2 dln a + dln rho_w")
    print("dln K_q     = dln Z_q + 2 dln c_s - 2 dln a")
    print("dln P_0     = 5 (dln c_s - dln a)")

    subbanner("XXXVII.2 — Exact inversion")

    sol = sp.solve([eq1, eq2, eq3, eq4], [drho, da, dcs, dZq], dict=True)[0]
    drho_sol = sp.simplify(sol[drho])
    da_sol = sp.simplify(sol[da])
    dcs_sol = sp.simplify(sol[dcs])
    dZq_sol = sp.simplify(sol[dZq])

    print("dln rho_w =")
    sp.pprint(drho_sol)
    print("dln a =")
    sp.pprint(da_sol)
    print("dln c_s =")
    sp.pprint(dcs_sol)
    print("dln Z_q =")
    sp.pprint(dZq_sol)

    expect_zero("Theta-law recovered", sp.simplify(dTheta - 2 * drho_sol))
    expect_zero("K_s law recovered", sp.simplify(dKs - (2 * da_sol + drho_sol)))
    expect_zero("K_q law recovered", sp.simplify(dKq - (dZq_sol + 2 * dcs_sol - 2 * da_sol)))
    expect_zero("P_0 law recovered", sp.simplify(dP0 - 5 * (dcs_sol - da_sol)))

    subbanner("XXXVII.3 — Full-bundle form using P_0 = N_0 / D_0")

    dcs_bundle = sp.simplify(dcs_sol.subs(dP0, dN0 - dD0))
    dZq_bundle = sp.simplify(dZq_sol.subs(dP0, dN0 - dD0))

    print("dln c_s in full-bundle variables =")
    sp.pprint(dcs_bundle)
    print("dln Z_q in full-bundle variables =")
    sp.pprint(dZq_bundle)
    expect_zero("bundle identity for dln c_s", sp.simplify(dcs_bundle - (sp.Rational(1,2)*dKs - sp.Rational(1,4)*dTheta + sp.Rational(1,5)*(dN0 - dD0))))
    expect_zero("bundle identity for dln Z_q", sp.simplify(dZq_bundle - (dKq - sp.Rational(2,5)*(dN0 - dD0))))

    subbanner("XXXVII.4 — Frozen-wall corollary")

    frozen = {dTheta: 0}
    drho_frozen = sp.simplify(drho_sol.subs(frozen))
    da_frozen = sp.simplify(da_sol.subs(frozen))
    dcs_frozen = sp.simplify(dcs_sol.subs(frozen))
    dZq_frozen = sp.simplify(dZq_sol.subs(frozen))

    print("If the explicit Family-1 wall datum is frozen at first order (dln Theta_w = 0):")
    print("dln rho_w =")
    sp.pprint(drho_frozen)
    print("dln a =")
    sp.pprint(da_frozen)
    print("dln c_s =")
    sp.pprint(dcs_frozen)
    print("dln Z_q =")
    sp.pprint(dZq_frozen)

    banner("STEP 37 LEDGER")
    print("The last four irreducible lower-branch drifts are exact algebraic images of")
    print("the bundle observables (Theta_w, K_s, K_q, P_0):")
    print("  dln rho_w = 1/2 dln Theta_w")
    print("  dln a     = 1/2 dln K_s - 1/4 dln Theta_w")
    print("  dln c_s   = 1/2 dln K_s - 1/4 dln Theta_w + 1/5 dln P_0")
    print("  dln Z_q   = dln K_q - 2/5 dln P_0")
    print()
    print("Since P_0 = N_0/D_0, the grouped bundle enters only through the isotropic")
    print("static normalization ratio N_0/D_0.")
    print()
    print("If the explicit Family-1 wall datum is frozen at first order, then")
    print("  dln rho_w = 0,")
    print("  dln a = 1/2 dln K_s,")
    print("  dln c_s = 1/2 dln K_s + 1/5 dln P_0,")
    print("  dln Z_q = dln K_q - 2/5 dln P_0.")


if __name__ == "__main__":
    main()
