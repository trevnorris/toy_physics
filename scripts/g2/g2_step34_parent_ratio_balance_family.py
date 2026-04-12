#!/usr/bin/env python3
"""
Step 34 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Rewrites the balanced two-channel outlet core in terms of the normalized
   parent overlap ratios
       r = lambda/sqrt(K_s K_q),
       g = g_q sqrt(K_s)/(g_s sqrt(K_q)).
2. Verifies that the exact compensation surface collapses to the one-parameter
   parent family
       1 + r^2 = 4 (g - r)^2.
3. Derives the exact D/N side-tube law in the same variables.
4. Rewrites sigma_c and the anomaly-matching odd coefficient directly in parent
   ratio language.

Interpretation
--------------
After this step the surviving outlet is no longer a five-parameter core tuning.
It is a one-parameter parent compensation family plus one density/loading scalar.
The electron anomaly becomes a tiny odd detuning of that family rather than a
new structural mechanism.
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
    banner("STEP 34 — THE BALANCED OUTLET AS A ONE-PARAMETER PARENT COMPENSATION FAMILY")

    K_s, K_q, lam = sp.symbols("K_s K_q lambda", positive=True, real=True)
    g_s, g_q = sp.symbols("g_s g_q", positive=True, real=True)
    rho_c = sp.symbols("rho_c", positive=True, real=True)
    sigma_c = sp.symbols("sigma_c", positive=True, real=True)
    r, g = sp.symbols("mathfrak_r mathfrak_g", real=True)
    a, L_W = sp.symbols("a L_W", positive=True, real=True)
    gamma_c, gamma_0 = sp.symbols("gamma_c gamma_0", real=True)
    x = sp.symbols("x", positive=True, real=True)  # x = Lambda_1 f

    subbanner("XXXIV.1 — Parent-ratio rewrite of the exact core-balance theorem")

    r_expr = sp.simplify(lam / sp.sqrt(K_s * K_q))
    g_expr = sp.simplify(g_q * sp.sqrt(K_s) / (g_s * sp.sqrt(K_q)))
    rho_c_expr = sp.simplify(g_s**2 / K_s)
    sigma_c_expr = sp.simplify((K_s * g_q - lam * g_s) ** 2 / (K_s**2 * K_q * (1 + lam**2 / (K_s * K_q))))

    balance_eq = sp.expand(g_s**2 * (K_s * K_q + lam**2) - 4 * (K_s * g_q - lam * g_s) ** 2)
    balance_ratio = sp.simplify(balance_eq / (g_s**2 * K_s * K_q))
    balance_ratio_sub = sp.simplify(
        balance_ratio.subs(
            {
                lam: r * sp.sqrt(K_s * K_q),
                g_q: g * g_s * sp.sqrt(K_q) / sp.sqrt(K_s),
            }
        )
    )

    print("r =")
    sp.pprint(sp.Eq(sp.Symbol("r"), r_expr))
    print("g =")
    sp.pprint(sp.Eq(sp.Symbol("g"), g_expr))
    print("Normalized core-balance equation:")
    sp.pprint(sp.Eq(balance_ratio_sub, 0))
    expect_zero("balance law -> 1 + r^2 - 4(g-r)^2", sp.simplify(balance_ratio_sub - (1 + r**2 - 4 * (g - r) ** 2)))

    sol_g = sp.solve(sp.Eq(1 + r**2, 4 * (g - r) ** 2), g)
    print("Exact parent compensation branches g(r):")
    for sol in sol_g:
        sp.pprint(sp.Eq(g, sp.simplify(sol)))

    subbanner("XXXIV.2 — Exact sigma_c collapse and balanced-even reduction")

    sigma_ratio = sp.simplify(
        sigma_c_expr.subs(
            {
                lam: r * sp.sqrt(K_s * K_q),
                g_q: g * g_s * sp.sqrt(K_q) / sp.sqrt(K_s),
            }
        )
    )
    sigma_ratio_expected = sp.simplify((g_s**2 / K_s) * (g - r) ** 2 / (1 + r**2))
    expect_zero("sigma_c parent-ratio form", sp.simplify(sigma_ratio - sigma_ratio_expected))

    sigma_bal_minus = sp.simplify(sigma_ratio_expected.subs(g, sol_g[0]))
    sigma_bal_plus = sp.simplify(sigma_ratio_expected.subs(g, sol_g[1]))
    print("Balanced sigma_c on either compensation branch:")
    sp.pprint(sp.Eq(sp.Symbol("sigma_c"), sigma_bal_minus))
    expect_zero("balanced branches agree", sp.simplify(sigma_bal_minus - sigma_bal_plus))
    expect_zero("sigma_bal - rho_c/4", sp.simplify(sigma_bal_minus - rho_c_expr / 4))

    subbanner("XXXIV.3 — D/N side-tube selection in parent variables")

    kappa_0 = sp.symbols("kappa_0", real=True)
    kappa0_tube = sp.simplify(4 * L_W**2 / (sp.pi**2 * a**2))
    kappa0_bal = sp.simplify((1 + r**2) / 3)
    L_W_sol = sp.solve(sp.Eq(kappa0_tube, kappa0_bal), L_W)

    print("Balanced-even bare mixed coefficient:")
    sp.pprint(sp.Eq(kappa_0, kappa0_bal))
    print("Exact D/N side-tube length:")
    sp.pprint(sp.Eq(L_W, sp.simplify(L_W_sol[0])))
    expect_zero(
        "tube law reproduces kappa_0 = (1+r^2)/3",
        sp.simplify(kappa0_tube.subs(L_W, L_W_sol[0]) - kappa0_bal),
    )

    subbanner("XXXIV.4 — Electron-anomaly odd coefficient in parent-ratio variables")

    chi_target = sp.simplify(1 / (1 + x))
    chi_bal_core = sp.simplify((1 - 9 * sigma_c * gamma_c) / (1 - sigma_c))
    gamma_c_req = sp.simplify(sp.solve(sp.Eq(chi_bal_core, chi_target), gamma_c)[0])
    gamma_c_parent = sp.simplify(gamma_c_req.subs(sigma_c, rho_c / 4))
    gamma_0_parent = sp.simplify((1 + r**2) * gamma_c_parent)

    print("Target normalized outlet law:")
    sp.pprint(sp.Eq(sp.Symbol("chi_e"), chi_target))
    print("Required effective odd coefficient on the balanced parent family:")
    sp.pprint(sp.Eq(gamma_c, gamma_c_parent))
    print("Required bare odd coefficient:")
    sp.pprint(sp.Eq(gamma_0, gamma_0_parent))

    gamma_c_shift = sp.simplify(gamma_c_parent - sp.Rational(1, 9))
    gamma_0_shift = sp.simplify(gamma_0_parent - (1 + r**2) / 9)
    print("Deviation from the canonical balanced values:")
    sp.pprint(sp.Eq(sp.Symbol("delta_gamma_c"), gamma_c_shift))
    sp.pprint(sp.Eq(sp.Symbol("delta_gamma_0"), gamma_0_shift))

    expect_zero(
        "exact parent-family anomaly match",
        sp.simplify(chi_bal_core.subs({sigma_c: rho_c / 4, gamma_c: gamma_c_parent}) - chi_target),
    )

    subbanner("XXXIV.5 — Representative benchmark")

    Lambda1_num = sp.Float("0.279605891931464")
    f_num = sp.Float("0.001161409732093")
    x_num = sp.N(Lambda1_num * f_num, 25)

    gamma_c_rho1 = sp.N(gamma_c_parent.subs({rho_c: 1, x: x_num}), 25)
    gamma_c_shift_rho1 = sp.N(gamma_c_shift.subs({rho_c: 1, x: x_num}), 25)

    print("At rho_c = 1 and the carried electron x = Lambda_1 f:")
    print("gamma_c =")
    sp.pprint(gamma_c_rho1)
    print("gamma_c - 1/9 =")
    sp.pprint(gamma_c_shift_rho1)

    banner("STEP 34 LEDGER")
    print("The balanced outlet is now a one-parameter parent compensation family:")
    print("  r = lambda/sqrt(K_s K_q)")
    print("  g = g_q sqrt(K_s)/(g_s sqrt(K_q))")
    print("  1 + r^2 = 4 (g - r)^2")
    print()
    print("The D/N side-tube geometry is fixed by the same parent ratio:")
    print("  L_W = (pi a / 2) * sqrt((1 + r^2)/3)")
    print()
    print("On that balanced family,")
    print("  sigma_c = rho_c/4")
    print("with rho_c = g_s^2/K_s.")
    print()
    print("The electron anomaly requires only a tiny odd detuning:")
    print("  gamma_c = (rho_c + 4x)/(9 rho_c (1 + x))")
    print("  gamma_0 = (1 + r^2) (rho_c + 4x)/(9 rho_c (1 + x))")
    print("where x = Lambda_1 f.")


if __name__ == "__main__":
    main()
