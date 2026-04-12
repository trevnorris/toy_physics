#!/usr/bin/env python3
"""
Step 33 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Starts from the concrete two-channel throat-core outlet model in the
   moving-throat notes:
       - static shell compliance s,
       - mixed side-channel q,
       - hybridization lambda.
2. Re-derives the exact reduced outlet data
       rho_c, sigma_c, kappa_c, gamma_c
   and the exact compensation surface that preserves the canonical even branch.
3. Imposes the exact electron-anomaly target from Step 31/32 on that same
   compensated-even core and solves for the required odd mixed-channel
   coefficient.
4. Rewrites the answer in terms of the concrete core couplings and the D/N
   auxiliary-tube realization.

Interpretation
--------------
After this step the surviving compensated Robin–mixed outlet is no longer just a
reduced DtN class. The electron-anomaly sliver becomes one tiny odd detuning of
an otherwise canonical balanced two-channel throat core.
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
    banner("STEP 33 — THE ELECTRON ANOMALY AS A TINY ODD DETUNING OF A BALANCED CORE")

    K_s, K_q, lam = sp.symbols("K_s K_q lambda", positive=True, real=True)
    g_s, g_q = sp.symbols("g_s g_q", real=True)
    kappa0, gamma0 = sp.symbols("kappa_0 gamma_0", real=True)
    x = sp.symbols("x", positive=True, real=True)  # x = Lambda_1 f
    a, L_W = sp.symbols("a L_W", positive=True, real=True)

    subbanner("XXXIII.1 — Exact reduced outlet data of the concrete two-channel core")

    r_c = sp.simplify(lam**2 / (K_s * K_q))
    rho_c = sp.simplify(g_s**2 / K_s)
    sigma_c = sp.simplify((K_s * g_q - lam * g_s)**2 / (K_s**2 * K_q * (1 + r_c)))
    kappa_c = sp.simplify(kappa0 / (1 + r_c))
    gamma_c_expr = sp.simplify(gamma0 / (1 + r_c))

    print("r_c =")
    sp.pprint(r_c)
    print("rho_c =")
    sp.pprint(rho_c)
    print("sigma_c =")
    sp.pprint(sigma_c)
    print("kappa_c =")
    sp.pprint(kappa_c)
    print("gamma_c =")
    sp.pprint(gamma_c_expr)

    subbanner("XXXIII.2 — Exact compensation surface preserving the canonical even branch")

    balance_eq = sp.expand(g_s**2 * (K_s * K_q + lam**2) - 4 * (K_s * g_q - lam * g_s)**2)
    sol_gq = sp.solve(sp.Eq(balance_eq, 0), g_q)

    print("Coupling-balance equation:")
    sp.pprint(sp.Eq(balance_eq, 0))
    print("Exact mixed-coupling branches:")
    for sol in sol_gq:
        sp.pprint(sp.Eq(g_q, sp.simplify(sol)))

    sigma_bal = sp.simplify(sigma_c.subs(g_q, sol_gq[0]))
    expect_zero("sigma_c is the same on both balance branches", sp.simplify(sigma_bal - sigma_c.subs(g_q, sol_gq[1])))
    print("On the balance surface:")
    sp.pprint(sp.Eq(sp.Symbol("sigma_*"), sigma_bal))
    expect_zero("rho_c - 4*sigma_*", sp.simplify(rho_c - 4 * sigma_bal))

    # Even-preserving condition from Stage 98/99.
    kappa0_even = sp.simplify((1 + r_c) / 3)
    print("Exact even-preserving bare mixed-channel coefficient:")
    sp.pprint(sp.Eq(kappa0, kappa0_even))
    expect_zero("kappa_c - 1/3 under the even-preserving law", sp.simplify(kappa_c.subs(kappa0, kappa0_even) - sp.Rational(1, 3)))

    subbanner("XXXIII.3 — Exact electron-anomaly law on the balanced core")

    sigma_eff, gamma_eff = sp.symbols("sigma_eff gamma_eff", positive=True, real=True)
    chi_target = sp.simplify(1 / (1 + x))
    chi_bal_eff = sp.simplify((1 - 9 * sigma_eff * gamma_eff) / (1 - sigma_eff))
    gamma_eff_req = sp.simplify(sp.solve(sp.Eq(chi_bal_eff, chi_target), gamma_eff)[0])

    gamma_c_req = sp.simplify(gamma_eff_req.subs(sigma_eff, sigma_bal))
    gamma0_req = sp.simplify((1 + r_c) * gamma_c_req)

    print("Balanced-core normalization law:")
    sp.pprint(sp.Eq(sp.Symbol("chi_Q^core"), chi_bal_eff.subs(sigma_eff, sigma_bal).subs(gamma_eff, sp.Symbol("gamma_c"))))
    print("Required effective odd coefficient on the balanced core:")
    sp.pprint(sp.Eq(sp.Symbol("gamma_c"), gamma_c_req))
    print("Required bare odd coefficient:")
    sp.pprint(sp.Eq(gamma0, gamma0_req))

    gamma_c_shift = sp.simplify(gamma_c_req - sp.Rational(1, 9))
    gamma0_shift = sp.simplify(gamma0_req - (1 + r_c) / 9)

    print("Shift away from the canonical effective odd coefficient 1/9:")
    sp.pprint(sp.Eq(sp.Symbol(r"\delta\gamma_c"), gamma_c_shift))
    print("Shift away from the canonical bare odd coefficient (1+r_c)/9:")
    sp.pprint(sp.Eq(sp.Symbol(r"\delta\gamma_0"), gamma0_shift))

    expect_zero(
        "exact balanced-core anomaly match",
        sp.simplify(chi_bal_eff.subs({sigma_eff: sigma_bal, gamma_eff: gamma_c_req}) - chi_target),
    )

    subbanner("XXXIII.4 — Auxiliary D/N tube realization of the even-preserving branch")

    kappa0_tube = sp.simplify(4 * L_W**2 / (sp.pi**2 * a**2))
    L_W_even = sp.simplify(sp.solve(sp.Eq(kappa0_tube, kappa0_even), L_W)[0])

    print("Bare D/N mixed-tube even coefficient:")
    sp.pprint(sp.Eq(kappa0, kappa0_tube))
    print("Even-preserving auxiliary-tube length:")
    sp.pprint(sp.Eq(L_W, L_W_even))
    expect_zero(
        "tube length indeed gives kappa_0 = (1+r_c)/3",
        sp.simplify(kappa0_tube.subs(L_W, L_W_even) - kappa0_even),
    )

    subbanner("XXXIII.5 — Numerical example on the carried benchmark")

    Lambda1_num = sp.Float("0.279605891931464")
    f_num = sp.Float("0.001161409732093")
    x_num = sp.N(Lambda1_num * f_num, 25)

    gamma_c_sigma_half = sp.N(gamma_eff_req.subs({x: x_num, sigma_eff: sp.Rational(1, 2)}), 25)
    gamma_c_shift_sigma_half = sp.N((gamma_eff_req - sp.Rational(1, 9)).subs({x: x_num, sigma_eff: sp.Rational(1, 2)}), 25)

    print("At a representative balanced loading sigma_* = 1/2:")
    print("effective odd coefficient gamma_c =")
    sp.pprint(gamma_c_sigma_half)
    print("shift away from 1/9 =")
    sp.pprint(gamma_c_shift_sigma_half)

    banner("STEP 33 LEDGER")
    print("Concrete two-channel core data:")
    print("  rho_c   = g_s^2 / K_s")
    print("  sigma_c = (K_s g_q - lambda g_s)^2 / [K_s^2 K_q (1 + r_c)]")
    print("  kappa_c = kappa_0 / (1 + r_c)")
    print("  gamma_c = gamma_0 / (1 + r_c)")
    print("  r_c     = lambda^2 / (K_s K_q)")
    print()
    print("Exact balance surface preserving the canonical even branch:")
    print("  g_s^2 (K_s K_q + lambda^2) = 4 (K_s g_q - lambda g_s)^2")
    print("  kappa_0 = (1 + r_c)/3")
    print()
    print("On that balanced-even core, the electron anomaly requires")
    print("  gamma_c = (sigma_* + Lambda_1*f) / (9*sigma_*(1 + Lambda_1*f))")
    print("and hence")
    print("  gamma_0 = (1 + r_c)*(sigma_* + Lambda_1*f) / (9*sigma_*(1 + Lambda_1*f))")
    print()
    print("So the quartic g-2 sliver is now a tiny odd detuning of an otherwise")
    print("canonical balanced core. The even-preserving D/N auxiliary-tube length")
    print("relation is unchanged:")
    print("  L_W = (pi*a/2) * sqrt((1 + r_c)/3)")


if __name__ == "__main__":
    main()
