#!/usr/bin/env python3
"""
Step 31 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Starts from the exact Stage-90/91 isotropic DtN deformation law

       chi_Q = 3*(S*beta^5 + 9*Sigma_5) / (3*S - Sigma_0),

   after the canonical even fingerprint has been enforced.
2. Inserts the exact Step-29 electron-point outgoing-normalization target

       chi_e = 1 / (1 + Lambda_1*f),

   so the quartic g-2 sliver becomes one exact finite-f branch-selection surface.
3. Solves the exact anomaly surface for Sigma_5 and derives three exact pure
   realizations:
       - compensated argument deformation,
       - pure static additive core,
       - pure odd l=2 core outlet.
4. Verifies that the exact finite-f surface reduces to the linear Step-30 law
       5 b + a_0/3 + 9 a_5 = -Lambda_1
   near the canonical branch.

Interpretation
--------------
After this step the quartic anomaly layer is no longer only a tangent constraint.
It becomes one exact finite-f isotropic DtN branch-selection surface. The next
question is therefore not whether a deformation exists, but which explicit
moving-throat outlet class can realize that exact surface while keeping the
conservative even l=2 fingerprint intact.
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
    banner("STEP 31 — THE EXACT FINITE-f ISOTROPIC DtN ANOMALY SURFACE")

    S, beta, Sigma0, Sigma5 = sp.symbols("S beta Sigma_0 Sigma_5", real=True)
    Lambda1, f = sp.symbols("Lambda_1 f", positive=True, real=True)
    x = sp.symbols("x", positive=True, real=True)  # x = Lambda_1 f

    chi_exact = sp.simplify(3 * (S * beta**5 + 9 * Sigma5) / (3 * S - Sigma0))
    chi_target = sp.simplify(1 / (1 + x))

    subbanner("XXXI.1 — Exact electron-point outgoing-normalization target")

    print("Exact isotropic DtN deformation law:")
    sp.pprint(sp.Eq(sp.Symbol("chi_Q"), chi_exact))
    print("Exact electron-point target from Step 29:")
    sp.pprint(sp.Eq(sp.Symbol("chi_e"), chi_target))

    exact_surface = sp.expand((1 + x) * 3 * (S * beta**5 + 9 * Sigma5) - (3 * S - Sigma0))
    print("Exact anomaly surface:")
    sp.pprint(sp.Eq(exact_surface, 0))

    Sigma5_surface = sp.simplify(sp.solve(sp.Eq(exact_surface, 0), Sigma5)[0])
    print("Solved for Sigma_5:")
    sp.pprint(sp.Eq(Sigma5, Sigma5_surface))

    subbanner("XXXI.2 — Canonical-preservation submanifold as a special case")

    Sigma5_preserve = sp.simplify(sp.solve(sp.Eq(chi_exact, 1), Sigma5)[0])
    print("Exact preservation submanifold chi_Q = 1:")
    sp.pprint(sp.Eq(Sigma5, Sigma5_preserve))

    expect_zero(
        "canonical-preservation law recovered by setting x = 0 in the anomaly surface",
        sp.simplify(Sigma5_surface.subs(x, 0) - Sigma5_preserve),
    )

    subbanner("XXXI.3 — Exact pure isotropic realizations of the electron anomaly")

    beta_pure = sp.simplify(sp.solve(sp.Eq(chi_exact.subs({Sigma0: 0, Sigma5: 0}), chi_target), beta)[0])
    Sigma0_pure = sp.simplify(sp.solve(sp.Eq(chi_exact.subs({beta: 1, Sigma5: 0}), chi_target), Sigma0)[0])
    Sigma5_pure = sp.simplify(sp.solve(sp.Eq(chi_exact.subs({beta: 1, Sigma0: 0}), chi_target), Sigma5)[0])

    print("Pure compensated argument-deformation branch:")
    sp.pprint(sp.Eq(beta, beta_pure))
    print("Pure static additive-core branch:")
    sp.pprint(sp.Eq(Sigma0, Sigma0_pure))
    print("Pure odd l=2 core-outlet branch:")
    sp.pprint(sp.Eq(Sigma5, Sigma5_pure))

    expect_zero(
        "pure-beta branch hits the exact anomaly surface",
        sp.simplify(exact_surface.subs({Sigma0: 0, Sigma5: 0, beta: beta_pure})),
    )
    expect_zero(
        "pure-Sigma0 branch hits the exact anomaly surface",
        sp.simplify(exact_surface.subs({beta: 1, Sigma5: 0, Sigma0: Sigma0_pure})),
    )
    expect_zero(
        "pure-Sigma5 branch hits the exact anomaly surface",
        sp.simplify(exact_surface.subs({beta: 1, Sigma0: 0, Sigma5: Sigma5_pure})),
    )

    subbanner("XXXI.4 — Linearized reduction back to the Step-30 tangent law")

    eps, s, b, a0, a5 = sp.symbols("eps s b a_0 a_5", real=True)
    chi_lin = sp.expand(
        sp.series(
            chi_exact.subs({
                S: 1 + eps * s,
                beta: 1 + eps * b,
                Sigma0: eps * a0,
                Sigma5: eps * a5,
            }),
            eps,
            0,
            2,
        ).removeO()
    )
    chi_tgt_lin = sp.expand(sp.series(1 / (1 + eps * Lambda1), eps, 0, 2).removeO())
    delta_lin = sp.simplify((chi_lin - chi_tgt_lin).coeff(eps, 1))
    print("First-order mismatch coefficient around the canonical branch:")
    sp.pprint(delta_lin)
    expect_zero(
        "linearized exact surface -> Step-30 tangent law",
        sp.simplify(delta_lin - (5 * b + a0 / 3 + 9 * a5 + Lambda1)),
    )

    subbanner("XXXI.5 — Numerical values on the carried benchmark")

    Lambda1_num = sp.Float("0.279605891931464")
    f_num = sp.Float("0.001161409732093")
    x_num = sp.N(Lambda1_num * f_num, 25)

    chi_num = sp.N(chi_target.subs(x, x_num), 25)
    beta_num = sp.N(beta_pure.subs(x, x_num), 25)
    Sigma0_num = sp.N(Sigma0_pure.subs({x: x_num, S: 1}), 25)
    Sigma5_num = sp.N(Sigma5_pure.subs({x: x_num, S: 1}), 25)

    print("x = Lambda_1 f =")
    sp.pprint(x_num)
    print("chi_e =")
    sp.pprint(chi_num)
    print("Pure-beta exact branch (S arbitrary, Sigma_0 = Sigma_5 = 0):")
    sp.pprint(beta_num)
    print("Pure-static exact branch at S = 1:")
    sp.pprint(Sigma0_num)
    print("Pure-odd exact branch at S = 1:")
    sp.pprint(Sigma5_num)

    banner("STEP 31 LEDGER")
    print("Exact finite-f target:")
    print("  chi_e = 1 / (1 + Lambda_1 f)")
    print()
    print("Exact isotropic DtN anomaly surface:")
    print("  3*(S*beta^5 + 9*Sigma_5)*(1 + Lambda_1 f) = 3*S - Sigma_0")
    print("or equivalently")
    print("  3*S*((1 + Lambda_1 f)*beta^5 - 1) + Sigma_0 + 27*(1 + Lambda_1 f)*Sigma_5 = 0")
    print()
    print("Exact pure branches:")
    print("  pure beta   : beta   = (1 + Lambda_1 f)^(-1/5)")
    print("  pure Sigma0 : Sigma_0 = -3*S*Lambda_1*f")
    print("  pure Sigma5 : Sigma_5 = -S*Lambda_1*f / (9*(1 + Lambda_1 f))")
    print()
    print("Carried-benchmark values:")
    print(f"  chi_e        = {chi_num}")
    print(f"  pure beta    = {beta_num}")
    print(f"  pure Sigma_0 = {Sigma0_num}   (at S = 1)")
    print(f"  pure Sigma_5 = {Sigma5_num}   (at S = 1)")
    print()
    print("Interpretation:")
    print("  The quartic g-2 sliver is now an exact finite-f isotropic DtN surface, not")
    print("  only a tangent. The next task is to determine which explicit outlet class")
    print("  can realize this surface without spoiling the already-fixed conservative")
    print("  even l=2 fingerprint.")


if __name__ == "__main__":
    main()
