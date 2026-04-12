#!/usr/bin/env python3
"""
Step 32 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Tests the exact finite-f anomaly surface from Step 31 against the first
   explicit isotropic outlet classes in the moving-throat notes:
       - pure Robin core,
       - standalone mixed side-channel pole,
       - hybrid Robin–mixed outlet.
2. Proves that:
       - a pure Robin core can reproduce the exact electron outgoing defect, but
         it necessarily deforms the conservative even l=2 fingerprint;
       - a standalone mixed side-channel pole cannot preserve the canonical even
         branch unless it is absent;
       - a nontrivial hybrid Robin–mixed compensation branch exists and is the
         first explicit outlet class compatible with the fixed even branch.
3. Solves the exact electron-anomaly law on that compensated branch.

Interpretation
--------------
After this step the outlet problem is no longer “some deformed isotropic DtN
branch.” The first explicit audit says the only serious surviving candidate is a
compensated Robin–mixed outlet family. Pure Robin is too blunt, and a naive
standalone mixed pole is too rigid.
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
    banner("STEP 32 — EXPLICIT OUTLET-CLASS AUDIT OF THE EXACT g-2 DtN SURFACE")

    x = sp.symbols("x", positive=True, real=True)  # x = Lambda_1 f
    chi_e = sp.simplify(1 / (1 + x))

    subbanner("XXXII.1 — Pure isotropic Robin outlet")

    rho_R = sp.symbols("rho_R", real=True)
    chi_R = sp.simplify(3 / (3 - rho_R))
    rho_R_sol = sp.simplify(sp.solve(sp.Eq(chi_R, chi_e), rho_R)[0])

    print("Pure Robin outlet normalization:")
    sp.pprint(sp.Eq(sp.Symbol("chi_Q^R"), chi_R))
    print("Exact Robin value that matches the electron outgoing defect:")
    sp.pprint(sp.Eq(rho_R, rho_R_sol))
    expect_zero("pure Robin exact anomaly match", sp.simplify(chi_R.subs(rho_R, rho_R_sol) - chi_e))

    subbanner("XXXII.2 — Standalone mixed side-channel pole")

    sigma_W, kappa_W, gamma_W = sp.symbols("sigma_W kappa_W gamma_W", real=True)

    L0_mix = -(3 + sigma_W)
    L2_mix = sp.Rational(1, 3) - sigma_W * kappa_W
    L4_mix = sp.Rational(1, 9) - sigma_W * kappa_W**2
    L5_mix = sp.Rational(1, 9) - sigma_W * gamma_W

    even1_mix = sp.factor(-L2_mix / L0_mix - sp.Rational(1, 9))
    even2_mix = sp.factor(L2_mix**2 / L0_mix**2 - L4_mix / L0_mix - sp.Rational(4, 81))
    chi_mix = sp.simplify((-L5_mix / L0_mix) / (sp.Rational(1, 27)))

    print("Standalone mixed-pole even constraints:")
    print("first even condition =")
    sp.pprint(even1_mix)
    print("second even condition =")
    sp.pprint(even2_mix)

    expect_zero(
        "second even condition after imposing kappa_W = -1/9",
        sp.simplify(even2_mix.subs(kappa_W, -sp.Rational(1, 9)) + 4 * sigma_W / (81 * (sigma_W + 3))),
    )

    print("Standalone mixed-pole normalization:")
    sp.pprint(sp.Eq(sp.Symbol("chi_Q^mix"), sp.factor(chi_mix)))

    subbanner("XXXII.3 — Hybrid Robin–mixed outlet and exact canonical-even branches")

    rho_R, sigma_W, kappa_W, gamma_W = sp.symbols("rho_R sigma_W kappa_W gamma_W", real=True)

    L0_h = -3 + rho_R - sigma_W
    L2_h = sp.Rational(1, 3) - sigma_W * kappa_W
    L4_h = sp.Rational(1, 9) - sigma_W * kappa_W**2
    L5_h = sp.Rational(1, 9) - sigma_W * gamma_W

    even1_h = sp.factor(-L2_h / L0_h - sp.Rational(1, 9))
    even2_h = sp.factor(L2_h**2 / L0_h**2 - L4_h / L0_h - sp.Rational(4, 81))

    branches = sp.solve(
        [sp.Eq(sp.factor(sp.together(even1_h).as_numer_denom()[0]), 0),
         sp.Eq(sp.factor(sp.together(even2_h).as_numer_denom()[0]), 0)],
        [rho_R, kappa_W],
        dict=True,
    )

    # Discard the singular branch rho_R = sigma_W + 3.
    kept_branches = [b for b in branches if sp.simplify(b[rho_R] - (sigma_W + 3)) != 0]

    print("Canonical-even hybrid branches:")
    for br in kept_branches:
        sp.pprint(br)

    nontrivial_branch = {rho_R: 4 * sigma_W, kappa_W: sp.Rational(1, 3)}
    chi_hyb = sp.simplify(((-L5_h / L0_h) / (sp.Rational(1, 27))).subs(nontrivial_branch))
    print("Nontrivial compensated-branch normalization:")
    sp.pprint(sp.Eq(sp.Symbol("chi_Q^hyb"), sp.factor(chi_hyb)))

    gamma_can = sp.simplify(sp.solve(sp.Eq(chi_hyb, 1), gamma_W)[0])
    print("Canonical-outgoing value on the compensated branch:")
    sp.pprint(sp.Eq(gamma_W, gamma_can))
    expect_zero("canonical compensated branch", sp.simplify(chi_hyb.subs(gamma_W, gamma_can) - 1))

    subbanner("XXXII.4 — Exact electron-anomaly law on the compensated branch")

    gamma_e = sp.simplify(sp.solve(sp.Eq(chi_hyb, chi_e), gamma_W)[0])
    print("Exact compensated-branch anomaly family:")
    sp.pprint(sp.Eq(gamma_W, gamma_e))

    # Useful equivalent form relative to the canonical value 1/9.
    gamma_shift = sp.simplify(gamma_e - sp.Rational(1, 9))
    print("Shift away from the canonical odd outlet:")
    sp.pprint(sp.Eq(sp.Symbol(r"\delta\gamma_W"), gamma_shift))

    expect_zero("exact compensated branch matches the electron target", sp.simplify(chi_hyb.subs(gamma_W, gamma_e) - chi_e))

    subbanner("XXXII.5 — Numerical values on the carried benchmark")

    Lambda1_num = sp.Float("0.279605891931464")
    f_num = sp.Float("0.001161409732093")
    x_num = sp.N(Lambda1_num * f_num, 25)

    rho_R_num = sp.N(rho_R_sol.subs(x, x_num), 25)
    gamma_e_sigma_half = sp.N(gamma_e.subs({x: x_num, sigma_W: sp.Rational(1, 2)}), 25)
    gamma_shift_sigma_half = sp.N(gamma_shift.subs({x: x_num, sigma_W: sp.Rational(1, 2)}), 25)

    print("Pure Robin exact match:")
    sp.pprint(rho_R_num)
    print("Example compensated-branch odd coefficient at sigma_W = 1/2:")
    sp.pprint(gamma_e_sigma_half)
    print("Its shift away from 1/9:")
    sp.pprint(gamma_shift_sigma_half)

    banner("STEP 32 LEDGER")
    print("Pure Robin core:")
    print("  chi_Q^R = 3/(3 - rho_R)")
    print("  exact electron match -> rho_R = -3*Lambda_1*f")
    print("  but the raw Robin outlet also deforms the even l=2 fingerprint.")
    print()
    print("Standalone mixed side-channel pole:")
    print("  cannot preserve the canonical even branch unless sigma_W = 0.")
    print()
    print("Hybrid Robin–mixed outlet:")
    print("  exact canonical-even branches are")
    print("    (i)  rho_R = sigma_W,  kappa_W = 0      [trivial cancellation]")
    print("    (ii) rho_R = 4 sigma_W, kappa_W = 1/3   [nontrivial compensated branch]")
    print()
    print("On the nontrivial compensated branch:")
    print("  chi_Q^hyb = (1 - 9*sigma_W*gamma_W)/(1 - sigma_W)")
    print("  canonical outgoing branch -> gamma_W = 1/9")
    print("  electron anomaly branch   -> gamma_W = (sigma_W + Lambda_1*f)/(9*sigma_W*(1 + Lambda_1*f))")
    print()
    print("Interpretation:")
    print("  The first explicit outlet audit leaves only one serious candidate class")
    print("  for the exact g-2 DtN surface: a compensated Robin–mixed outlet. Pure")
    print("  Robin is too blunt, and a naive standalone mixed pole is too rigid.")


if __name__ == "__main__":
    main()
