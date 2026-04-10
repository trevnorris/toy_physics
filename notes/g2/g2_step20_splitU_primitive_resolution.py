#!/usr/bin/env python3
"""
Step 20 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Starts from the Step-19 coherent reduced balance pair
       Lambda_1  = -q_Lambda + q_Z + 2 q_chi + sigma q_eps,
       Delta_mix =  q_Lambda - (sigma/2) q_eps.
2. Resolves the split-blocking drift into the primitive coherent-kernel drifts
       q_eps = q_W - beta q_U,
   where beta = 2 / (11 + 9 delta_U).
3. Verifies the exact primitive coherent law
       Lambda_1  = -q_Lambda + q_Z + 2 q_chi + sigma q_W - sigma beta q_U,
       Delta_mix =  q_Lambda - (sigma/2) q_W + (sigma beta/2) q_U.
4. Solves that law exactly for (q_Lambda, q_W) in terms of the remaining
   primitive drifts and computes the minimum-Euclidean-norm representative in
   the primitive variable set
       x = (q_Lambda, q_Z, q_chi, q_W, q_U).
5. Verifies the exact minimum-norm identity
       q_U^min = - beta q_W^min,
   so the split-U drift is a suppressed counter-drift of the wall-blocking drift.
6. Uses the constructive-branch bound beta < 2/11 to show that, on the minimum-
   norm primitive closure, the quartic repair is dominated by q_W rather than q_U.

Interpretation
--------------
This is the first primitive microscopic answer to the question “what actually
moves?” on the coherent branch. The reduced Step-19 variable q_eps is not a new
mysterious slot: it resolves into a dominant wall-blocking drift q_W together
with a smaller, opposite-sign U-splitting counter-drift q_U. So the primitive
quartic repair is naturally concentrated in the wall-blocking / outgoing /
overlap / interference microledger, not in large axial-splitting motion.
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
    banner("STEP 20 — SPLIT-U PRIMITIVE RESOLUTION OF THE QUARTIC LAYER")

    # Primitive coherent-kernel variables.
    deltaU = sp.symbols("delta_U", positive=True, real=True)
    eps = sp.symbols("epsilon", positive=True, real=True)
    sigma = sp.simplify(2 * eps / (1 - eps))
    beta = sp.simplify(2 / (11 + 9 * deltaU))

    qLam, qZ, qchi, qW, qU = sp.symbols("q_Lambda q_Z q_chi q_W q_U", real=True)
    qeps = sp.symbols("q_eps", real=True)
    Lambda1, Delta_mix = sp.symbols("Lambda_1 Delta_mix", real=True)

    subbanner("XX.1 — Primitive split-U resolution of q_eps")

    print("beta =")
    sp.pprint(beta)
    expect_zero("beta - 2/(11+9 delta_U)", sp.simplify(beta - 2 / (11 + 9 * deltaU)))
    print("Constructive branch bound: 0 < beta < 2/11 for delta_U > 0.")

    qeps_split = sp.simplify(qW - beta * qU)
    print("q_eps =")
    sp.pprint(qeps_split)

    subbanner("XX.2 — Exact primitive coherent balance law")

    Lambda1_law = sp.simplify(-qLam + qZ + 2 * qchi + sigma * qeps_split)
    Delta_mix_law = sp.simplify(qLam - sp.Rational(1, 2) * sigma * qeps_split)

    print("Lambda_1 =")
    sp.pprint(Lambda1_law)
    print("Delta_mix =")
    sp.pprint(Delta_mix_law)

    expect_zero(
        "Primitive Lambda_1 law - [-q_Lambda + q_Z + 2 q_chi + sigma q_W - sigma beta q_U]",
        sp.simplify(Lambda1_law - (-qLam + qZ + 2 * qchi + sigma * qW - sigma * beta * qU)),
    )
    expect_zero(
        "Primitive Delta_mix law - [q_Lambda - (sigma/2) q_W + (sigma beta/2) q_U]",
        sp.simplify(Delta_mix_law - (qLam - sp.Rational(1, 2) * sigma * qW + sp.Rational(1, 2) * sigma * beta * qU)),
    )

    subbanner("XX.3 — Exact solution for q_Lambda and q_W")

    sol = sp.solve(
        [sp.Eq(Lambda1, Lambda1_law), sp.Eq(Delta_mix, Delta_mix_law)],
        [qLam, qW],
        dict=True,
    )[0]
    qLam_sol = sp.simplify(sol[qLam])
    qW_sol = sp.simplify(sol[qW])

    print("q_Lambda =")
    sp.pprint(qLam_sol)
    print("q_W =")
    sp.pprint(qW_sol)

    expect_zero(
        "Recovered primitive Lambda_1 constraint",
        sp.simplify(Lambda1_law.subs({qLam: qLam_sol, qW: qW_sol}) - Lambda1),
    )
    expect_zero(
        "Recovered primitive Delta_mix constraint",
        sp.simplify(Delta_mix_law.subs({qLam: qLam_sol, qW: qW_sol}) - Delta_mix),
    )

    print("So for any chosen (q_Z, q_chi, q_U), the exact primitive repair is")
    print("  q_Lambda = Lambda_1 + 2 Delta_mix - q_Z - 2 q_chi")
    print("  q_W      = beta q_U + 2 (Lambda_1 + Delta_mix - q_Z - 2 q_chi) / sigma")

    subbanner("XX.4 — Useful primitive special branches")

    # Pure wall-blocking branch qU = 0.
    qW_pure = sp.simplify(qW_sol.subs(qU, 0))
    print("Pure wall-blocking realization (q_U = 0):")
    sp.pprint(qW_pure)

    # Pure split-U realization qW = 0.
    qU_pure = sp.simplify(sp.solve(sp.Eq(qW_sol, 0), qU)[0])
    print("Pure split-U realization (q_W = 0) requires")
    sp.pprint(qU_pure)

    subbanner("XX.5 — Minimum-norm primitive representative")

    sigma_s, beta_s = sp.symbols("sigma beta", positive=True, real=True)
    x = sp.Matrix([qLam, qZ, qchi, qW, qU])
    A = sp.Matrix([
        [-1, 1, 2, sigma_s, -sigma_s * beta_s],
        [1, 0, 0, -sigma_s / 2, sigma_s * beta_s / 2],
    ])
    b = sp.Matrix([Lambda1, Delta_mix])
    x_min = sp.simplify(A.T * (A * A.T).inv() * b)

    print("x_min = (q_Lambda, q_Z, q_chi, q_W, q_U)^T =")
    sp.pprint(x_min)
    expect_zero("A x_min - b", sp.simplify(A * x_min - b))

    qLam_min, qZ_min, qchi_min, qW_min, qU_min = x_min
    expect_zero("q_U^min + beta q_W^min", sp.simplify(qU_min + beta_s * qW_min))

    print("Exact minimum-norm identity:")
    print("  q_U^min = - beta q_W^min")
    print("So the U-split drift is a suppressed counter-drift of the wall-blocking drift.")

    x_support = sp.simplify(x_min.subs(Delta_mix, 0))
    print("Minimum-norm support-carried-baseline representative (Delta_mix = 0) =")
    sp.pprint(x_support)

    subbanner("XX.6 — Primitive dominance consequence")

    print("Because beta = 2/(11 + 9 delta_U) and delta_U > 0 on the constructive branch,")
    print("we have 0 < beta < 2/11.")
    print("Therefore |q_U^min| = beta |q_W^min| < (2/11) |q_W^min|.")
    print("So the minimum-norm primitive quartic repair is always dominated by q_W rather")
    print("than q_U: wall-blocking motion carries the main split-blocking burden, while")
    print("the U-splitting drift is a smaller opposite-sign companion.")


if __name__ == "__main__":
    main()
