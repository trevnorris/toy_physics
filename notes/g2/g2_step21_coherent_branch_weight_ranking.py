#!/usr/bin/env python3
"""
Step 21 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Starts from the coherent-branch data carried by the moving-throat notes:
       sigma = 2 epsilon_*/(1-epsilon_*),
       beta  = 2/(11+9 delta_U_*),
   with constructive-branch inequalities
       0 < epsilon_* < 1,
       0 < beta < 2/11.
2. Takes the Step-20 support-carried baseline closure (Delta_mix = 0) and derives
   the exact minimum-norm primitive representative
       x = (q_Lambda, q_Z, q_chi, q_W, q_U).
3. Rewrites the primitive weights directly as rational functions of
       (epsilon_*, beta).
4. Proves the exact identities and ordering facts
       q_chi = 2 q_Z,
       q_Z > q_Lambda,
       q_Z > q_W > |q_U|,
   on the whole constructive branch.
5. Derives the only nontrivial ranking thresholds:
       q_W = q_Lambda  at  epsilon_* = 1/(2+beta^2),
       |q_U| = q_Lambda at epsilon_* = beta/(1+beta+beta^2).
6. Splits the constructive branch into three exact ranking regimes and prints a
   few representative samples.

Interpretation
--------------
The moving-throat notes do sharpen the coherent branch enough to constrain the
primitive quartic layer strongly, but they do not yet pin down a unique numeric
pair (epsilon_*, delta_U_*). So the honest next move is a branch-value audit and
phase diagram rather than fake numerical specialization. The outcome is already
very informative: the quartic repair is always dominated first by q_chi and then
by q_Z, while q_U is never a leading primitive carrier.
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


def rank_string(weight_map: dict[str, float]) -> str:
    ordered = sorted(weight_map.items(), key=lambda kv: (-kv[1], kv[0]))
    return " > ".join(name for name, _ in ordered)


def main() -> None:
    banner("STEP 21 — COHERENT-BRANCH VALUE AUDIT AND PRIMITIVE WEIGHT RANKING")

    # Coherent-branch parameters.
    eps = sp.symbols("epsilon_star", positive=True, real=True)
    beta = sp.symbols("beta", positive=True, real=True)
    deltaU = sp.symbols("delta_U_star", positive=True, real=True)
    sigma = sp.simplify(2 * eps / (1 - eps))
    beta_from_delta = sp.simplify(2 / (11 + 9 * deltaU))

    Lambda1 = sp.symbols("Lambda_1", positive=True, real=True)
    qLam, qZ, qchi, qW, qU = sp.symbols("q_Lambda q_Z q_chi q_W q_U", real=True)

    subbanner("XXI.1 — What the coherent branch actually fixes")

    print("sigma =")
    sp.pprint(sigma)
    print("beta(delta_U_*) =")
    sp.pprint(beta_from_delta)
    print("Constructive-branch inequalities carried by the notes:")
    print("  0 < epsilon_* < 1")
    print("  0 < delta_U_*  and therefore  0 < beta < 2/11")

    subbanner("XXI.2 — Step-20 minimum-norm primitive closure at Delta_mix = 0")

    A = sp.Matrix([
        [-1, 1, 2, sigma, -sigma * beta],
        [1, 0, 0, -sigma / 2, sigma * beta / 2],
    ])
    b = sp.Matrix([Lambda1, 0])
    x_min = sp.simplify(A.T * (A * A.T).inv() * b)

    print("x_min = (q_Lambda, q_Z, q_chi, q_W, q_U)^T =")
    sp.pprint(x_min)
    expect_zero("A x_min - b", sp.simplify(A * x_min - b))

    qLam_min, qZ_min, qchi_min, qW_min, qU_min = x_min

    subbanner("XXI.3 — Exact primitive weights in (epsilon_*, beta)")

    Den = sp.simplify(5 - 10 * eps + (11 + 6 * beta**2) * eps**2)
    print("Denominator N(epsilon_*, beta) =")
    sp.pprint(Den)
    print("Useful positivity rewrite:")
    sp.pprint(sp.expand(5 * (1 - eps) ** 2 + 6 * eps**2 * (1 + beta**2)))
    expect_zero(
        "N - [5(1-epsilon_*)^2 + 6 epsilon_*^2 (1+beta^2)]",
        sp.simplify(Den - (5 * (1 - eps) ** 2 + 6 * eps**2 * (1 + beta**2))),
    )

    wLam = sp.simplify(eps**2 * (1 + beta**2) / Den)
    wZ = sp.simplify((1 - 2 * eps + (2 + beta**2) * eps**2) / Den)
    wchi = sp.simplify(2 * (1 - 2 * eps + (2 + beta**2) * eps**2) / Den)
    wW = sp.simplify(eps * (1 - eps) / Den)
    wUabs = sp.simplify(beta * eps * (1 - eps) / Den)

    print("Normalized primitive weights (q_i = Lambda_1 * w_i, except q_U = -Lambda_1 * |w_U|):")
    print("w_Lambda =")
    sp.pprint(wLam)
    print("w_Z =")
    sp.pprint(wZ)
    print("w_chi =")
    sp.pprint(wchi)
    print("w_W =")
    sp.pprint(wW)
    print("|w_U| =")
    sp.pprint(wUabs)

    expect_zero("q_Lambda^min - Lambda_1 w_Lambda", sp.simplify(qLam_min - Lambda1 * wLam))
    expect_zero("q_Z^min - Lambda_1 w_Z", sp.simplify(qZ_min - Lambda1 * wZ))
    expect_zero("q_chi^min - Lambda_1 w_chi", sp.simplify(qchi_min - Lambda1 * wchi))
    expect_zero("q_W^min - Lambda_1 w_W", sp.simplify(qW_min - Lambda1 * wW))
    expect_zero("q_U^min + Lambda_1 |w_U|", sp.simplify(qU_min + Lambda1 * wUabs))

    subbanner("XXI.4 — Exact identities and always-true ordering facts")

    expect_zero("w_chi - 2 w_Z", sp.simplify(wchi - 2 * wZ))
    expect_zero("w_Z - w_Lambda - (1-epsilon_*)^2/N", sp.simplify(wZ - wLam - (1 - eps) ** 2 / Den))
    expect_zero(
        "w_Z - w_W - [beta^2 epsilon_*^2 + 3 (epsilon_* - 1/2)^2 + 1/4]/N",
        sp.simplify(wZ - wW - (beta**2 * eps**2 + 3 * (eps - sp.Rational(1, 2)) ** 2 + sp.Rational(1, 4)) / Den),
    )
    expect_zero("w_W - |w_U| - beta-free factor", sp.simplify(wW - wUabs - eps * (1 - eps) * (1 - beta) / Den))
    expect_zero(
        "w_chi - w_W - [2 beta^2 epsilon_*^2 + 5 (epsilon_* - 1/2)^2 + 3/4]/N",
        sp.simplify(wchi - wW - (2 * beta**2 * eps**2 + 5 * (eps - sp.Rational(1, 2)) ** 2 + sp.Rational(3, 4)) / Den),
    )

    print("Consequences on the whole constructive branch:")
    print("  w_chi = 2 w_Z")
    print("  w_chi > w_Z > w_W > |w_U| > 0  is only partly true: w_Z > w_W > |w_U| always,")
    print("  and w_Z > w_Lambda always, but the ordering of w_Lambda relative to w_W and |w_U|")
    print("  depends on epsilon_*.")

    subbanner("XXI.5 — Exact ranking thresholds")

    eps_cross = sp.simplify(1 / (2 + beta**2))
    eps_small = sp.simplify(beta / (1 + beta + beta**2))

    expect_zero("w_W - w_Lambda - epsilon_* [1-(2+beta^2)epsilon_*]/N", sp.factor(sp.simplify(wW - wLam - eps * (1 - (2 + beta**2) * eps) / Den)))
    expect_zero("w_Lambda - |w_U| - epsilon_*[(1+beta+beta^2)epsilon_* - beta]/N", sp.factor(sp.simplify(wLam - wUabs - eps * ((1 + beta + beta**2) * eps - beta) / Den)))

    print("w_W = w_Lambda at")
    sp.pprint(eps_cross)
    print("|w_U| = w_Lambda at")
    sp.pprint(eps_small)

    beta_max = sp.Rational(2, 11)
    eps_small_max = sp.simplify(beta_max / (1 + beta_max + beta_max**2))
    eps_cross_min = sp.simplify(1 / (2 + beta_max**2))

    print("Using 0 < beta < 2/11, the thresholds satisfy")
    print("  0 < epsilon_small(beta) <")
    sp.pprint(eps_small_max)
    print("  and")
    print("  ")
    sp.pprint(eps_cross_min)
    print("< epsilon_cross(beta) < 1/2")

    subbanner("XXI.6 — Exact constructive-branch ranking regimes")

    print("Region I  (very weak blocking):")
    print("  0 < epsilon_* < beta/(1+beta+beta^2)")
    print("  =>  q_chi > q_Z > q_W > |q_U| > q_Lambda")
    print()
    print("Region II (intermediate blocking):")
    print("  beta/(1+beta+beta^2) < epsilon_* < 1/(2+beta^2)")
    print("  =>  q_chi > q_Z > q_W > q_Lambda > |q_U|")
    print()
    print("Region III (strong blocking):")
    print("  epsilon_* > 1/(2+beta^2)")
    print("  =>  q_chi > q_Z > q_Lambda > q_W > |q_U|")

    subbanner("XXI.7 — Limiting behavior")

    print("epsilon_* -> 0^+:")
    print("  w_Lambda ->")
    sp.pprint(sp.simplify(sp.limit(wLam, eps, 0, dir="+")))
    print("  w_Z ->")
    sp.pprint(sp.simplify(sp.limit(wZ, eps, 0, dir="+")))
    print("  w_chi ->")
    sp.pprint(sp.simplify(sp.limit(wchi, eps, 0, dir="+")))
    print("  w_W ->")
    sp.pprint(sp.simplify(sp.limit(wW, eps, 0, dir="+")))
    print("  |w_U| ->")
    sp.pprint(sp.simplify(sp.limit(wUabs, eps, 0, dir="+")))

    print("epsilon_* -> 1^-:")
    print("  w_Lambda ->")
    sp.pprint(sp.simplify(sp.limit(wLam, eps, 1, dir="-")))
    print("  w_Z ->")
    sp.pprint(sp.simplify(sp.limit(wZ, eps, 1, dir="-")))
    print("  w_chi ->")
    sp.pprint(sp.simplify(sp.limit(wchi, eps, 1, dir="-")))
    print("  w_W ->")
    sp.pprint(sp.simplify(sp.limit(wW, eps, 1, dir="-")))
    print("  |w_U| ->")
    sp.pprint(sp.simplify(sp.limit(wUabs, eps, 1, dir="-")))

    subbanner("XXI.8 — Representative constructive-branch samples")

    samples = [
        (sp.Rational(1, 2), sp.Rational(1, 20)),
        (sp.Rational(1, 2), sp.Rational(1, 5)),
        (sp.Rational(1, 2), sp.Rational(4, 5)),
        (sp.Integer(1), sp.Rational(1, 20)),
        (sp.Integer(1), sp.Rational(3, 10)),
        (sp.Integer(1), sp.Rational(4, 5)),
        (sp.Integer(5), sp.Rational(1, 20)),
        (sp.Integer(5), sp.Rational(3, 10)),
        (sp.Integer(5), sp.Rational(4, 5)),
    ]

    print("delta_U_*   epsilon_*   beta        ranking")
    for delta_val, eps_val in samples:
        beta_val = sp.N(beta_from_delta.subs(deltaU, delta_val), 12)
        vals = {
            "q_chi": float(sp.N(wchi.subs({eps: eps_val, beta: beta_val}), 20)),
            "q_Z": float(sp.N(wZ.subs({eps: eps_val, beta: beta_val}), 20)),
            "q_Lambda": float(sp.N(wLam.subs({eps: eps_val, beta: beta_val}), 20)),
            "q_W": float(sp.N(wW.subs({eps: eps_val, beta: beta_val}), 20)),
            "|q_U|": float(sp.N(wUabs.subs({eps: eps_val, beta: beta_val}), 20)),
        }
        print(f"{str(delta_val):>8}   {str(eps_val):>10}   {beta_val:>8}   {rank_string(vals)}")


if __name__ == "__main__":
    main()
