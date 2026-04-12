#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 249 — explicit Family-1 threshold window and the reduction to one remaining
wall-depth datum on the reference branch.
"""


def main() -> None:
    banner("STAGE 249 — FAMILY-1 THRESHOLD WINDOW")

    Pe_req, Theta_w = sp.symbols("Pe_req Theta_w", positive=True, real=True)
    kappa = sp.Rational(12321, 5)
    eta = sp.Integer(37)
    Lambda_ell = sp.Integer(37)
    alpha_r = sp.Integer(10)

    alpha = sp.sqrt(kappa)
    Delta0 = sp.simplify(
        eta * (sp.cosh(alpha) - 1)
        / (kappa * (eta * sp.cosh(alpha) + alpha * sp.sinh(alpha)))
    )
    DeltaInf = sp.simplify(
        (eta * sp.sinh(alpha) + alpha * (sp.cosh(alpha) - 1))
        / (alpha * (eta * sp.cosh(alpha) + alpha * sp.sinh(alpha)))
    )

    Upsilon_fail = sp.simplify(Pe_req / (Lambda_ell**2 * DeltaInf))
    Upsilon_suff = sp.simplify(Pe_req / (Lambda_ell**2 * Delta0))
    Xi_fail = sp.simplify(Pe_req / DeltaInf)
    Xi_suff = sp.simplify(Pe_req / Delta0)

    subbanner("I. Exact operator thresholds on the explicit branch")
    print("Delta_0(12321/5,37) ~", sp.N(Delta0, 30))
    print("Delta_inf(12321/5,37) ~", sp.N(DeltaInf, 30))
    print()
    print("Upsilon_fail / Pe_req ~", sp.N(Upsilon_fail / Pe_req, 30))
    print("Upsilon_suff / Pe_req ~", sp.N(Upsilon_suff / Pe_req, 30))
    print()
    print("Xi_fail / Pe_req ~", sp.N(Xi_fail / Pe_req, 30))
    print("Xi_suff / Pe_req ~", sp.N(Xi_suff / Pe_req, 30))

    subbanner("II. Large-alpha interpretation")
    alpha_large = sp.simplify(111 / sp.sqrt(5))
    expect_zero("alpha = 111/sqrt(5)", sp.simplify(alpha - alpha_large))
    print("alpha ~", sp.N(alpha, 20))

    subbanner("III. Reduction to one remaining wall-depth datum")
    Upsilon_w = sp.simplify(alpha_r**2 * Theta_w)
    expect_zero("Upsilon_w = 100 Theta_w", sp.simplify(Upsilon_w - 100 * Theta_w))

    Theta_fail = sp.simplify(Upsilon_fail / alpha_r**2)
    Theta_suff = sp.simplify(Upsilon_suff / alpha_r**2)

    print("Theta_fail / Pe_req ~", sp.N(Theta_fail / Pe_req, 30))
    print("Theta_suff / Pe_req ~", sp.N(Theta_suff / Pe_req, 30))

    banner("STAGE 249 LEDGER")
    print("1. On the explicit Family-1 / healing-locked reference branch, the operator scales are")
    print("      Delta_0 ~ 1.7330e-4,   Delta_inf ~ 2.0145e-2.")
    print("2. Therefore the explicit wall-loading window is")
    print("      Upsilon_fail ~ 0.03626056 Pe_req,")
    print("      Upsilon_suff ~ 4.21495342 Pe_req.")
    print("3. Equivalently, the fixed-point coupling window is")
    print("      Xi_fail ~ 49.6407091 Pe_req,")
    print("      Xi_suff ~ 5770.2712261 Pe_req.")
    print("4. Writing the Family-1 wall depth as V0 = alpha_r mu_* with alpha_r = 10 gives")
    print("      Upsilon_w = 100 Theta_w.")
    print("5. So the entire explicit reference branch is now reduced to one remaining microscopic datum:")
    print("      Theta_w.")
    print("   Its explicit thresholds are")
    print("      Theta_fail ~ 3.62605618e-4 Pe_req,")
    print("      Theta_suff ~ 4.21495342e-2 Pe_req.")


if __name__ == "__main__":
    main()
