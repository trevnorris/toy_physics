#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 246 — explicit branch placement map and threshold surfaces for the canonical
thin-wall/tanh parent support branch.
"""


def main() -> None:
    banner("STAGE 246 — EXPLICIT BRANCH THRESHOLD SURFACES")

    chi_s, Lambda_ell = sp.symbols("chi_s Lambda_ell", positive=True, real=True)
    Upsilon_w, Pe_req = sp.symbols("Upsilon_w Pe_req", positive=True, real=True)

    kappa = sp.simplify(4 * chi_s**2 + sp.Rational(4, 5) * Lambda_ell**2)
    eta = sp.simplify(Lambda_ell)
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

    subbanner("I. Exact explicit-branch threshold surfaces")
    print("kappa =")
    sp.pprint(kappa)
    print()
    print("eta =")
    sp.pprint(eta)
    print()
    print("Upsilon_fail =")
    sp.pprint(Upsilon_fail)
    print()
    print("Upsilon_suff =")
    sp.pprint(Upsilon_suff)

    subbanner("II. The exact explicit-branch theorem")
    print("The first explicit branch satisfies:")
    print("  Upsilon_w <= Upsilon_fail  -> fail")
    print("  Upsilon_w >= Upsilon_suff  -> succeed")
    print("  only the intermediate band needs the full fixed-point solve.")

    subbanner("III. Shell-gradient dominated asymptotics")
    Lam = sp.symbols("Lam", positive=True, real=True)
    kappa_grad = sp.simplify(sp.Rational(4, 5) * Lam**2)
    eta_grad = Lam
    alpha_grad = sp.sqrt(kappa_grad)
    D0_grad = sp.simplify(
        eta_grad * (sp.cosh(alpha_grad) - 1)
        / (kappa_grad * (eta_grad * sp.cosh(alpha_grad) + alpha_grad * sp.sinh(alpha_grad)))
    )
    Dinf_grad = sp.simplify(
        (eta_grad * sp.sinh(alpha_grad) + alpha_grad * (sp.cosh(alpha_grad) - 1))
        / (alpha_grad * (eta_grad * sp.cosh(alpha_grad) + alpha_grad * sp.sinh(alpha_grad)))
    )
    Ufail_grad = sp.simplify(Pe_req / (Lam**2 * Dinf_grad))
    Usuff_grad = sp.simplify(Pe_req / (Lam**2 * D0_grad))

    expect_zero(
        "gradient-dominated fail asymptotic",
        sp.simplify(sp.limit(Lam * Ufail_grad, Lam, sp.oo) - 2 * Pe_req / sp.sqrt(5)),
    )
    expect_zero(
        "gradient-dominated suff asymptotic",
        sp.simplify(sp.limit(Usuff_grad, Lam, sp.oo) - sp.Rational(4, 5) * (1 + 2 / sp.sqrt(5)) * Pe_req),
    )

    subbanner("IV. Compression dominated asymptotics")
    chi = sp.symbols("chi", positive=True, real=True)
    Lam2 = sp.symbols("Lam2", positive=True, real=True)
    kappa_comp = sp.simplify(4 * chi**2)
    eta_comp = Lam2
    alpha_comp = sp.sqrt(kappa_comp)
    D0_comp = sp.simplify(
        eta_comp * (sp.cosh(alpha_comp) - 1)
        / (kappa_comp * (eta_comp * sp.cosh(alpha_comp) + alpha_comp * sp.sinh(alpha_comp)))
    )
    Dinf_comp = sp.simplify(
        (eta_comp * sp.sinh(alpha_comp) + alpha_comp * (sp.cosh(alpha_comp) - 1))
        / (alpha_comp * (eta_comp * sp.cosh(alpha_comp) + alpha_comp * sp.sinh(alpha_comp)))
    )
    Ufail_comp = sp.simplify(Pe_req / (Lam2**2 * Dinf_comp))
    Usuff_comp = sp.simplify(Pe_req / (Lam2**2 * D0_comp))

    expect_zero(
        "compression-dominated fail asymptotic",
        sp.simplify(sp.limit(Ufail_comp / (2 * Pe_req * chi / Lam2**2), chi, sp.oo) - 1),
    )
    expect_zero(
        "compression-dominated suff asymptotic",
        sp.simplify(
            sp.limit(
                Usuff_comp / (4 * Pe_req * chi**2 * (Lam2 + 2 * chi) / Lam2**3),
                chi,
                sp.oo,
            )
            - 1
        ),
    )
    expect_zero(
        "strong-compression suff asymptotic",
        sp.simplify(sp.limit(Usuff_comp / (8 * Pe_req * chi**3 / Lam2**3), chi, sp.oo) - 1),
    )

    banner("STAGE 246 LEDGER")
    print("1. The canonical explicit parent branch is placed by")
    print("      kappa = 4 chi_s^2 + (4/5) Lambda_ell^2,   eta = Lambda_ell,   W_wall = Upsilon_w Lambda_ell^2.")
    print("2. Its exact support theorem is")
    print("      Upsilon_w <= Upsilon_fail(chi_s,Lambda_ell)  -> fail,")
    print("      Upsilon_w >= Upsilon_suff(chi_s,Lambda_ell)  -> succeed.")
    print("3. In the shell-gradient dominated regime,")
    print("      Upsilon_fail ~ 2 Pe_req / (sqrt(5) Lambda_ell),")
    print("      Upsilon_suff -> (4/5)(1 + 2/sqrt(5)) Pe_req.")
    print("4. In the compression dominated regime,")
    print("      Upsilon_fail ~ 2 Pe_req chi_s / Lambda_ell^2,")
    print("      Upsilon_suff ~ 4 Pe_req chi_s^2 (Lambda_ell + 2 chi_s) / Lambda_ell^3.")
    print("5. So the next PDE-side task is now a direct branch-placement problem in (chi_s, Lambda_ell, Upsilon_w).")


if __name__ == "__main__":
    main()
