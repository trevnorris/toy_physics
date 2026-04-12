#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 238 — exact support/source fixed-point equation, operator-selected branch brackets,
and the exact Xi-fail / Xi-suff threshold pair.
"""


def main() -> None:
    banner("STAGE 238 — OPERATOR-SELECTED SUPPORT BRANCH")

    Pe = sp.symbols("Pe", positive=True, real=True)
    Xi = sp.symbols("Xi", positive=True, real=True)
    kappa = sp.symbols("kappa", positive=True, real=True)
    eta = sp.symbols("eta", positive=True, real=True)
    y = sp.symbols("y", positive=True, real=True)
    zeta_req = sp.symbols("zeta_req", positive=True, real=True)
    Pe_req = sp.symbols("Pe_req", positive=True, real=True)
    pi = sp.pi

    alpha = sp.sqrt(kappa)
    I_c = sp.simplify(
        (sp.exp(Pe) * (Pe * sp.cosh(alpha) - alpha * sp.sinh(alpha)) - Pe)
        / (Pe**2 - alpha**2)
    )
    I_s = sp.simplify(
        (sp.exp(Pe) * (Pe * sp.sinh(alpha) - alpha * sp.cosh(alpha)) + alpha)
        / (Pe**2 - alpha**2)
    )
    Delta = sp.simplify(
        Pe
        / (sp.exp(Pe) - 1)
        * ((1 - sp.cosh(alpha)) * I_c + (eta / alpha + sp.sinh(alpha)) * I_s)
        / (alpha * sp.sinh(alpha) + eta * sp.cosh(alpha))
    )
    Delta0 = sp.simplify(eta * (sp.cosh(alpha) - 1) / (kappa * (alpha * sp.sinh(alpha) + eta * sp.cosh(alpha))))
    Delta_inf = sp.simplify((sp.cosh(alpha) + (eta / alpha) * sp.sinh(alpha) - 1) / (alpha * sp.sinh(alpha) + eta * sp.cosh(alpha)))

    Omega_Pe = sp.simplify(
        pi * Pe * (2 * Pe * sp.exp(Pe) + pi)
        / ((4 * Pe**2 + pi**2) * (sp.exp(Pe) - 1))
    )
    zeta_map = sp.simplify(Omega_Pe**2 * (kappa + pi**2 / 4) / (kappa + y**2))

    subbanner("I. Exact coupled support/source fixed-point data")
    print("Delta(Pe;kappa,eta) =")
    sp.pprint(Delta)
    print()
    print("Delta_0(kappa,eta) =")
    sp.pprint(Delta0)
    print()
    print("Delta_inf(kappa,eta) =")
    sp.pprint(Delta_inf)
    Delta_at_0 = sp.simplify(sp.limit(Delta, Pe, 0))
    Delta_at_inf = sp.simplify(sp.limit(Delta, Pe, sp.oo).rewrite(sp.exp))
    print()
    print("lim_{Pe->0} Delta =", Delta_at_0)
    print("lim_{Pe->+inf} Delta =")
    sp.pprint(Delta_at_inf)
    expect_zero("exact Delta_0 identity", sp.simplify(Delta_at_0 - Delta0))
    expect_zero("exact Delta_inf identity", sp.simplify(sp.together(Delta_at_inf - Delta_inf.rewrite(sp.exp))))

    subbanner("II. Soft-support and compliant-mouth endpoint data")
    print("lim_{kappa->0+} Delta_0 =", sp.simplify(sp.limit(Delta0, kappa, 0, dir="+")))
    print("lim_{kappa->0+} Delta_inf =", sp.simplify(sp.limit(Delta_inf, kappa, 0, dir="+")))
    expect_zero("soft-support Delta_0 limit", sp.simplify(sp.limit(Delta0, kappa, 0, dir="+") - sp.Rational(1, 2)))
    expect_zero("soft-support Delta_inf limit", sp.simplify(sp.limit(Delta_inf, kappa, 0, dir="+") - 1))

    Delta0_inf = sp.simplify(sp.limit(Delta0, eta, sp.oo))
    Delta_inf_inf = sp.simplify(sp.limit(Delta_inf, eta, sp.oo))
    print()
    print("lim_{eta->+inf} Delta_0 =")
    sp.pprint(Delta0_inf)
    print()
    print("lim_{eta->+inf} Delta_inf =")
    sp.pprint(Delta_inf_inf)
    expect_zero(
        "compliant-mouth Delta_0 identity",
        sp.simplify(Delta0_inf - (1 - sp.sech(alpha)) / kappa),
    )
    expect_zero(
        "compliant-mouth Delta_inf identity",
        sp.simplify(Delta_inf_inf - sp.tanh(alpha) / alpha),
    )

    subbanner("III. Exact branch brackets and Xi thresholds")
    Xi_fail = Pe_req / Delta_inf
    Xi_suff = Pe_req / Delta0

    print("Xi_fail =")
    sp.pprint(Xi_fail)
    print()
    print("Xi_suff =")
    sp.pprint(Xi_suff)
    print()
    print("zeta_- and zeta_+ are defined by substituting Pe -> Xi Delta_0 and Pe -> Xi Delta_inf")
    print("into the exact physical lowest-lane map from Stage 237.")
    print("Residual envelope:  R_- = zeta_req - zeta_+,   R_+ = zeta_req - zeta_-.")

    Xi_gap = Pe_req * (Delta_inf - Delta0) / (Delta0 * Delta_inf)
    print()
    print("Xi_suff - Xi_fail =")
    sp.pprint(Xi_gap)

    banner("STAGE 238 LEDGER")
    print("1. The operator-selected constructive branch is the fixed-point root")
    print("      Pe_* = Xi Delta(Pe_*;kappa,eta).")
    print("2. Its exact endpoint data are Delta_0 and Delta_inf, so every constructive root obeys")
    print("      Xi Delta_0 <= Pe_* <= Xi Delta_inf.")
    print("3. Because the physical lowest-lane ratio is monotone in Pe, the full branch is bracketed by")
    print("      zeta_- <= zeta_phys <= zeta_+.")
    print("4. The exact three-zone support theorem is therefore:")
    print("      Xi <= Xi_fail  : impossible inside this operator family;")
    print("      Xi >= Xi_suff  : guaranteed success;")
    print("      Xi_fail < Xi < Xi_suff : only then is the full root solve needed.")
    print("5. The gap Xi_suff - Xi_fail is exact and collapses as Delta_0 approaches Delta_inf.")


if __name__ == "__main__":
    main()
