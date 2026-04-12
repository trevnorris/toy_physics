#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 278 — the selected moving-throat outgoing branch makes chi_Q a direct continuum observable.

What this script does
---------------------
1. Starts from the exact compensated Robin–mixed outgoing branch

       chi_Q = (1 - 9 sigma_W gamma_W)/(1 - sigma_W).

2. Extracts the direct outgoing defect

       Delta_Q = chi_Q - 1 = sigma_W (1 - 9 gamma_W)/(1 - sigma_W).

3. Inverts the relation to solve for the physical odd mixed-channel coefficient
   gamma_W in terms of the observable pair (sigma_W, chi_Q).
4. Records the natural-source-map conservative normalization required by the exact
   odd factorization,

       N_Q^(req) = 1/chi_Q,

   so the actual selected outgoing branch supplies chi_Q directly and therefore
   determines the last reduced conservative normalization requirement.
5. Gives the two exact DtN-gauge observables on the compensated branch:
      - core gauge  : Sigma_0 = 3 sigma_W,  Sigma_5 = - sigma_W gamma_W
      - scale gauge : S = 1-sigma_W,        Sigma_5 = sigma_W (1/9 - gamma_W)
"""


def main() -> None:
    banner("STAGE 278 — SELECTED OUTGOING BRANCH AS DIRECT chi_Q OBSERVABLE")

    sigmaW, gammaW = sp.symbols("sigma_W gamma_W", positive=True, real=True)
    chiQ = sp.symbols("chi_Q", positive=True, real=True)

    subbanner("I. Exact outgoing normalization on the compensated branch")
    chi_hyb = sp.simplify((1 - 9 * sigmaW * gammaW) / (1 - sigmaW))
    DeltaQ = sp.simplify(chi_hyb - 1)

    print("chi_Q^(hyb) =")
    sp.pprint(chi_hyb)
    print("Delta_Q = chi_Q - 1 =")
    sp.pprint(sp.factor(DeltaQ))

    expect_zero(
        "Delta_Q - sigma_W(1 - 9 gamma_W)/(1 - sigma_W)",
        sp.simplify(DeltaQ - sigmaW * (1 - 9 * gammaW) / (1 - sigmaW)),
    )

    subbanner("II. Exact inversion for the odd mixed-channel coefficient")
    gamma_from_chi = sp.solve(sp.Eq(chiQ, chi_hyb), gammaW)[0]
    print("gamma_W =")
    sp.pprint(sp.factor(gamma_from_chi))
    expect_zero(
        "chi_Q inversion check",
        sp.simplify(chi_hyb.subs(gammaW, gamma_from_chi) - chiQ),
    )

    subbanner("III. Natural source-map reduction")
    NQ_req = sp.simplify(1 / chi_hyb)
    print("If the natural source-map odd normalization is imposed, the required conservative factor is")
    print("N_Q^(req) =")
    sp.pprint(sp.factor(NQ_req))

    expect_zero(
        "N_Q^(req) - (1 - sigma_W)/(1 - 9 sigma_W gamma_W)",
        sp.simplify(NQ_req - (1 - sigmaW) / (1 - 9 * sigmaW * gammaW)),
    )

    subbanner("IV. Exact gauge-visible DtN observables")
    Sigma0_core = sp.simplify(3 * sigmaW)
    Sigma5_core = sp.simplify(-sigmaW * gammaW)
    Sigma5_scale = sp.simplify(sigmaW * (sp.Rational(1, 9) - gammaW))

    print("Core gauge:  Sigma_0 =")
    sp.pprint(Sigma0_core)
    print("             Sigma_5 =")
    sp.pprint(Sigma5_core)

    print("Scale gauge: Sigma_5 =")
    sp.pprint(sp.factor(Sigma5_scale))
    expect_zero(
        "scale-gauge Sigma_5 - (chi_Q - 1)(1 - sigma_W)/9",
        sp.simplify(Sigma5_scale - (chi_hyb - 1) * (1 - sigmaW) / 9),
    )

    subbanner("V. Canonical outgoing point")
    chi_can = sp.simplify(chi_hyb.subs(gammaW, sp.Rational(1, 9)))
    Delta_can = sp.simplify(DeltaQ.subs(gammaW, sp.Rational(1, 9)))
    print("gamma_W = 1/9 gives chi_Q =")
    sp.pprint(chi_can)
    print("and Delta_Q =")
    sp.pprint(Delta_can)
    expect_zero("canonical chi_Q - 1", chi_can - 1)
    expect_zero("canonical Delta_Q", Delta_can)

    print("\nInterpretation:")
    print("- The compensated selected moving-throat outgoing branch makes chi_Q a direct continuum observable.")
    print("- The only branch data needed are the static mixed loading sigma_W and the odd mixed outlet gamma_W.")
    print("- Once chi_Q is known, the natural-source-map branch fixes the required conservative normalization N_Q = 1/chi_Q.")
    print("- In the core gauge the same outgoing defect is seen directly as Sigma_0 = 3 sigma_W and Sigma_5 = - sigma_W gamma_W.")
    print("- In the scale gauge it is seen as a pure mouth renormalization S = 1-sigma_W plus the odd slot Sigma_5 = sigma_W(1/9 - gamma_W).")


if __name__ == "__main__":
    main()
