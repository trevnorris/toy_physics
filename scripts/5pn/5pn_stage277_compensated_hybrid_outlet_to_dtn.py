#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 277 — exact compensated Robin–mixed outlet dictionary into the DtN deformation gauges.

What this script does
---------------------
1. Rebuilds the explicit isotropic Robin–mixed outlet

       Lambda_hyb = Lambda_out + rho_R - sigma_W/(1 - kappa_W z^2 - i gamma_W z^5).

2. Restricts to the exact compensated canonical-even branch

       rho_R = 4 sigma_W,
       kappa_W = 1/3,

   and verifies the exact outgoing factor

       chi_Q = (1 - 9 sigma_W gamma_W)/(1 - sigma_W).

3. Shows that the same explicit outlet admits two exact DtN gauge embeddings:
      - core gauge  : S = 1,       beta = 1, Sigma_0 = 3 sigma_W, Sigma_5 = - sigma_W gamma_W
      - scale gauge : S = 1-sigma, beta = 1, Sigma_0 = 0,         Sigma_5 = sigma_W(1/9 - gamma_W)

4. Verifies that the canonical compensated branch gamma_W = 1/9 is pure-scale in
   the scale gauge, i.e. Sigma_5 = 0, which is exactly the Stage-91 robustness class.
"""


def main() -> None:
    banner("STAGE 277 — COMPENSATED HYBRID OUTLET TO DTN")

    z = sp.symbols("z", real=True)
    sigmaW, gammaW = sp.symbols("sigma_W gamma_W", positive=True, real=True)
    rhoR, kappaW = sp.symbols("rho_R kappa_W", real=True)
    S, beta = sp.symbols("S beta", positive=True, real=True)
    Sigma0, Sigma2, Sigma4, Sigma5 = sp.symbols("Sigma_0 Sigma_2 Sigma_4 Sigma_5", real=True)

    Lambda_out = -3 + z**2 / 3 + z**4 / 9 + sp.I * z**5 / 9

    subbanner("I. Exact compensated hybrid outlet")
    Lambda_hyb = sp.expand(Lambda_out + rhoR - sigmaW / (1 - kappaW * z**2 - sp.I * gammaW * z**5))
    Lambda_hyb_series = sp.series(Lambda_hyb, z, 0, 6).removeO()
    print("Lambda_hyb(z) =")
    sp.pprint(Lambda_hyb_series)

    Lambda_comp = sp.simplify(Lambda_hyb_series.subs({rhoR: 4 * sigmaW, kappaW: sp.Rational(1, 3)}))
    print("Compensated branch rho_R = 4 sigma_W, kappa_W = 1/3 gives")
    sp.pprint(Lambda_comp)

    Y_comp = sp.series(((-3 + 3 * sigmaW) / Lambda_comp), z, 0, 6).removeO()
    chi_comp = sp.simplify(27 * sp.expand(Y_comp).coeff(z, 5) / sp.I)

    print("Normalized compensated response Yhat_hyb(z) =")
    sp.pprint(Y_comp)
    print("chi_Q^(hyb) =")
    sp.pprint(sp.factor(chi_comp))
    expect_zero(
        "chi_Q^(hyb) - (1 - 9 sigma_W gamma_W)/(1 - sigma_W)",
        sp.simplify(chi_comp - (1 - 9 * sigmaW * gammaW) / (1 - sigmaW)),
    )

    subbanner("II. Core gauge embedding")
    Lambda_def_core = sp.expand(
        S * Lambda_out.subs(z, beta * z) + Sigma0 + Sigma2 * z**2 + Sigma4 * z**4 + sp.I * Sigma5 * z**5
    )
    core_subs = {
        S: 1,
        beta: 1,
        Sigma0: 3 * sigmaW,
        Sigma2: -sigmaW / 3,
        Sigma4: -sigmaW / 9,
        Sigma5: -sigmaW * gammaW,
    }
    Lambda_core = sp.expand(Lambda_def_core.subs(core_subs))

    print("Core gauge image =")
    sp.pprint(Lambda_core)
    expect_zero("core gauge matches compensated outlet", sp.simplify(Lambda_core - Lambda_comp))

    subbanner("III. Scale gauge embedding")
    scale_subs = {
        S: 1 - sigmaW,
        beta: 1,
        Sigma0: 0,
        Sigma2: 0,
        Sigma4: 0,
        Sigma5: sigmaW * (sp.Rational(1, 9) - gammaW),
    }
    Lambda_scale = sp.expand(Lambda_def_core.subs(scale_subs))
    print("Scale gauge image =")
    sp.pprint(Lambda_scale)
    expect_zero("scale gauge matches compensated outlet", sp.simplify(Lambda_scale - Lambda_comp))

    subbanner("IV. Canonical compensated outgoing branch")
    chi_can = sp.simplify(chi_comp.subs(gammaW, sp.Rational(1, 9)))
    print("On gamma_W = 1/9, chi_Q^(hyb) =")
    sp.pprint(chi_can)
    expect_zero("canonical compensated chi_Q - 1", sp.simplify(chi_can - 1))

    Sigma5_scale_can = sp.simplify(scale_subs[Sigma5].subs(gammaW, sp.Rational(1, 9)))
    print("Scale-gauge odd core slot on the canonical compensated branch =")
    sp.pprint(Sigma5_scale_can)
    expect_zero("canonical compensated scale-gauge Sigma_5", Sigma5_scale_can)

    print("\nInterpretation:")
    print("- The compensated Robin–mixed outlet is an exact explicit moving-throat outgoing branch.")
    print("- It admits a core gauge with (beta, Sigma_0, Sigma_5) = (1, 3 sigma_W, -sigma_W gamma_W).")
    print("- It also admits a scale gauge with S = 1 - sigma_W and Sigma_5 = sigma_W(1/9 - gamma_W).")
    print("- On the canonical compensated branch gamma_W = 1/9, the scale gauge becomes a pure-scale deformation,")
    print("  so the outgoing normalization stays exactly canonical by the Stage-91 robustness theorem.")


if __name__ == "__main__":
    main()
