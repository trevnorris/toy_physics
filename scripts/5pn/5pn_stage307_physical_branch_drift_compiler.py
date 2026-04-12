#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage300_302_common import banner, subbanner, expect_zero

"""
Stage 307 — exact first-order drift compiler from the physical branch observables.

This stage linearizes the Stage-306 direct coherent-branch observables and proves
that their logarithmic drifts reproduce the earlier physical defect laws exactly.

Main outputs:
  d ln R_tr,
  d ln R_target,
  d ln epsilon_eta,
all directly in the physical coherent variables
  (chi_0, delta_U, epsilon_eta, Z_W, Omega_W^2, epsilon).
"""

if __name__ == "__main__":
    banner("STAGE 307 — PHYSICAL BRANCH DRIFT COMPILER")

    chi0, deltaU = sp.symbols("chi_0 delta_U", positive=True, real=True)
    epsEta, eps = sp.symbols("epsilon_eta epsilon", positive=True, real=True)
    ZW, OmW2, Lambda0 = sp.symbols("Z_W Omega_W2 Lambda_0", positive=True, real=True)

    dlnchi0, dlndeltaU = sp.symbols("dlnchi_0 dlndelta_U", real=True)
    dlnepsEta, dlneps = sp.symbols("dlnepsilon_eta dlnepsilon", real=True)
    dlnZW, dlnOmW2 = sp.symbols("dlnZ_W dlnOmega_W2", real=True)
    epslin = sp.symbols("epslin", real=True)

    Rtr = sp.simplify((1 + chi0 / (1 + deltaU)) / (1 + chi0))
    Rtarget = sp.simplify(Lambda0 * OmW2 * (1 - epsEta) * (1 - eps) ** 2 / (ZW * (1 + chi0) ** 2))

    lin_subs = {
        chi0: chi0 * sp.exp(epslin * dlnchi0),
        deltaU: deltaU * sp.exp(epslin * dlndeltaU),
        epsEta: epsEta * sp.exp(epslin * dlnepsEta),
        eps: eps * sp.exp(epslin * dlneps),
        ZW: ZW * sp.exp(epslin * dlnZW),
        OmW2: OmW2 * sp.exp(epslin * dlnOmW2),
    }

    dlnRtr = sp.simplify(sp.diff(sp.log(Rtr.subs(lin_subs)), epslin).subs({epslin: 0}))
    dlnRtarget = sp.simplify(sp.diff(sp.log(Rtarget.subs(lin_subs)), epslin).subs({epslin: 0}))
    dlnepsEta_exact = dlnepsEta

    subbanner("I. Exact logarithmic drifts of the physical branch observables")
    print("d ln R_tr =")
    sp.pprint(dlnRtr)
    print("d ln R_target =")
    sp.pprint(dlnRtarget)
    print("d ln epsilon_eta =")
    sp.pprint(dlnepsEta_exact)

    Ctr = sp.simplify(chi0 * deltaU / ((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU)))
    dlnRtr_expected = sp.simplify(-Ctr * ((1 + deltaU) * dlnchi0 + (1 + chi0) * dlndeltaU))
    dlnRtarget_expected = sp.simplify(
        dlnOmW2
        - dlnZW
        - 2 * chi0 * dlnchi0 / (1 + chi0)
        - 2 * eps * dlneps / (1 - eps)
        - epsEta * dlnepsEta / (1 - epsEta)
    )

    expect_zero("d ln R_tr closed form", sp.simplify(dlnRtr - dlnRtr_expected))
    expect_zero("d ln R_target closed form", sp.simplify(dlnRtarget - dlnRtarget_expected))

    subbanner("II. Exact physical defect triplet in direct drift language")
    Theta1 = sp.simplify(dlnRtr)
    Xi1 = sp.simplify(-dlnRtarget - epsEta * dlnepsEta / (1 - epsEta))
    R1 = sp.simplify(dlnRtarget)

    print("Theta_1 =")
    sp.pprint(Theta1)
    print("Xi_1 =")
    sp.pprint(Xi1)
    print("R_1 =")
    sp.pprint(R1)

    expect_zero(
        "Xi_1 - [dlnZ_W - dlnOmega_W^2 + 2 chi_0 dlnchi_0/(1+chi_0) + 2 epsilon dlnepsilon/(1-epsilon)]",
        sp.simplify(Xi1 - (dlnZW - dlnOmW2 + 2 * chi0 * dlnchi0 / (1 + chi0) + 2 * eps * dlneps / (1 - eps))),
    )
    expect_zero("R_1 + Xi_1 + epsilon_eta dlnepsilon_eta/(1-epsilon_eta)", sp.simplify(R1 + Xi1 + epsEta * dlnepsEta / (1 - epsEta)))

    subbanner("III. Exact match to the earlier physical-branch slope formulas")
    zetaZ, omegaW = sp.symbols("zeta_Z omega_W", real=True)
    chi1, deltaU1 = sp.symbols("chi_1 delta_U1", real=True)
    eps1, eta1 = sp.symbols("epsilon_1 eta_1", real=True)

    Theta1_phys = sp.simplify(Theta1.subs({dlnchi0: chi1 / chi0, dlndeltaU: deltaU1 / deltaU}))
    Xi1_phys = sp.simplify(
        Xi1.subs({
            dlnZW: zetaZ,
            dlnOmW2: omegaW,
            dlnchi0: chi1 / chi0,
            dlneps: eps1 / eps,
            dlnepsEta: eta1 / epsEta,
        })
    )
    R1_phys = sp.simplify(
        R1.subs({
            dlnZW: zetaZ,
            dlnOmW2: omegaW,
            dlnchi0: chi1 / chi0,
            dlneps: eps1 / eps,
            dlnepsEta: eta1 / epsEta,
        })
    )

    Theta1_expected = sp.simplify(
        -(chi0 * (1 + chi0) * deltaU1 + deltaU * (1 + deltaU) * chi1)
        / ((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU))
    )
    Xi1_expected = sp.simplify(zetaZ - omegaW + 2 * chi1 / (1 + chi0) + 2 * eps1 / (1 - eps))
    R1_expected = sp.simplify(omegaW - eta1 / (1 - epsEta) - zetaZ - 2 * chi1 / (1 + chi0) - 2 * eps1 / (1 - eps))

    print("Theta_1 in the earlier slope notation =")
    sp.pprint(Theta1_phys)
    print("Xi_1 in the earlier slope notation =")
    sp.pprint(Xi1_phys)
    print("R_1 in the earlier slope notation =")
    sp.pprint(R1_phys)

    expect_zero("Theta_1 slope match", sp.simplify(Theta1_phys - Theta1_expected))
    expect_zero("Xi_1 slope match", sp.simplify(Xi1_phys - Xi1_expected))
    expect_zero("R_1 slope match", sp.simplify(R1_phys - R1_expected))

    banner("STAGE 307 LEDGER")
    print("1. The physical branch observables have exact drift laws")
    print("      d ln R_tr     = -C_tr[(1+delta_U)dlnchi_0 + (1+chi_0)dlndelta_U],")
    print("      d ln R_target = dlnOmega_W^2 - dlnZ_W - 2chi_0 dlnchi_0/(1+chi_0)")
    print("                       - 2epsilon dlnepsilon/(1-epsilon)")
    print("                       - epsilon_eta dlnepsilon_eta/(1-epsilon_eta),")
    print("      d ln epsilon_eta = dlnepsilon_eta.")
    print("2. Composed with the Stage-304 defect map, these reproduce the earlier")
    print("   physical-branch formulas for (Theta_1, Xi_1, R_1) exactly.")
    print("3. So the first-order orbit-lock test is now written directly in the physical")
    print("   coherent variables rather than only in the abstract quotient packet.")
