#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage300_302_common import banner, subbanner, expect_zero

"""
Stage 304 — exact first-order defect compiler in direct branch-observable language.

What this script does
---------------------
1. Linearizes the Stage-303 finite quotient coordinates in the logarithmic drifts of
      R_tr, R_target, epsilon_eta.
2. Composes that linearization with the exact Stage-302 quotient-to-defect map.
3. Shows that the physical first-order defect triplet is triangular in those direct
   branch-observable drifts.

Interpretation
--------------
On the coherent local tracking branch, the actual physical first-order defect test
is even smaller than the Stage-303 finite quotient packet suggests.  At linear order,
it is just a triangular map of
    (delta ln R_tr, delta ln R_target, delta ln epsilon_eta).
"""

if __name__ == "__main__":
    banner("STAGE 304 — DIRECT BRANCH DEFECT MAP")

    chi0s, deltaUs, epsEtas = sp.symbols(
        "chi0_star deltaU_star epsilon_eta_star", positive=True, real=True
    )
    dlnRtr, dlnRtarget, dlnepsEta = sp.symbols(
        "dlnR_tr dlnR_target dlnepsilon_eta", real=True
    )

    Cstar = sp.simplify(
        (1 + chi0s) * (1 + deltaUs) * (1 + chi0s + deltaUs) / (chi0s * deltaUs)
    )
    Ctr = sp.simplify(1 / Cstar)
    Bstar = sp.simplify(2 * (1 + chi0s + deltaUs) / deltaUs)
    Atr = sp.simplify(Bstar * Ctr)

    qtr = sp.simplify(-Cstar * dlnRtr)
    qnt = sp.simplify(Bstar * dlnRtr - epsEtas * dlnepsEta / (1 - epsEtas) - dlnRtarget)
    qeta = dlnepsEta

    subbanner("I. Exact quotient residuals in direct branch-observable drifts")
    print("q_tr =")
    sp.pprint(qtr)
    print("q_nt =")
    sp.pprint(qnt)
    print("q_eta =")
    sp.pprint(qeta)

    Theta1 = sp.simplify(-Ctr * qtr)
    Xi1 = sp.simplify(Atr * qtr + qnt)
    R1 = sp.simplify(-epsEtas * qeta / (1 - epsEtas) - Xi1)

    subbanner("II. Exact physical defect triplet")
    print("Theta_1 =")
    sp.pprint(Theta1)
    print("Xi_1 =")
    sp.pprint(Xi1)
    print("R_1 =")
    sp.pprint(R1)

    expect_zero("Theta_1 - dlnR_tr", sp.simplify(Theta1 - dlnRtr))
    expect_zero(
        "Xi_1 + dlnR_target + epsilon_eta_* dlnepsilon_eta/(1-epsilon_eta_*)",
        sp.simplify(Xi1 + dlnRtarget + epsEtas * dlnepsEta / (1 - epsEtas)),
    )
    expect_zero("R_1 - dlnR_target", sp.simplify(R1 - dlnRtarget))

    subbanner("III. Exact inverse map")
    dlnRtr_inv = sp.simplify(Theta1)
    dlnRtarget_inv = sp.simplify(R1)
    dlnepsEta_inv = sp.simplify(-(1 - epsEtas) * (R1 + Xi1) / epsEtas)

    print("dln R_tr =")
    sp.pprint(dlnRtr_inv)
    print("dln R_target =")
    sp.pprint(dlnRtarget_inv)
    print("dln epsilon_eta =")
    sp.pprint(dlnepsEta_inv)

    expect_zero("dlnR_tr inverse", sp.simplify(dlnRtr_inv - dlnRtr))
    expect_zero("dlnR_target inverse", sp.simplify(dlnRtarget_inv - dlnRtarget))
    expect_zero("dlnepsilon_eta inverse", sp.simplify(dlnepsEta_inv - dlnepsEta))

    banner("STAGE 304 LEDGER")
    print("1. The direct branch-observable drifts map to the physical defects by")
    print("      Theta_1 = delta ln R_tr,")
    print("      Xi_1    = - delta ln R_target - [epsilon_eta_*/(1-epsilon_eta_*)] delta ln epsilon_eta,")
    print("      R_1     = delta ln R_target.")
    print("2. The inverse map is exact and triangular, so the actual first-order defect")
    print("   test can be carried entirely in the three direct observable drifts")
    print("      (delta ln R_tr, delta ln R_target, delta ln epsilon_eta).")
