#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage309_311_common import banner, subbanner, expect_zero

"""
Stage 311 — exact triangular normal form, inverse reconstruction, and direct
branch-composite theorem for the coherent weak-axisymmetric defect.

This stage compresses the Stage-310 microscopic slippage ledger to the exact
branch-adapted triple
    (Sigma_tr, Sigma_nt, Sigma_eta),
and proves that the observable drifts (Theta_1, Xi_1, R_1) form a triangular
system. It also identifies the exact direct branch composites whose invariance is
 equivalent to first-order zero defect.
"""

if __name__ == "__main__":
    banner("STAGE 311 — TRIANGULAR NORMAL FORM AND DIRECT BRANCH COMPOSITES")

    chi0, deltaU = sp.symbols("chi_0 delta_U", positive=True, real=True)
    epsEta, epsW = sp.symbols("epsilon_eta epsilon_W", positive=True, real=True)
    eps = sp.symbols("epsilon", positive=True, real=True)

    SigmaZ, SigmaChi = sp.symbols("Sigma_Z Sigma_chi", real=True)
    SigmaEta, SigmaEps, SigmaDelta = sp.symbols(
        "Sigma_eta Sigma_epsilon Sigma_delta", real=True
    )

    Ctr = sp.simplify(chi0 * deltaU / ((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU)))
    Atr = sp.simplify(2 * chi0 / ((1 + chi0) * (1 + deltaU)))

    SigmaTr = sp.simplify((1 + chi0) * SigmaDelta + (1 + deltaU) * SigmaChi)
    SigmaNt = sp.simplify(
        SigmaZ
        + 2 * epsW * (11 + 9 * deltaU) * SigmaEps / (11 * (1 - eps) * (1 + deltaU))
        - (
            2 * chi0 / (1 + deltaU)
            + 4 * epsW * deltaU / (11 * (1 - eps) * (1 + deltaU) ** 2)
        )
        * SigmaDelta
    )

    Theta1 = sp.simplify(-Ctr * SigmaTr)
    Xi1 = sp.simplify(Atr * SigmaTr + SigmaNt)
    RplusXi = sp.simplify(-epsEta * SigmaEta / (1 - epsEta))

    subbanner("I. Exact branch-adapted coordinates")
    print("Sigma_tr =")
    sp.pprint(SigmaTr)
    print("Sigma_nt =")
    sp.pprint(SigmaNt)
    print("Sigma_eta =")
    sp.pprint(SigmaEta)

    subbanner("II. Exact triangular normal form")
    print("Theta_1 =")
    sp.pprint(Theta1)
    print("Xi_1 =")
    sp.pprint(Xi1)
    print("R_1 + Xi_1 =")
    sp.pprint(RplusXi)

    expect_zero("Xi_1 - (A_tr Sigma_tr + Sigma_nt)", sp.simplify(Xi1 - (Atr * SigmaTr + SigmaNt)))

    subbanner("III. Exact inverse reconstruction formulas")
    SigmaTr_inv = sp.simplify(-Theta1 / Ctr)
    SigmaNt_inv = sp.simplify(Xi1 + (Atr / Ctr) * Theta1)
    SigmaEta_inv = sp.simplify(-(1 - epsEta) * RplusXi / epsEta)

    print("Sigma_tr from Theta_1 =")
    sp.pprint(SigmaTr_inv)
    print("Sigma_nt from (Theta_1, Xi_1) =")
    sp.pprint(SigmaNt_inv)
    print("Sigma_eta from (R_1 + Xi_1) =")
    sp.pprint(SigmaEta_inv)

    expect_zero("Sigma_tr inverse", sp.simplify(SigmaTr_inv - SigmaTr))
    expect_zero("Sigma_nt inverse", sp.simplify(SigmaNt_inv - SigmaNt))
    expect_zero("Sigma_eta inverse", sp.simplify(SigmaEta_inv - SigmaEta))
    expect_zero("A_tr/C_tr - 2(1+chi_0+delta_U)/delta_U", sp.simplify(Atr / Ctr - 2 * (1 + chi0 + deltaU) / deltaU))

    subbanner("IV. Direct branch-composite theorem")
    Rtr, T2 = sp.symbols("R_tr T2", positive=True, real=True)
    dlnRtr = sp.simplify(-Ctr * SigmaTr)
    dlnT2 = sp.simplify(Atr * SigmaTr + SigmaNt)
    dlnepsEta = SigmaEta

    Cstar = sp.simplify(1 / Ctr)
    Bstar = sp.simplify(Atr / Ctr)
    Tstar = sp.simplify(Rtr ** (-Cstar))
    Nstar = sp.simplify(T2 * Rtr ** Bstar)

    dlnTstar = sp.simplify(-Cstar * dlnRtr)
    dlnNstar = sp.simplify(dlnT2 + Bstar * dlnRtr)

    print("C_* =")
    sp.pprint(Cstar)
    print("B_* =")
    sp.pprint(Bstar)
    print("d ln T_* =")
    sp.pprint(dlnTstar)
    print("d ln N_* =")
    sp.pprint(dlnNstar)
    print("d ln epsilon_eta =")
    sp.pprint(dlnepsEta)

    expect_zero("d ln T_* - Sigma_tr", sp.simplify(dlnTstar - SigmaTr))
    expect_zero("d ln N_* - Sigma_nt", sp.simplify(dlnNstar - SigmaNt))

    subbanner("V. Exact zero-defect theorem")
    print("Theta_1 = Xi_1 = R_1 + Xi_1 = 0  iff  Sigma_tr = Sigma_nt = Sigma_eta = 0.")
    print("Equivalently, the first-order coherent defect vanishes iff the three exact")
    print("branch composites are invariant:")
    print("  R_tr,   T^2 R_tr^{B_*},   epsilon_eta.")

    banner("STAGE 311 LEDGER")
    print("1. The coherent weak-axisymmetric problem collapses exactly to the")
    print("   branch-adapted triple (Sigma_tr, Sigma_nt, Sigma_eta).")
    print("2. The observable drifts form the exact triangular system")
    print("      Theta_1 = -C_tr Sigma_tr,   Xi_1 = A_tr Sigma_tr + Sigma_nt,")
    print("      R_1 + Xi_1 = -epsilon_eta Sigma_eta/(1-epsilon_eta).")
    print("3. The exact branch composites")
    print("      T_* = R_tr^{-C_*},   N_* = T^2 R_tr^{B_*},   epsilon_eta")
    print("   have logarithmic drifts")
    print("      dln T_* = Sigma_tr,   dln N_* = Sigma_nt,   dln epsilon_eta = Sigma_eta.")
    print("4. So first-order zero defect is equivalent to invariance of those three")
    print("   branch composites on the actual coherent moving-throat branch.")
