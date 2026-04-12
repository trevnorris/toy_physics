#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage309_311_common import banner, subbanner, expect_zero

"""
Stage 310 — microscopic coherent-kernel slippage decomposition.

This stage pushes the Stage-309 physical defect law down to the actual coherent
kernel couplings. It identifies the four mixed/outgoing slippages carrying the
weak-axisymmetric grouped defect, plus the selected-branch dressing slippage.

Main outputs:
  Sigma_Z, Sigma_chi, Sigma_epsilon, Sigma_delta, Sigma_eta,
  Xi_1 in microscopic slippage form,
  and the exact tracking combination Sigma_tr.
"""

if __name__ == "__main__":
    banner("STAGE 310 — MICROSCOPIC SLIPPAGE DECOMPOSITION")

    chi0, deltaU = sp.symbols("chi_0 delta_U", positive=True, real=True)
    epsEta, epsW = sp.symbols("epsilon_eta epsilon_W", positive=True, real=True)
    split = sp.simplify(1 - sp.Rational(2, 11) * deltaU / (1 + deltaU))
    eps = sp.simplify(epsW * split)

    lam1, c1, gam1 = sp.symbols("lambda_1 c_1 gamma_1", real=True)
    kapU, kapEta, kapW = sp.symbols("kappa_U kappa_eta kappa_W", real=True)
    mu1, tau1 = sp.symbols("mu_1 tau_1", real=True)

    subbanner("I. Exact physical-branch drifts from the microscopic variables")
    zetaZ = sp.simplify(2 * lam1 - kapEta - kapW)
    omegaW = sp.simplify(kapW - mu1)
    chi1 = sp.simplify(chi0 * (gam1 + c1 - kapU))
    eta1 = sp.simplify(epsEta * (2 * c1 - kapU - kapEta))
    epsW1 = sp.simplify(epsW * (2 * gam1 + 2 * lam1 - kapU - kapW))
    deltaU1 = sp.simplify(deltaU * (tau1 - kapU))

    print("zeta_Z =")
    sp.pprint(zetaZ)
    print("omega_W =")
    sp.pprint(omegaW)
    print("chi_1 =")
    sp.pprint(chi1)
    print("eta_1 =")
    sp.pprint(eta1)
    print("epsilonW_1 =")
    sp.pprint(epsW1)
    print("deltaU_1 =")
    sp.pprint(deltaU1)

    subbanner("II. Exact split-blocking drift and physical defect law")
    tau = sp.symbols("tau", real=True)
    eps_t = sp.simplify(
        (epsW + tau * epsW1)
        * (1 - sp.Rational(2, 11) * (deltaU + tau * deltaU1) / (1 + deltaU + tau * deltaU1))
    )
    eps1 = sp.simplify(sp.diff(eps_t, tau).subs({tau: 0}))
    eps1_expected = sp.simplify(
        split * epsW1 - sp.Rational(2, 11) * epsW * deltaU1 / (1 + deltaU) ** 2
    )
    expect_zero("epsilon_1 transport law", sp.simplify(eps1 - eps1_expected))
    print("epsilon_1 =")
    sp.pprint(eps1)

    Xi1_phys = sp.simplify(zetaZ - omegaW + 2 * chi1 / (1 + chi0) + 2 * eps1 / (1 - eps))
    print("Xi_1 =")
    sp.pprint(Xi1_phys)

    subbanner("III. Exact microscopic slippage variables")
    SigmaZ = sp.simplify(2 * lam1 + mu1 - kapEta - 2 * kapW)
    SigmaChi = sp.simplify(gam1 + c1 - kapU)
    SigmaEta = sp.simplify(2 * c1 - kapU - kapEta)
    SigmaEps = sp.simplify(2 * gam1 + 2 * lam1 - kapU - kapW)
    SigmaDelta = sp.simplify(tau1 - kapU)

    print("Sigma_Z =")
    sp.pprint(SigmaZ)
    print("Sigma_chi =")
    sp.pprint(SigmaChi)
    print("Sigma_eta =")
    sp.pprint(SigmaEta)
    print("Sigma_epsilon =")
    sp.pprint(SigmaEps)
    print("Sigma_delta =")
    sp.pprint(SigmaDelta)

    expect_zero("Sigma_Z - (zeta_Z - omega_W)", sp.simplify(SigmaZ - (zetaZ - omegaW)))
    expect_zero("Sigma_chi - chi_1/chi_0", sp.simplify(SigmaChi - chi1 / chi0))
    expect_zero("Sigma_eta - eta_1/epsilon_eta", sp.simplify(SigmaEta - eta1 / epsEta))
    expect_zero("Sigma_epsilon - epsilonW_1/epsilon_W", sp.simplify(SigmaEps - epsW1 / epsW))
    expect_zero("Sigma_delta - deltaU_1/delta_U", sp.simplify(SigmaDelta - deltaU1 / deltaU))

    Xi1_slip = sp.simplify(
        SigmaZ
        + 2 * chi0 * SigmaChi / (1 + chi0)
        + 2 * epsW / (1 - eps)
        * (
            (11 + 9 * deltaU) * SigmaEps / (11 * (1 + deltaU))
            - 2 * deltaU * SigmaDelta / (11 * (1 + deltaU) ** 2)
        )
    )
    print("Xi_1 in slippage form =")
    sp.pprint(Xi1_slip)
    expect_zero("Xi_1 microscopic slippage decomposition", sp.simplify(Xi1_phys - Xi1_slip))

    subbanner("IV. Exact tracking slippage combination")
    SigmaTr = sp.simplify((1 + chi0) * SigmaDelta + (1 + deltaU) * SigmaChi)
    Theta1 = sp.simplify(
        -(chi0 * deltaU) * SigmaTr / ((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU))
    )
    print("Sigma_tr =")
    sp.pprint(SigmaTr)
    print("Theta_1 =")
    sp.pprint(Theta1)

    Theta1_phys = sp.simplify(
        -(chi0 * (1 + chi0) * deltaU1 + deltaU * (1 + deltaU) * chi1)
        / ((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU))
    )
    expect_zero("Theta_1 tracking split", sp.simplify(Theta1 - Theta1_phys))

    banner("STAGE 310 LEDGER")
    print("1. The coherent weak-axisymmetric grouped defect depends only on four")
    print("   microscopic mixed/outgoing slippages")
    print("      (Sigma_Z, Sigma_chi, Sigma_epsilon, Sigma_delta),")
    print("   plus the selected-branch dressing slippage Sigma_eta.")
    print("2. In exact microscopic form,")
    print("      Xi_1 = Sigma_Z + 2 chi_0 Sigma_chi/(1+chi_0)")
    print("             + 2 epsilon_W/(1-epsilon) [ (11+9 delta_U) Sigma_epsilon/(11(1+delta_U))")
    print("                                        - 2 delta_U Sigma_delta/(11(1+delta_U)^2) ].")
    print("3. The tracking-factor drift depends only on")
    print("      Sigma_tr = (1+chi_0) Sigma_delta + (1+delta_U) Sigma_chi.")
    print("4. So exact tracking rigidity is necessary but not sufficient: the grouped")
    print("   defect also depends on the nontracking mixed/outgoing slippages.")
