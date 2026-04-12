#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage312_314_common import *

"""
Stage 312 — direct microscopic monomial coordinates.

This stage takes the Stage-311 three-composite zero-defect theorem and rewrites it
in the direct microscopic coherent-kernel variables
    (lambda_W, c_etaU, gamma, K_U, K_etaeff, K_Weff, mu_W, T_U).

The point is to show that the branch-adapted defect triple
    (Sigma_tr, Sigma_nt, Sigma_eta)
is already the logarithmic drift of three exact microscopic monomials, so the
actual moving-throat branch can be tested without first passing through the
intermediate observable composites.
"""

if __name__ == "__main__":
    banner("STAGE 312 — DIRECT MICROSCOPIC MONOMIAL COORDINATES")

    # Actual microscopic kernel variables.
    lambda_W, c_etaU, gamma, K_U, K_etaeff, K_Weff, mu_W, T_U = sp.symbols(
        "lambda_W c_etaU gamma K_U K_etaeff K_Weff mu_W T_U",
        positive=True,
        real=True,
    )
    # Reference branch values.
    lambda_Ws, c_etaUs, gamma_s, K_Us, K_etaeffs, K_Weffs, mu_Ws, T_Us = sp.symbols(
        "lambda_Ws c_etaUs gamma_s K_Us K_etaeffs K_Weffs mu_Ws T_Us",
        positive=True,
        real=True,
    )

    L, sigma = sp.symbols("L sigma", positive=True, real=True)
    chi0_star, deltaU_star = sp.symbols("chi0_star deltaU_star", positive=True, real=True)
    epsW_star, eps_star = sp.symbols("epsilonW_star epsilon_star", positive=True, real=True)

    mono = direct_monomials(
        lambda_W,
        c_etaU,
        gamma,
        K_U,
        K_etaeff,
        K_Weff,
        mu_W,
        T_U,
        L,
        sigma,
        chi0_star,
        deltaU_star,
        epsW_star,
        eps_star,
    )

    subbanner("I. Direct microscopic monomials on the coherent branch")
    print("C_tr,* =")
    sp.pprint(mono["Ctr"])
    print("C_nt,* =")
    sp.pprint(mono["Cnt"])
    print("epsilon_eta =")
    sp.pprint(mono["epsEta"])

    subbanner("II. Exact logarithmic microscopic drifts")
    eps = sp.symbols("eps", real=True)
    lambda_1, c_1, gamma_1, kappa_U, kappa_eta, kappa_W, mu_1, tau_1 = sp.symbols(
        "lambda_1 c_1 gamma_1 kappa_U kappa_eta kappa_W mu_1 tau_1", real=True
    )

    lin_subs = {
        lambda_W: lambda_Ws * sp.exp(eps * lambda_1),
        c_etaU: c_etaUs * sp.exp(eps * c_1),
        gamma: gamma_s * sp.exp(eps * gamma_1),
        K_U: K_Us * sp.exp(eps * kappa_U),
        K_etaeff: K_etaeffs * sp.exp(eps * kappa_eta),
        K_Weff: K_Weffs * sp.exp(eps * kappa_W),
        mu_W: mu_Ws * sp.exp(eps * mu_1),
        T_U: T_Us * sp.exp(eps * tau_1),
    }

    dlnCtr = sp.simplify(sp.diff(sp.log(mono["Ctr"].subs(lin_subs)), eps).subs({eps: 0}))
    dlnCnt = sp.simplify(sp.diff(sp.log(mono["Cnt"].subs(lin_subs)), eps).subs({eps: 0}))
    dlnepsEta = sp.simplify(sp.diff(sp.log(mono["epsEta"].subs(lin_subs)), eps).subs({eps: 0}))

    slips = microscopic_slippages(lambda_1, c_1, gamma_1, kappa_U, kappa_eta, kappa_W, mu_1, tau_1)
    branch = branch_adapted_slippages(
        chi0_star,
        deltaU_star,
        epsW_star,
        eps_star,
        slips["Sigma_chi"],
        slips["Sigma_delta"],
        slips["Sigma_eta"],
        slips["Sigma_Z"],
        slips["Sigma_eps"],
    )

    print("d ln C_tr,* =")
    sp.pprint(dlnCtr)
    print("d ln C_nt,* =")
    sp.pprint(dlnCnt)
    print("d ln epsilon_eta =")
    sp.pprint(dlnepsEta)

    expect_zero("dln C_tr,* - Sigma_tr", sp.simplify(dlnCtr - branch["Sigma_tr"]))
    expect_zero("dln C_nt,* - Sigma_nt", sp.simplify(dlnCnt - branch["Sigma_nt"]))
    expect_zero("dln epsilon_eta - Sigma_eta", sp.simplify(dlnepsEta - branch["Sigma_eta"]))

    subbanner("III. Explicit microscopic form of the branch-adapted coordinates")
    print("Sigma_tr =")
    sp.pprint(branch["Sigma_tr"])
    print("Sigma_nt =")
    sp.pprint(branch["Sigma_nt"])
    print("Sigma_eta =")
    sp.pprint(branch["Sigma_eta"])

    subbanner("IV. Exact zero-defect equivalence in monomial language")
    print("Theta_1 = Xi_1 = R_1 = 0  iff")
    print("  d ln C_tr,* = 0,")
    print("  d ln C_nt,* = 0,")
    print("  d ln epsilon_eta = 0.")
    print()
    print("So the actual moving-throat branch is now tested directly by microscopic")
    print("monomial rigidity in")
    print("  C_tr,*,   C_nt,*,   epsilon_eta,")
    print("rather than by a larger slippage ledger.")
