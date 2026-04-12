#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage319_322_common import *

"""
Stage 320 — exact microscopic actual-branch tester.

This stage composes the Stage-319 microscopic->reduced drift extractor with the
Stage-317 reduced orbit tester.  The result is an exact map from the microscopic
kernel drifts directly to
    (Sigma_tr, Sigma_nt, Sigma_eta)
and therefore to
    (Theta_1, Xi_1, R_1).

It also solves the exact zero-defect compatibility ledger in microscopic drift
variables, showing that the microscopic zero-defect kernel is 5-dimensional.
"""

if __name__ == "__main__":
    banner("STAGE 320 — MICROSCOPIC ACTUAL-BRANCH TESTER")

    lambda1, c1, gamma1, kappaU, kappaEta, kappaW, mu1, tau1 = sp.symbols(
        "lambda_1 c_1 gamma_1 kappa_U kappa_eta kappa_W mu_1 tau_1", real=True
    )
    chi0_star, deltaU_star, epsW_star, eps_star, epsEta_star = sp.symbols(
        "chi0_star deltaU_star epsilonW_star epsilon_star epsilon_eta_star",
        positive=True,
        real=True,
    )

    red = microscopic_reduced_drifts(lambda1, c1, gamma1, kappaU, kappaEta, kappaW, mu1, tau1)
    sigmas = continuum_monomial_drifts(
        red["dln_chi0"],
        red["dln_deltaU"],
        red["dln_ZhatW"],
        red["dln_epsW"],
        red["dln_epsEta"],
        chi0_star,
        deltaU_star,
        epsW_star,
        eps_star,
    )
    defects = defect_map_from_sigmas(
        sigmas["Sigma_tr"], sigmas["Sigma_nt"], sigmas["Sigma_eta"],
        chi0_star, deltaU_star, epsEta_star,
    )

    mic = sp.Matrix([lambda1, c1, gamma1, kappaU, kappaEta, kappaW, mu1, tau1])
    Mmicro = sp.Matrix([
        [sp.diff(sigmas["Sigma_tr"], x) for x in mic],
        [sp.diff(sigmas["Sigma_nt"], x) for x in mic],
        [sp.diff(sigmas["Sigma_eta"], x) for x in mic],
    ])

    subbanner("I. Exact microscopic monomial-drift tester")
    print("Sigma_tr =")
    sp.pprint(sigmas["Sigma_tr"])
    print("Sigma_nt =")
    sp.pprint(sigmas["Sigma_nt"])
    print("Sigma_eta =")
    sp.pprint(sigmas["Sigma_eta"])
    print("M_micro =")
    sp.pprint(Mmicro)
    print("rank(M_micro) =", Mmicro.rank())

    subbanner("II. Exact microscopic defect packet")
    print("Theta_1 =")
    sp.pprint(defects["Theta1"])
    print("Xi_1 =")
    sp.pprint(defects["Xi1"])
    print("R_1 =")
    sp.pprint(defects["R1"])

    subbanner("III. Exact microscopic zero-defect compatibility ledger")
    sol_tau1 = sp.simplify(
        kappaU - (1 + deltaU_star) * (gamma1 + c1 - kappaU) / (1 + chi0_star)
    )
    sol_kappaEta = sp.simplify(2 * c1 - kappaU)
    E_star, F_star = monomial_exponents(chi0_star, deltaU_star, epsW_star, eps_star)
    sol_mu1 = sp.simplify(
        2 * c1 - kappaU + 2 * kappaW - 2 * lambda1
        - E_star * (2 * gamma1 + 2 * lambda1 - kappaU - kappaW)
        - F_star * (1 + deltaU_star) * (gamma1 + c1 - kappaU) / (1 + chi0_star)
    )

    print("tau_1 =")
    sp.pprint(sol_tau1)
    print("kappa_eta =")
    sp.pprint(sol_kappaEta)
    print("mu_1 =")
    sp.pprint(sol_mu1)

    sigmas_zero = continuum_monomial_drifts(
        red["dln_chi0"],
        red["dln_deltaU"].subs({tau1: sol_tau1}),
        red["dln_ZhatW"].subs({kappaEta: sol_kappaEta, mu1: sol_mu1}),
        red["dln_epsW"],
        red["dln_epsEta"].subs({kappaEta: sol_kappaEta}),
        chi0_star,
        deltaU_star,
        epsW_star,
        eps_star,
    )
    defects_zero = defect_map_from_sigmas(
        sigmas_zero["Sigma_tr"], sigmas_zero["Sigma_nt"], sigmas_zero["Sigma_eta"],
        chi0_star, deltaU_star, epsEta_star,
    )

    expect_zero("Sigma_tr on microscopic zero-defect kernel", sigmas_zero["Sigma_tr"])
    expect_zero("Sigma_nt on microscopic zero-defect kernel", sigmas_zero["Sigma_nt"])
    expect_zero("Sigma_eta on microscopic zero-defect kernel", sigmas_zero["Sigma_eta"])
    expect_zero("Theta_1 on microscopic zero-defect kernel", defects_zero["Theta1"])
    expect_zero("Xi_1 on microscopic zero-defect kernel", defects_zero["Xi1"])
    expect_zero("R_1 on microscopic zero-defect kernel", defects_zero["R1"])

    subbanner("IV. Codimension statement")
    print("The microscopic actual-branch tester has rank 3 inside the eight-drift kernel")
    print("space, so its exact zero-defect kernel is 5-dimensional.  This is the same")
    print("codimension-3 similarity-orbit closure seen earlier, now written directly in")
    print("the Stage-318 reduced actual-branch language.")
