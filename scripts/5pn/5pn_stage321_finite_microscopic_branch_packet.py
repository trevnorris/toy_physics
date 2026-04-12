#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage319_322_common import *

"""
Stage 321 — exact finite microscopic branch packet.

This stage evaluates the finite Stage-318 quotient packet directly on the
microscopic coherent-kernel state
    (lambda_W, c_etaU, gamma, K_U, K_etaeff, K_Weff, mu_W, T_U).

It proves that the finite packet (q_tr,q_nt,q_eta) is exactly the finite log-ratio
packet of the three direct microscopic monomials
    (C_tr, C_nt, epsilon_eta),
and writes the physical branch observables (R_tr, R_target, epsilon_eta)
directly in microscopic variables.
"""

if __name__ == "__main__":
    banner("STAGE 321 — FINITE MICROSCOPIC BRANCH PACKET")

    lambdaW, cetaU, gamma, KU, Ketaeff, KWeff, muW, TU = sp.symbols(
        "lambda_W c_etaU gamma K_U K_etaeff K_Weff mu_W T_U", positive=True, real=True
    )
    lambdaW_ref, cetaU_ref, gamma_ref, KU_ref, Ketaeff_ref, KWeff_ref, muW_ref, TU_ref = sp.symbols(
        "lambda_W_ref c_etaU_ref gamma_ref K_U_ref K_etaeff_ref K_Weff_ref mu_W_ref T_U_ref",
        positive=True,
        real=True,
    )
    L, sigma, Lambda0 = sp.symbols("L sigma Lambda_0", positive=True, real=True)
    chi0_star, deltaU_star, epsW_star, eps_star = sp.symbols(
        "chi0_star deltaU_star epsilonW_star epsilon_star", positive=True, real=True
    )

    state = microscopic_state_to_reduced(lambdaW, cetaU, gamma, KU, Ketaeff, KWeff, muW, TU, L, sigma)
    state_ref = microscopic_state_to_reduced(lambdaW_ref, cetaU_ref, gamma_ref, KU_ref, Ketaeff_ref, KWeff_ref, muW_ref, TU_ref, L, sigma)

    q = continuum_monomial_quotients(
        state["chi0"], state["deltaU"], state["ZhatW"], state["epsW"], state["epsEta"],
        state_ref["chi0"], state_ref["deltaU"], state_ref["ZhatW"], state_ref["epsW"], state_ref["epsEta"],
        chi0_star, deltaU_star, epsW_star, eps_star,
    )

    mono = direct_microscopic_monomials(
        lambdaW, cetaU, gamma, KU, Ketaeff, KWeff, muW, TU, L, sigma,
        chi0_star, deltaU_star, epsW_star, eps_star,
    )
    mono_ref = direct_microscopic_monomials(
        lambdaW_ref, cetaU_ref, gamma_ref, KU_ref, Ketaeff_ref, KWeff_ref, muW_ref, TU_ref, L, sigma,
        chi0_star, deltaU_star, epsW_star, eps_star,
    )

    obs = physical_branch_observables_from_microscopic(
        lambdaW, cetaU, gamma, KU, Ketaeff, KWeff, muW, TU, L, sigma, Lambda0,
    )

    subbanner("I. Direct microscopic monomials")
    print("C_tr,* =")
    sp.pprint(mono["Ctr"])
    print("C_nt,* =")
    sp.pprint(mono["Cnt"])
    print("epsilon_eta =")
    sp.pprint(mono["epsEta"])

    subbanner("II. Finite quotient packet from microscopic monomials")
    qtr_mono = sp.simplify(sp.log(mono["Ctr"] / mono_ref["Ctr"]))
    qnt_mono = sp.simplify(sp.log(mono["Cnt"] / mono_ref["Cnt"]))
    qeta_mono = sp.simplify(sp.log(mono["epsEta"] / mono_ref["epsEta"]))

    print("q_tr =")
    sp.pprint(q["qtr"])
    print("q_nt =")
    sp.pprint(q["qnt"])
    print("q_eta =")
    sp.pprint(q["qeta"])

    expect_zero("q_tr monomial quotient equivalence", sp.simplify(sp.expand_log(q["qtr"], force=True) - sp.expand_log(qtr_mono, force=True)))
    expect_zero("q_nt monomial quotient equivalence", sp.simplify(sp.expand_log(q["qnt"], force=True) - sp.expand_log(qnt_mono, force=True)))
    expect_zero("q_eta monomial quotient equivalence", sp.simplify(sp.expand_log(q["qeta"], force=True) - sp.expand_log(qeta_mono, force=True)))

    subbanner("III. Direct branch observables in microscopic variables")
    print("chi0 =")
    sp.pprint(obs["chi0"])
    print("deltaU =")
    sp.pprint(obs["deltaU"])
    print("ZhatW =")
    sp.pprint(obs["ZhatW"])
    print("epsilonW =")
    sp.pprint(obs["epsW"])
    print("epsilon =")
    sp.pprint(obs["eps"])
    print("R_tr =")
    sp.pprint(obs["Rtr"])
    print("R_target =")
    sp.pprint(obs["Rtarget"])
    print("epsilon_eta =")
    sp.pprint(obs["epsEta"])

    expect_zero(
        "R_target * ZhatW * (1+chi0)^2 - Lambda_0 (1-epsilon_eta)(1-epsilon)^2",
        sp.simplify(obs["Rtarget"] * obs["ZhatW"] * (1 + obs["chi0"]) ** 2 - Lambda0 * (1 - obs["epsEta"]) * (1 - obs["eps"]) ** 2),
    )

    subbanner("IV. Interpretation")
    print("The finite Stage-318 packet is now evaluated directly on the microscopic")
    print("coherent-kernel state.  So the actual moving-throat branch does not need a")
    print("separate reduced-state postulate before the finite quotient packet can be")
    print("compiled: it is already a direct log-ratio packet of three microscopic")
    print("monomials and three direct microscopic branch observables.")
