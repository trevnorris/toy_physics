#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage326_328_common import *
from fivepn_stage319_322_common import microscopic_state_to_reduced, coherent_support_ratio

"""
Stage 326 — exact microscopic coherent placement-state extractor.

This stage closes the remaining gap between the older microscopic coherent-kernel
state and the newer actual coherent placement variables used by the Stage 323–325
orbit/support packet compilers.

The actual coherent branch variables are
    (chi0, deltaU, Z_W, epsilon_W, epsilon_eta, Lambda, zeta),
not just the reduced Stage-318 state
    (chi0, deltaU, ZhatW, epsilon_W, epsilon_eta).

The new point is that the actual placement state is itself an exact microscopic
kernel quotient, and the older reduced variable is recovered by the exact identity
    ZhatW = Z_W * Lambda0 / Lambda.
"""

if __name__ == "__main__":
    banner("STAGE 326 — MICROSCOPIC COHERENT PLACEMENT-STATE EXTRACTOR")

    lambdaW, cetaU, gamma = sp.symbols("lambda_W c_etaU gamma", positive=True, real=True)
    KU, Ketaeff, KWeff, muW = sp.symbols("K_U K_etaeff K_Weff mu_W", positive=True, real=True)
    TU, L, sigma = sp.symbols("T_U L sigma", positive=True, real=True)
    lambdaPhi, KPhieff = sp.symbols("lambda_phi K_phi_eff", positive=True, real=True)
    G, cs, a, c = sp.symbols("G c_s a c", positive=True, real=True)

    Lambda0 = lambda0_constant(G, cs, a, c)
    state = microscopic_coherent_placement_state(
        lambdaW, cetaU, gamma, KU, Ketaeff, KWeff, muW, TU, L, sigma, Lambda0, lambdaPhi, KPhieff
    )

    red_old = microscopic_state_to_reduced(lambdaW, cetaU, gamma, KU, Ketaeff, KWeff, muW, TU, L, sigma)
    zeta_old = coherent_support_ratio(lambdaPhi, KPhieff, lambdaW, KWeff)
    obs = coherent_branch_observables(state["chi0"], state["deltaU"], state["ZW"], state["epsW"], state["epsEta"], state["Lambda"])
    supp = support_packet(state["chi0"], state["deltaU"], state["ZW"], state["epsW"], state["epsEta"], state["Lambda"], state["zeta"])

    subbanner("I. Exact microscopic -> actual coherent placement state")
    print("chi0 =")
    sp.pprint(state["chi0"])
    print("deltaU =")
    sp.pprint(state["deltaU"])
    print("Z_W =")
    sp.pprint(state["ZW"])
    print("epsilon_W =")
    sp.pprint(state["epsW"])
    print("epsilon_eta =")
    sp.pprint(state["epsEta"])
    print("Lambda =")
    sp.pprint(state["Lambda"])
    print("zeta =")
    sp.pprint(state["zeta"])

    subbanner("II. Exact bridge back to the older reduced state")
    print("Lambda0 =")
    sp.pprint(Lambda0)
    print("ZhatW from the actual placement state =")
    sp.pprint(state["ZhatW"])
    print("ZhatW from the older microscopic reduction =")
    sp.pprint(red_old["ZhatW"])
    expect_zero("ZhatW equivalence", state["ZhatW"] - red_old["ZhatW"])
    expect_zero("chi0 equivalence", state["chi0"] - red_old["chi0"])
    expect_zero("deltaU equivalence", state["deltaU"] - red_old["deltaU"])
    expect_zero("epsilon_W equivalence", state["epsW"] - red_old["epsW"])
    expect_zero("epsilon_eta equivalence", state["epsEta"] - red_old["epsEta"])
    expect_zero("zeta equivalence", state["zeta"] - zeta_old)

    subbanner("III. Exact physical branch observables in microscopic kernel variables")
    print("epsilon =")
    sp.pprint(obs["eps"])
    print("R_tr =")
    sp.pprint(obs["Rtr"])
    print("R_target =")
    sp.pprint(obs["Rtarget"])

    subbanner("IV. Exact support packet in microscopic kernel variables")
    print("M_mix =")
    sp.pprint(supp["Mmix"])
    print("S(zeta;epsilon) =")
    sp.pprint(supp["S"])
    print("M_tr =")
    sp.pprint(supp["Mtr"])
    print("R_target * M_mix =")
    sp.pprint(supp["product"])
    expect_zero(
        "mixed-only product law",
        sp.simplify(supp["product"] - 8 * state["Lambda"] * (1 - supp["eps"]) / sp.pi**2),
    )

    subbanner("V. Interpretation")
    print("The actual coherent local branch is now completely writable as a microscopic")
    print("kernel quotient state. The orbit-lock packet depends only on")
    print("  (chi0, deltaU, Z_W, epsilon_W, epsilon_eta, Lambda),")
    print("while the support packet adds only the single extra support ratio zeta.")
    print("So the older reduced variable ZhatW is no longer fundamental; it is exactly")
    print("the quotient Z_W * Lambda0 / Lambda of two actual branch observables.")
