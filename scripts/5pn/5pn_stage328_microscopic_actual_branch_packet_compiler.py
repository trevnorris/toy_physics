#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage326_328_common import *

"""
Stage 328 — exact microscopic actual-branch packet compiler.

This stage closes the loop between the earlier microscopic monomial theorem and the
newer actual coherent placement-map packet compiler.

It shows three exact equivalences:

1. the finite orbit packet (q_tr,q_nt,q_eta) can be computed directly from the
   microscopic kernel state using the actual placement variables,
2. the infinitesimal orbit defect packet (Theta_1, Xi_1, R_1) computed from the
   actual placement drifts is exactly the same packet as the older Stage-313/314
   microscopic compatibility ledger,
3. the full actual branch splits cleanly into
      - an orbit-lock packet blind to zeta,
      - a support packet that depends on zeta but not on the orbit tester.
"""

if __name__ == "__main__":
    banner("STAGE 328 — MICROSCOPIC ACTUAL-BRANCH PACKET COMPILER")

    # Microscopic state symbols
    lambdaW, cetaU, gamma = sp.symbols("lambda_W c_etaU gamma", positive=True, real=True)
    KU, Ketaeff, KWeff, muW = sp.symbols("K_U K_etaeff K_Weff mu_W", positive=True, real=True)
    TU, L, sigma = sp.symbols("T_U L sigma", positive=True, real=True)
    lambdaPhi, KPhieff = sp.symbols("lambda_phi K_phi_eff", positive=True, real=True)
    G, cs, a, c = sp.symbols("G c_s a c", positive=True, real=True)
    Lambda0 = lambda0_constant(G, cs, a, c)

    lambdaW_ref, cetaU_ref, gamma_ref = sp.symbols("lambda_W_ref c_etaU_ref gamma_ref", positive=True, real=True)
    KU_ref, Ketaeff_ref, KWeff_ref, muW_ref = sp.symbols("K_U_ref K_etaeff_ref K_Weff_ref mu_W_ref", positive=True, real=True)
    TU_ref = sp.symbols("T_U_ref", positive=True, real=True)
    lambdaPhi_ref, KPhieff_ref = sp.symbols("lambda_phi_ref K_phi_eff_ref", positive=True, real=True)

    state = microscopic_coherent_placement_state(
        lambdaW, cetaU, gamma, KU, Ketaeff, KWeff, muW, TU, L, sigma, Lambda0, lambdaPhi, KPhieff
    )
    state_ref = microscopic_coherent_placement_state(
        lambdaW_ref, cetaU_ref, gamma_ref, KU_ref, Ketaeff_ref, KWeff_ref, muW_ref, TU_ref, L, sigma, Lambda0, lambdaPhi_ref, KPhieff_ref
    )

    # Microscopic drift symbols
    lambda1, c1, gamma1 = sp.symbols("lambda_1 c_1 gamma_1", real=True)
    kappaU, kappaEta, kappaW = sp.symbols("kappa_U kappa_eta kappa_W", real=True)
    mu1, tau1 = sp.symbols("mu_1 tau_1", real=True)
    phi1, kappaPhi = sp.symbols("phi_1 kappa_phi", real=True)
    drifts = microscopic_coherent_placement_drifts(
        lambda1, c1, gamma1, kappaU, kappaEta, kappaW, mu1, tau1, phi1, kappaPhi
    )

    # Reference-branch/star data
    chi0_star, deltaU_star, epsW_star, epsEta_star = sp.symbols(
        "chi0_star deltaU_star epsilon_W_star epsilon_eta_star",
        positive=True,
        real=True,
    )
    eps_star = split_epsilon(epsW_star, deltaU_star)

    q = microscopic_finite_packet(state, state_ref, chi0_star, deltaU_star, epsW_star)
    direct_mono = direct_monomials(
        lambdaW, cetaU, gamma, KU, Ketaeff, KWeff, muW, TU, L, sigma,
        chi0_star, deltaU_star, epsW_star, eps_star,
    )
    direct_mono_ref = direct_monomials(
        lambdaW_ref, cetaU_ref, gamma_ref, KU_ref, Ketaeff_ref, KWeff_ref, muW_ref, TU_ref, L, sigma,
        chi0_star, deltaU_star, epsW_star, eps_star,
    )

    packet = microscopic_packet_from_placement_drifts(
        drifts,
        chi0_star,
        deltaU_star,
        epsW_star,
        epsEta_star,
    )
    direct_sigmas = microscopic_direct_sigmas(
        lambda1, c1, gamma1, kappaU, kappaEta, kappaW, mu1, tau1,
        chi0_star, deltaU_star, epsW_star, epsEta_star,
    )
    direct_defects = defect_map_from_sigmas(
        direct_sigmas["Sigma_tr"],
        direct_sigmas["Sigma_nt"],
        direct_sigmas["Sigma_eta"],
        chi0_star,
        deltaU_star,
        epsEta_star,
    )

    obs = coherent_branch_observables(state["chi0"], state["deltaU"], state["ZW"], state["epsW"], state["epsEta"], state["Lambda"])
    supp = support_packet(state["chi0"], state["deltaU"], state["ZW"], state["epsW"], state["epsEta"], state["Lambda"], state["zeta"])

    subbanner("I. Exact finite orbit packet from the microscopic actual branch state")
    print("q_tr =")
    sp.pprint(q["qtr"])
    print("q_nt =")
    sp.pprint(q["qnt"])
    print("q_eta =")
    sp.pprint(q["qeta"])

    qtr_direct = sp.simplify(sp.log(direct_mono["Ctr"] / direct_mono_ref["Ctr"]))
    qnt_direct = sp.simplify(sp.log(direct_mono["Cnt"] / direct_mono_ref["Cnt"]))
    qeta_direct = sp.simplify(sp.log(direct_mono["epsEta"] / direct_mono_ref["epsEta"]))

    expect_zero(
        "q_tr exact microscopic equivalence",
        sp.simplify(sp.expand_log(q["qtr"], force=True) - sp.expand_log(qtr_direct, force=True)),
    )
    expect_zero(
        "q_nt exact microscopic equivalence",
        sp.simplify(sp.expand_log(q["qnt"], force=True) - sp.expand_log(qnt_direct, force=True)),
    )
    expect_zero(
        "q_eta exact microscopic equivalence",
        sp.simplify(sp.expand_log(q["qeta"], force=True) - sp.expand_log(qeta_direct, force=True)),
    )

    subbanner("II. Exact orbit defect packet from microscopic actual-branch drifts")
    print("Sigma_tr =")
    sp.pprint(packet["Sigma_tr"])
    print("Sigma_nt =")
    sp.pprint(packet["Sigma_nt"])
    print("Sigma_eta =")
    sp.pprint(packet["Sigma_eta"])
    print("Theta_1 =")
    sp.pprint(packet["Theta1"])
    print("Xi_1 =")
    sp.pprint(packet["Xi1"])
    print("R_1 =")
    sp.pprint(packet["R1"])

    expect_zero("Sigma_tr exact microscopic equivalence", packet["Sigma_tr"] - direct_sigmas["Sigma_tr"])
    expect_zero("Sigma_nt exact microscopic equivalence", packet["Sigma_nt"] - direct_sigmas["Sigma_nt"])
    expect_zero("Sigma_eta exact microscopic equivalence", packet["Sigma_eta"] - direct_sigmas["Sigma_eta"])
    expect_zero("Theta_1 exact microscopic equivalence", packet["Theta1"] - direct_defects["Theta1"])
    expect_zero("Xi_1 exact microscopic equivalence", packet["Xi1"] - direct_defects["Xi1"])
    expect_zero("R_1 exact microscopic equivalence", packet["R1"] - direct_defects["R1"])

    subbanner("III. Exact direct-observable orbit-lock form on the microscopic actual branch")
    dln_Rtr = sp.simplify(packet["Theta1"])
    dln_Rtarget = sp.simplify(packet["R1"])
    print("R_tr =")
    sp.pprint(obs["Rtr"])
    print("R_target =")
    sp.pprint(obs["Rtarget"])
    print("epsilon_eta =")
    sp.pprint(state["epsEta"])
    print("Orbit lock is exactly:")
    print("  dln R_tr = 0,")
    print("  dln R_target = 0,")
    print("  dln epsilon_eta = 0.")

    subbanner("IV. Exact support packet is separate and blind to the orbit tester")
    print("zeta =")
    sp.pprint(state["zeta"])
    print("M_mix =")
    sp.pprint(supp["Mmix"])
    print("S(zeta;epsilon) =")
    sp.pprint(supp["S"])
    print("M_tr =")
    sp.pprint(supp["Mtr"])
    expect_zero("partial_{lambda_phi} q_tr", sp.diff(q["qtr"], lambdaPhi))
    expect_zero("partial_{lambda_phi} q_nt", sp.diff(q["qnt"], lambdaPhi))
    expect_zero("partial_{lambda_phi} q_eta", sp.diff(q["qeta"], lambdaPhi))
    expect_zero("partial_{K_phi_eff} q_tr", sp.diff(q["qtr"], KPhieff))
    expect_zero("partial_{K_phi_eff} q_nt", sp.diff(q["qnt"], KPhieff))
    expect_zero("partial_{K_phi_eff} q_eta", sp.diff(q["qeta"], KPhieff))
    expect_zero("partial_{phi1} Theta_1", sp.diff(packet["Theta1"], phi1))
    expect_zero("partial_{phi1} Xi_1", sp.diff(packet["Xi1"], phi1))
    expect_zero("partial_{phi1} R_1", sp.diff(packet["R1"], phi1))
    expect_zero("partial_{kappa_phi} Theta_1", sp.diff(packet["Theta1"], kappaPhi))
    expect_zero("partial_{kappa_phi} Xi_1", sp.diff(packet["Xi1"], kappaPhi))
    expect_zero("partial_{kappa_phi} R_1", sp.diff(packet["R1"], kappaPhi))

    subbanner("V. Interpretation")
    print("The actual coherent moving-throat branch now ends at one exact microscopic")
    print("two-packet split:")
    print("  orbit packet  = (q_tr,q_nt,q_eta) or equivalently (Theta_1,Xi_1,R_1),")
    print("  support packet = (zeta; M_mix, S, M_tr).")
    print()
    print("The orbit packet compiled from the actual placement drifts is exactly the same")
    print("packet as the older Stage-313/314 microscopic compatibility ledger.  So there")
    print("is no hidden extra weak-axisymmetric obstruction between the microscopic kernel")
    print("data and the actual coherent branch variables.  What remains is to compute the")
    print("actual finite state and drift data from the completed moving-throat operator.")
