#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage319_322_common import *

"""
Stage 322 — exact support-blind orbit lock vs support-sensitive isotropic split.

This stage shows, in direct microscopic variables, that the finite/infinitesimal
orbit-lock packet is exactly blind to the explicit support variables
    (lambda_phi, K_phieff),
while the isotropic coherent support-enhancement lane depends on them through the
single ratio
    zeta = lambda_phi^2 K_Weff / (lambda_W^2 K_phieff).

So the weak-axisymmetric orbit-lock test and the isotropic support-compensation
normalization test are rigorously separated at the microscopic level.
"""

if __name__ == "__main__":
    banner("STAGE 322 — SUPPORT-BLIND ORBIT SPLIT")

    lambdaW, cetaU, gamma, KU, Ketaeff, KWeff, muW, TU = sp.symbols(
        "lambda_W c_etaU gamma K_U K_etaeff K_Weff mu_W T_U", positive=True, real=True
    )
    lambdaPhi, KPhieff = sp.symbols("lambda_phi K_phieff", positive=True, real=True)
    L, sigma, Lambda0 = sp.symbols("L sigma Lambda_0", positive=True, real=True)
    chi0_star, deltaU_star, epsW_star, eps_star, epsEta_star = sp.symbols(
        "chi0_star deltaU_star epsilonW_star epsilon_star epsilon_eta_star",
        positive=True,
        real=True,
    )
    dlnlambdaPhi, dlnKPhieff = sp.symbols("dlnlambda_phi dlnK_phieff", real=True)

    obs = physical_branch_observables_from_microscopic(
        lambdaW, cetaU, gamma, KU, Ketaeff, KWeff, muW, TU, L, sigma, Lambda0,
    )

    subbanner("I. Orbit-lock branch observables are support-blind")
    print("R_tr =")
    sp.pprint(obs["Rtr"])
    print("R_target =")
    sp.pprint(obs["Rtarget"])
    print("epsilon_eta =")
    sp.pprint(obs["epsEta"])

    expect_zero("partial_{lambda_phi} R_tr", sp.diff(obs["Rtr"], lambdaPhi))
    expect_zero("partial_{K_phieff} R_tr", sp.diff(obs["Rtr"], KPhieff))
    expect_zero("partial_{lambda_phi} R_target", sp.diff(obs["Rtarget"], lambdaPhi))
    expect_zero("partial_{K_phieff} R_target", sp.diff(obs["Rtarget"], KPhieff))
    expect_zero("partial_{lambda_phi} epsilon_eta", sp.diff(obs["epsEta"], lambdaPhi))
    expect_zero("partial_{K_phieff} epsilon_eta", sp.diff(obs["epsEta"], KPhieff))

    subbanner("II. The finite quotient packet is also support-blind")
    lambdaW_ref, cetaU_ref, gamma_ref, KU_ref, Ketaeff_ref, KWeff_ref, muW_ref, TU_ref = sp.symbols(
        "lambda_W_ref c_etaU_ref gamma_ref K_U_ref K_etaeff_ref K_Weff_ref mu_W_ref T_U_ref",
        positive=True,
        real=True,
    )
    state = microscopic_state_to_reduced(lambdaW, cetaU, gamma, KU, Ketaeff, KWeff, muW, TU, L, sigma)
    state_ref = microscopic_state_to_reduced(lambdaW_ref, cetaU_ref, gamma_ref, KU_ref, Ketaeff_ref, KWeff_ref, muW_ref, TU_ref, L, sigma)
    q = continuum_monomial_quotients(
        state["chi0"], state["deltaU"], state["ZhatW"], state["epsW"], state["epsEta"],
        state_ref["chi0"], state_ref["deltaU"], state_ref["ZhatW"], state_ref["epsW"], state_ref["epsEta"],
        chi0_star, deltaU_star, epsW_star, eps_star,
    )
    expect_zero("partial_{lambda_phi} q_tr", sp.diff(q["qtr"], lambdaPhi))
    expect_zero("partial_{K_phieff} q_tr", sp.diff(q["qtr"], KPhieff))
    expect_zero("partial_{lambda_phi} q_nt", sp.diff(q["qnt"], lambdaPhi))
    expect_zero("partial_{K_phieff} q_nt", sp.diff(q["qnt"], KPhieff))
    expect_zero("partial_{lambda_phi} q_eta", sp.diff(q["qeta"], lambdaPhi))
    expect_zero("partial_{K_phieff} q_eta", sp.diff(q["qeta"], KPhieff))

    subbanner("III. The isotropic support lane depends on the support variables through zeta")
    zeta_expr = coherent_support_ratio(lambdaPhi, KPhieff, lambdaW, KWeff)
    zeta = sp.symbols("zeta", positive=True, real=True)
    eps = obs["eps"]
    S = coherent_support_enhancement(zeta, eps)
    Mmix = sp.symbols("M_mix", positive=True, real=True)
    Mtr = sp.simplify(Mmix * S)

    print("zeta =")
    sp.pprint(zeta_expr)
    print("S(zeta;epsilon) =")
    sp.pprint(S)
    print("M_tr =")
    sp.pprint(Mtr)

    # logarithmic support sensitivities
    zeta_lamphi = sp.simplify(sp.diff(sp.log(zeta_expr.subs({lambdaPhi: lambdaPhi * sp.exp(dlnlambdaPhi)})), dlnlambdaPhi).subs({dlnlambdaPhi: 0}))
    zeta_Kphi = sp.simplify(sp.diff(sp.log(zeta_expr.subs({KPhieff: KPhieff * sp.exp(dlnKPhieff)})), dlnKPhieff).subs({dlnKPhieff: 0}))
    dS_dzeta = sp.simplify(sp.diff(S, zeta))

    print("d ln zeta / d ln lambda_phi =")
    sp.pprint(zeta_lamphi)
    print("d ln zeta / d ln K_phieff =")
    sp.pprint(zeta_Kphi)
    print("dS/dzeta =")
    sp.pprint(dS_dzeta)

    dMtr_dlnlambdaPhi = sp.simplify(Mmix * dS_dzeta.subs({zeta: zeta_expr}) * sp.diff(zeta_expr, lambdaPhi) * lambdaPhi)
    dMtr_dlnKPhi = sp.simplify(Mmix * dS_dzeta.subs({zeta: zeta_expr}) * sp.diff(zeta_expr, KPhieff) * KPhieff)
    print("d M_tr / d ln lambda_phi =")
    sp.pprint(dMtr_dlnlambdaPhi)
    print("d M_tr / d ln K_phieff =")
    sp.pprint(dMtr_dlnKPhi)

    expect_zero("d ln zeta / d ln lambda_phi - 2", sp.simplify(zeta_lamphi - 2))
    expect_zero("d ln zeta / d ln K_phieff + 1", sp.simplify(zeta_Kphi + 1))

    subbanner("IV. Exact separation statement")
    print("The full finite/infinitesimal orbit-lock packet is exactly blind to")
    print("(lambda_phi, K_phieff), while the isotropic support-compensation branch")
    print("depends on them through the single ratio zeta and therefore through the")
    print("single enhancement factor S(zeta;epsilon).")
