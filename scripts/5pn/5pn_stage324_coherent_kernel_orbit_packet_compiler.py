#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage323_325_common import *

"""
Stage 324 — exact orbit-packet compiler in coherent placement variables.

This stage compiles the finite and infinitesimal orbit-lock packet directly from the
coherent local-kernel placement variables

    (chi0, deltaU, Z_W, epsilon_W, epsilon_eta, Lambda),

rather than from the full microscopic state or the intermediate reduced variable ZhatW.

The exact bridge is

    ZhatW = Z_W * Lambda0 / Lambda,
    dln ZhatW = dln Z_W - dln Lambda.

So the whole weak-axisymmetric tester now acts directly on the physical placement map.
"""

if __name__ == "__main__":
    banner("STAGE 324 — COHERENT KERNEL ORBIT-PACKET COMPILER")

    chi0, deltaU, ZW, epsW, epsEta, Lambda = sp.symbols(
        "chi0 deltaU Z_W epsilon_W epsilon_eta Lambda", positive=True, real=True
    )
    chi0_ref, deltaU_ref, ZW_ref, epsW_ref, epsEta_ref, Lambda_ref = sp.symbols(
        "chi0_ref deltaU_ref Z_W_ref epsilon_W_ref epsilon_eta_ref Lambda_ref",
        positive=True,
        real=True,
    )
    chi0_star, deltaU_star, epsW_star, epsEta_star = sp.symbols(
        "chi0_star deltaU_star epsilon_W_star epsilon_eta_star",
        positive=True,
        real=True,
    )
    eps_star = split_epsilon(epsW_star, deltaU_star)
    dln_chi0, dln_deltaU, dln_ZW, dln_epsW, dln_epsEta, dln_Lambda = sp.symbols(
        "dlnchi0 dlndeltaU dlnZ_W dlnepsilon_W dlnepsilon_eta dlnLambda",
        real=True,
    )

    q = placement_finite_packet(
        chi0, deltaU, ZW, epsW, epsEta, Lambda,
        chi0_ref, deltaU_ref, ZW_ref, epsW_ref, epsEta_ref, Lambda_ref,
        chi0_star, deltaU_star, epsW_star, eps_star,
    )

    sigmas = placement_monomial_drifts(
        dln_chi0, dln_deltaU, dln_ZW, dln_epsW, dln_epsEta, dln_Lambda,
        chi0_star, deltaU_star, epsW_star, eps_star,
    )
    defects = defect_map_from_sigmas(
        sigmas["Sigma_tr"], sigmas["Sigma_nt"], sigmas["Sigma_eta"],
        chi0_star, deltaU_star, epsEta_star,
    )

    obs = coherent_branch_observables(chi0, deltaU, ZW, epsW, epsEta, Lambda)
    obs_star = coherent_branch_observables(chi0_star, deltaU_star, sp.Symbol("Z_W_star", positive=True), epsW_star, epsEta_star, sp.Symbol("Lambda_star", positive=True))

    subbanner("I. Exact finite quotient packet in coherent placement variables")
    print("q_tr =")
    sp.pprint(q["qtr"])
    print("q_nt =")
    sp.pprint(q["qnt"])
    print("q_eta =")
    sp.pprint(q["qeta"])

    subbanner("II. Exact infinitesimal monomial-drift packet")
    print("Sigma_tr =")
    sp.pprint(sigmas["Sigma_tr"])
    print("Sigma_nt =")
    sp.pprint(sigmas["Sigma_nt"])
    print("Sigma_eta =")
    sp.pprint(sigmas["Sigma_eta"])

    subbanner("III. Exact physical defect packet")
    print("Theta_1 =")
    sp.pprint(defects["Theta1"])
    print("Xi_1 =")
    sp.pprint(defects["Xi1"])
    print("R_1 =")
    sp.pprint(defects["R1"])

    subbanner("IV. Exact direct-observable drift bridge")
    var_to_dln = {
        chi0: dln_chi0,
        deltaU: dln_deltaU,
        ZW: dln_ZW,
        epsW: dln_epsW,
        epsEta: dln_epsEta,
        Lambda: dln_Lambda,
    }
    dln_Rtr = log_drift(obs["Rtr"], {chi0: dln_chi0, deltaU: dln_deltaU}).subs({chi0: chi0_star, deltaU: deltaU_star})
    dln_Rtarget = log_drift(obs["Rtarget"], var_to_dln).subs({chi0: chi0_star, deltaU: deltaU_star, epsW: epsW_star, epsEta: epsEta_star})

    print("d ln R_tr =")
    sp.pprint(dln_Rtr)
    print("d ln R_target =")
    sp.pprint(dln_Rtarget)

    expect_zero("Theta_1 - dln R_tr", sp.simplify(defects["Theta1"] - dln_Rtr))
    expect_zero(
        "Xi_1 + dln R_target + epsEta_* dln epsEta/(1-epsEta_*)",
        sp.simplify(defects["Xi1"] + dln_Rtarget + epsEta_star * dln_epsEta / (1 - epsEta_star)),
    )
    expect_zero("R_1 - dln R_target", sp.simplify(defects["R1"] - dln_Rtarget))

    subbanner("V. Reduced exact orbit-lock conditions on the coherent placement map")
    print("Exact weak-axisymmetric orbit lock on the coherent branch is equivalent to")
    print("  dln R_tr = 0,")
    print("  dln R_target = 0,")
    print("  dln epsilon_eta = 0.")
    print()
    print("In placement variables this is the same as")
    print("  (1+deltaU_*) dln chi0 + (1+chi0_*) dln deltaU = 0,")
    print("  dln Z_W - dln Lambda + E_* dln epsilon_W - F_* dln deltaU = 0,")
    print("  dln epsilon_eta = 0.")
