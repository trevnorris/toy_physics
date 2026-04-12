#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage329_332_common import *

"""
Stage 330 — actual coherent local D/N orbit packet.

What this stage does
--------------------
1. Compiles the finite quotient packet directly on the actual coherent local D/N branch
   in the exact branch observables (R_tr, R_target, epsilon_eta).
2. Compiles the infinitesimal defect packet on the same branch.
3. Verifies that the orbit packet is blind to the coherent support ratio zeta and
   therefore blind to the explicit D/N support-family choice.
"""

if __name__ == "__main__":
    banner("STAGE 330 — ACTUAL COHERENT LOCAL D/N ORBIT PACKET")

    chi0, deltaU, ZW, epsW, epsEta, Lambda, zeta = sp.symbols(
        "chi0 deltaU Z_W epsilon_W epsilon_eta Lambda zeta", positive=True, real=True
    )
    chi0_ref, deltaU_ref, ZW_ref, epsW_ref, epsEta_ref, Lambda_ref = sp.symbols(
        "chi0_ref deltaU_ref Z_W_ref epsilon_W_ref epsilon_eta_ref Lambda_ref",
        positive=True,
        real=True,
    )
    chi0_star, deltaU_star, epsEta_star = sp.symbols(
        "chi0_star deltaU_star epsilon_eta_star", positive=True, real=True
    )
    dln_chi0, dln_deltaU, dln_ZW, dln_epsW, dln_epsEta, dln_Lambda = sp.symbols(
        "dlnchi0 dlndeltaU dlnZ_W dlnepsilon_W dlnepsilon_eta dlnLambda", real=True
    )

    branch = coherent_tracking_state(chi0, deltaU, ZW, epsW, epsEta, Lambda, zeta)
    branch_ref = coherent_tracking_state(chi0_ref, deltaU_ref, ZW_ref, epsW_ref, epsEta_ref, Lambda_ref, zeta)

    q = direct_branch_packet(
        branch["Rtr"], branch["Rtarget"], epsEta,
        branch_ref["Rtr"], branch_ref["Rtarget"], epsEta_ref,
        chi0_star, deltaU_star,
    )

    var_to_dln = {
        chi0: dln_chi0,
        deltaU: dln_deltaU,
        ZW: dln_ZW,
        epsW: dln_epsW,
        epsEta: dln_epsEta,
        Lambda: dln_Lambda,
    }
    dln_Rtr = log_drift(branch["Rtr"], {chi0: dln_chi0, deltaU: dln_deltaU})
    dln_Rtarget = log_drift(branch["Rtarget"], var_to_dln)
    defects = direct_branch_defects(branch["Rtr"], branch["Rtarget"], epsEta, dln_Rtr, dln_Rtarget, dln_epsEta, epsEta_star)

    subbanner("I. Exact finite orbit packet on the actual coherent local D/N branch")
    print("R_tr =")
    sp.pprint(branch["Rtr"])
    print("R_target =")
    sp.pprint(branch["Rtarget"])
    print("q_tr =")
    sp.pprint(q["qtr"])
    print("q_nt =")
    sp.pprint(q["qnt"])
    print("q_eta =")
    sp.pprint(q["qeta"])

    subbanner("II. Exact infinitesimal defect packet on the same branch")
    print("d ln R_tr =")
    sp.pprint(dln_Rtr)
    print("d ln R_target =")
    sp.pprint(dln_Rtarget)
    print("Theta_1 =")
    sp.pprint(defects["Theta1"])
    print("Xi_1 =")
    sp.pprint(defects["Xi1"])
    print("R_1 =")
    sp.pprint(defects["R1"])
    expect_zero("Theta_1 - dln R_tr", sp.simplify(defects["Theta1"] - dln_Rtr))
    expect_zero(
        "Xi_1 + dln R_target + epsEta_* dln epsEta/(1-epsEta_*)",
        sp.simplify(defects["Xi1"] + dln_Rtarget + epsEta_star * dln_epsEta / (1 - epsEta_star)),
    )
    expect_zero("R_1 - dln R_target", sp.simplify(defects["R1"] - dln_Rtarget))

    subbanner("III. Exact support-blindness")
    expect_zero("partial_zeta R_tr", sp.diff(branch["Rtr"], zeta))
    expect_zero("partial_zeta R_target", sp.diff(branch["Rtarget"], zeta))
    expect_zero("partial_zeta q_tr", sp.diff(q["qtr"], zeta))
    expect_zero("partial_zeta q_nt", sp.diff(q["qnt"], zeta))
    expect_zero("partial_zeta q_eta", sp.diff(q["qeta"], zeta))
    expect_zero("partial_zeta Theta_1", sp.diff(defects["Theta1"], zeta))
    expect_zero("partial_zeta Xi_1", sp.diff(defects["Xi1"], zeta))
    expect_zero("partial_zeta R_1", sp.diff(defects["R1"], zeta))

    subbanner("IV. Interpretation")
    print("On the actual coherent local D/N branch, the orbit-lock packet is already")
    print("entirely carried by the three exact branch observables")
    print("  (R_tr, R_target, epsilon_eta),")
    print("or equivalently by the finite quotient packet (q_tr, q_nt, q_eta).")
    print()
    print("The coherent support lane does not enter this packet at all. So choosing the")
    print("physical D/N support harmonic changes the support theorem, but it does not")
    print("move the branch on or off the weak-axisymmetric similarity orbit.")
