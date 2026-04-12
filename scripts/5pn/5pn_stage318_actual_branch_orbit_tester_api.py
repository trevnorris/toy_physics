#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage315_317_common import *

"""
Stage 318 — packaged actual-branch orbit tester API.

This is the first ready-to-use branch tester assembled out of Stages 315–317.
It records exactly what data a future moving-throat PDE solve has to supply:

Finite branch state:
    (chi0, deltaU, ZhatW, epsilonW, epsilon_eta)
Reference branch state:
    (chi0_ref, deltaU_ref, ZhatW_ref, epsilonW_ref, epsilon_eta_ref)
Reference orbit exponents:
    (chi0_star, deltaU_star, epsilonW_star, epsilon_star)

From those data the tester returns either
    (q_tr, q_nt, q_eta)
for finite placement, or
    (Theta_1, Xi_1, R_1)
for linearized weak-axisymmetric defects.
"""

if __name__ == "__main__":
    banner("STAGE 318 — ACTUAL-BRANCH ORBIT TESTER API")

    chi0, deltaU, ZhatW, epsW, epsEta = sp.symbols(
        "chi0 deltaU ZhatW epsilonW epsilon_eta", positive=True, real=True
    )
    chi0_ref, deltaU_ref, ZhatW_ref, epsW_ref, epsEta_ref = sp.symbols(
        "chi0_ref deltaU_ref ZhatW_ref epsilonW_ref epsilon_eta_ref", positive=True, real=True
    )
    chi0_star, deltaU_star, epsW_star, eps_star, epsEta_star = sp.symbols(
        "chi0_star deltaU_star epsilonW_star epsilon_star epsilon_eta_star",
        positive=True,
        real=True,
    )

    subbanner("I. Finite actual-branch state -> quotient packet")
    q = continuum_monomial_quotients(
        chi0, deltaU, ZhatW, epsW, epsEta,
        chi0_ref, deltaU_ref, ZhatW_ref, epsW_ref, epsEta_ref,
        chi0_star, deltaU_star, epsW_star, eps_star,
    )
    print("q_tr =")
    sp.pprint(q["qtr"])
    print("q_nt =")
    sp.pprint(q["qnt"])
    print("q_eta =")
    sp.pprint(q["qeta"])

    subbanner("II. Infinitesimal drift -> defect packet")
    dln_chi0, dln_deltaU, dln_ZhatW, dln_epsW, dln_epsEta = sp.symbols(
        "dlnchi0 dlndeltaU dlnZhatW dlnepsilonW dlnepsilon_eta", real=True
    )
    sigmas = continuum_monomial_drifts(
        dln_chi0, dln_deltaU, dln_ZhatW, dln_epsW, dln_epsEta,
        chi0_star, deltaU_star, epsW_star, eps_star,
    )
    defects = defect_map_from_sigmas(
        sigmas["Sigma_tr"], sigmas["Sigma_nt"], sigmas["Sigma_eta"],
        chi0_star, deltaU_star, epsEta_star,
    )
    print("Theta_1 =")
    sp.pprint(defects["Theta1"])
    print("Xi_1 =")
    sp.pprint(defects["Xi1"])
    print("R_1 =")
    sp.pprint(defects["R1"])

    subbanner("III. Exact zero-defect test packet")
    a_free, b_free = sp.symbols("a_free b_free", real=True)
    basis = reduced_zero_defect_basis(chi0_star, deltaU_star, epsW_star, eps_star)
    drift_zero = sp.simplify(a_free * basis["v_chi"] + b_free * basis["v_epsW"])
    print("Generic zero-defect reduced drift vector =")
    sp.pprint(drift_zero)

    sigmas_zero = continuum_monomial_drifts(
        drift_zero[0], drift_zero[1], drift_zero[2], drift_zero[3], drift_zero[4],
        chi0_star, deltaU_star, epsW_star, eps_star,
    )
    defects_zero = defect_map_from_sigmas(
        sigmas_zero["Sigma_tr"], sigmas_zero["Sigma_nt"], sigmas_zero["Sigma_eta"],
        chi0_star, deltaU_star, epsEta_star,
    )
    expect_zero("Theta_1 on the exact reduced orbit plane", defects_zero["Theta1"])
    expect_zero("Xi_1 on the exact reduced orbit plane", defects_zero["Xi1"])
    expect_zero("R_1 on the exact reduced orbit plane", defects_zero["R1"])

    banner("STAGE 318 LEDGER")
    print("1. The finite actual-branch orbit tester needs only the five-state packet")
    print("      (chi0, deltaU, ZhatW, epsilonW, epsilon_eta),")
    print("   plus a reference state and the reference orbit exponents.")
    print("2. The infinitesimal tester acts on the five-drift packet")
    print("      (dln chi0, dln deltaU, dln ZhatW, dln epsilonW, dln epsilon_eta)")
    print("   and returns either the monomial drift triple")
    print("      (Sigma_tr, Sigma_nt, Sigma_eta)")
    print("   or the physical defect triple")
    print("      (Theta_1, Xi_1, R_1).")
    print("3. Exact orbit lock in the reduced continuum state space is equivalent to")
    print("   membership in the 2-plane spanned by the two basis vectors from Stage 317.")
