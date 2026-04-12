#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage323_325_common import *

"""
Stage 325 — exact two-packet compiler on the coherent local branch.

This stage packages the actual coherent branch into two rigorously separated packets:

1. Orbit-lock packet (finite or infinitesimal), compiled from
       (chi0, deltaU, Z_W, epsilon_W, epsilon_eta, Lambda).
2. Support/normalization packet, compiled from
       (chi0, deltaU, Z_W, epsilon_W, epsilon_eta, Lambda, zeta).

The exact separation is that zeta enters only the support-enhancement factor S(zeta;epsilon)
through the total baseline M_tr = M_mix S, while the orbit-lock packet is exactly blind to zeta.

So the physical coherent branch is already split into an orbit side and a support side before the
completed moving-throat PDE is solved numerically.
"""

if __name__ == "__main__":
    banner("STAGE 325 — COHERENT BRANCH TWO-PACKET COMPILER")

    chi0, deltaU, ZW, epsW, epsEta, Lambda, zeta = sp.symbols(
        "chi0 deltaU Z_W epsilon_W epsilon_eta Lambda zeta", positive=True, real=True
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
    supp = support_packet(chi0, deltaU, ZW, epsW, epsEta, Lambda, zeta)

    subbanner("I. Orbit-lock packet")
    print("Finite quotient packet =")
    sp.pprint(sp.Matrix([q["qtr"], q["qnt"], q["qeta"]]))
    print("Infinitesimal defect packet =")
    sp.pprint(sp.Matrix([defects["Theta1"], defects["Xi1"], defects["R1"]]))

    subbanner("II. Support/normalization packet")
    print("epsilon =")
    sp.pprint(supp["eps"])
    print("M_mix =")
    sp.pprint(supp["Mmix"])
    print("S(zeta;epsilon) =")
    sp.pprint(supp["S"])
    print("M_tr =")
    sp.pprint(supp["Mtr"])
    print("R_target =")
    sp.pprint(supp["Rtarget"])
    print("R_target * M_mix =")
    sp.pprint(supp["product"])

    subbanner("III. Exact mixed-only product law")
    expect_zero(
        "R_target * M_mix - 8 Lambda (1-epsilon)/pi^2",
        sp.simplify(supp["product"] - 8 * Lambda * (1 - supp["eps"]) / sp.pi**2),
    )

    subbanner("IV. Exact zeta separation")
    expect_zero("partial_zeta q_tr", sp.diff(q["qtr"], zeta))
    expect_zero("partial_zeta q_nt", sp.diff(q["qnt"], zeta))
    expect_zero("partial_zeta q_eta", sp.diff(q["qeta"], zeta))
    expect_zero("partial_zeta Theta_1", sp.diff(defects["Theta1"], zeta))
    expect_zero("partial_zeta Xi_1", sp.diff(defects["Xi1"], zeta))
    expect_zero("partial_zeta R_1", sp.diff(defects["R1"], zeta))
    expect_zero("partial_zeta M_mix", sp.diff(supp["Mmix"], zeta))
    expect_zero("partial_zeta R_target", sp.diff(supp["Rtarget"], zeta))
    dMtr_dzeta = sp.simplify(sp.diff(supp["Mtr"], zeta))
    print("d M_tr / d zeta =")
    sp.pprint(dMtr_dzeta)

    subbanner("V. Exact physical interpretation")
    print("The coherent local branch is already split into two exact packets.")
    print()
    print("Orbit side:")
    print("  depends only on (chi0, deltaU, Z_W, epsilon_W, epsilon_eta, Lambda)")
    print("  through (q_tr,q_nt,q_eta) or equivalently (Theta_1, Xi_1, R_1).")
    print()
    print("Support side:")
    print("  depends on the same six placement variables plus the single extra support ratio zeta,")
    print("  and zeta enters only through the support-enhancement factor S(zeta;epsilon).")
    print()
    print("So support compensation cannot rescue orbit-lock failure, and orbit lock does not")
    print("determine support enhancement.  The completed moving-throat PDE must satisfy the")
    print("two packets separately on the same physical branch.")
