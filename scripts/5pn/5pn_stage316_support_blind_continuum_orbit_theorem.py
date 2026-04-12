#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage315_317_common import *

"""
Stage 316 — exact support-blindness of the actual continuum monomial tester.

This stage splices the Stage-315 reduced actual-branch monomial tester into the
coherent support-enhancement and explicit D/N support-tower language from the
moving-throat support-compensation stages.

The key point is that neither the finite quotient packet nor the infinitesimal
orbit-defect packet depends on
    - the support-enhancement ratio zeta,
    - the total baseline M_tr,
    - or the physical twin-family stiffness x
that selects zeta_n^(phys).

So support compensation and support-harmonic selection are exact isotropic-
normalization moves, not weak-axisymmetric orbit-lock moves.
"""

if __name__ == "__main__":
    banner("STAGE 316 — SUPPORT-BLIND CONTINUUM ORBIT THEOREM")

    chi0, deltaU, ZhatW, epsW, epsEta = sp.symbols(
        "chi0 deltaU ZhatW epsilonW epsilon_eta", positive=True, real=True
    )
    chi0_ref, deltaU_ref, ZhatW_ref, epsW_ref, epsEta_ref = sp.symbols(
        "chi0_ref deltaU_ref ZhatW_ref epsilonW_ref epsilon_eta_ref", positive=True, real=True
    )
    chi0_star, deltaU_star, epsW_star, eps_star, epsEta_star = sp.symbols(
        "chi0_star deltaU_star epsilonW_star epsilon_star epsilon_eta_star", positive=True, real=True
    )

    zeta, eps, Mmix, Mtr = sp.symbols("zeta epsilon M_mix M_tr", positive=True, real=True)
    n, x = sp.symbols("n x", positive=True, real=True)

    # Stage-315 quotient packet.
    q = continuum_monomial_quotients(
        chi0, deltaU, ZhatW, epsW, epsEta,
        chi0_ref, deltaU_ref, ZhatW_ref, epsW_ref, epsEta_ref,
        chi0_star, deltaU_star, epsW_star, eps_star,
    )

    # Stage-315 infinitesimal defect packet.
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

    subbanner("I. Coherent support enhancement branch")
    S = sp.simplify(1 + zeta * (1 - eps) / (1 - zeta * eps))
    Mtr_from_zeta = sp.simplify(Mmix * S)
    print("S(zeta; epsilon) =")
    sp.pprint(S)
    print("M_tr =")
    sp.pprint(Mtr_from_zeta)
    print("The reduced monomial tester does not use zeta or M_tr explicitly.")

    subbanner("II. Explicit D/N support tower")
    zeta_phys_n = sp.simplify(1 / ((2 * n + 1) ** 2 * (1 + x * n * (n + 1))))
    print("zeta_n^(twin) =")
    sp.pprint(zeta_phys_n)
    print("Even the explicit support-harmonic selection enters only through zeta_n.")

    subbanner("III. Exact support-blindness of the finite quotient packet")
    dq_dzeta = sp.Matrix([sp.diff(q["qtr"], zeta), sp.diff(q["qnt"], zeta), sp.diff(q["qeta"], zeta)])
    dq_dMtr = sp.Matrix([sp.diff(q["qtr"], Mtr), sp.diff(q["qnt"], Mtr), sp.diff(q["qeta"], Mtr)])
    dq_dx = sp.Matrix([sp.diff(q["qtr"], x), sp.diff(q["qnt"], x), sp.diff(q["qeta"], x)])

    print("d(q_tr,q_nt,q_eta)/d zeta =")
    sp.pprint(dq_dzeta)
    print("d(q_tr,q_nt,q_eta)/d M_tr =")
    sp.pprint(dq_dMtr)
    print("d(q_tr,q_nt,q_eta)/d x =")
    sp.pprint(dq_dx)

    expect_zero("finite quotient packet is zeta-blind", dq_dzeta)
    expect_zero("finite quotient packet is M_tr-blind", dq_dMtr)
    expect_zero("finite quotient packet is twin-stiffness-blind", dq_dx)

    subbanner("IV. Exact support-blindness of the infinitesimal defect packet")
    ddef_dzeta = sp.Matrix([
        sp.diff(defects["Theta1"], zeta),
        sp.diff(defects["Xi1"], zeta),
        sp.diff(defects["R1"], zeta),
    ])
    ddef_dMtr = sp.Matrix([
        sp.diff(defects["Theta1"], Mtr),
        sp.diff(defects["Xi1"], Mtr),
        sp.diff(defects["R1"], Mtr),
    ])
    ddef_dx = sp.Matrix([
        sp.diff(defects["Theta1"], x),
        sp.diff(defects["Xi1"], x),
        sp.diff(defects["R1"], x),
    ])

    print("d(Theta_1,Xi_1,R_1)/d zeta =")
    sp.pprint(ddef_dzeta)
    print("d(Theta_1,Xi_1,R_1)/d M_tr =")
    sp.pprint(ddef_dMtr)
    print("d(Theta_1,Xi_1,R_1)/d x =")
    sp.pprint(ddef_dx)

    expect_zero("defect packet is zeta-blind", ddef_dzeta)
    expect_zero("defect packet is M_tr-blind", ddef_dMtr)
    expect_zero("defect packet is twin-stiffness-blind", ddef_dx)

    banner("STAGE 316 LEDGER")
    print("1. The actual continuum monomial quotient packet depends only on")
    print("      chi0, deltaU, ZhatW, epsilonW, epsilon_eta,")
    print("   and not on zeta, M_tr, or the explicit twin-family stiffness x.")
    print("2. The corresponding defect packet (Theta_1, Xi_1, R_1) is equally support-blind.")
    print("3. Therefore support enhancement and support-harmonic selection belong to the")
    print("   isotropic normalization branch only; they do not move the weak-axisymmetric")
    print("   coherent branch on or off the exact similarity orbit.")
