#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage315_317_common import *

"""
Stage 315 — exact actual-branch monomial drift map in reduced continuum variables.

This stage rewrites the Stage-312 direct microscopic monomials in the smallest
actual coherent-branch variable set that a moving-throat branch extractor would
naturally produce:
    (chi0, deltaU, ZhatW, epsilonW, epsilon_eta),
where
    ZhatW := lambda_W^2 mu_W / (K_etaeff K_Weff^2) = Z_W / Omega_W^2.

The payoff is that the actual-branch orbit-lock tester no longer needs the full
8-variable microscopic kernel state.  It needs only five branch observables or,
at infinitesimal level, five logarithmic drifts.
"""

if __name__ == "__main__":
    banner("STAGE 315 — CONTINUUM MONOMIAL DRIFT MAP")

    chi0, deltaU, ZhatW, epsW, epsEta = sp.symbols(
        "chi0 deltaU ZhatW epsilonW epsilon_eta", positive=True, real=True
    )
    chi0_ref, deltaU_ref, ZhatW_ref, epsW_ref, epsEta_ref = sp.symbols(
        "chi0_ref deltaU_ref ZhatW_ref epsilonW_ref epsilon_eta_ref", positive=True, real=True
    )
    chi0_star, deltaU_star, epsW_star, eps_star = sp.symbols(
        "chi0_star deltaU_star epsilonW_star epsilon_star", positive=True, real=True
    )

    mono = continuum_monomials(
        chi0, deltaU, ZhatW, epsW, epsEta,
        chi0_star, deltaU_star, epsW_star, eps_star,
    )
    q = continuum_monomial_quotients(
        chi0, deltaU, ZhatW, epsW, epsEta,
        chi0_ref, deltaU_ref, ZhatW_ref, epsW_ref, epsEta_ref,
        chi0_star, deltaU_star, epsW_star, eps_star,
    )

    subbanner("I. Exact direct monomials in reduced continuum variables")
    print("C_tr,* =")
    sp.pprint(mono["Ctr"])
    print("C_nt,* =")
    sp.pprint(mono["Cnt"])
    print("epsilon_eta =")
    sp.pprint(mono["epsEta"])

    subbanner("II. Exact finite quotient coordinates")
    print("q_tr =")
    sp.pprint(q["qtr"])
    print("q_nt =")
    sp.pprint(q["qnt"])
    print("q_eta =")
    sp.pprint(q["qeta"])

    qtr_expected = sp.simplify(
        (1 + deltaU_star) * sp.log(chi0 / chi0_ref)
        + (1 + chi0_star) * sp.log(deltaU / deltaU_ref)
    )
    qnt_expected = sp.simplify(
        sp.log(ZhatW / ZhatW_ref)
        + q["E_star"] * sp.log(epsW / epsW_ref)
        - q["F_star"] * sp.log(deltaU / deltaU_ref)
    )
    qeta_expected = sp.simplify(sp.log(epsEta / epsEta_ref))

    expect_zero("q_tr reduced-continuum form", sp.simplify(sp.expand(sp.expand_log(q["qtr"], force=True) - sp.expand_log(qtr_expected, force=True))))
    expect_zero("q_nt reduced-continuum form", sp.simplify(sp.expand(sp.expand_log(q["qnt"], force=True) - sp.expand_log(qnt_expected, force=True))))
    expect_zero("q_eta reduced-continuum form", sp.simplify(sp.expand(sp.expand_log(q["qeta"], force=True) - sp.expand_log(qeta_expected, force=True))))

    subbanner("III. Exact infinitesimal drift map")
    dln_chi0, dln_deltaU, dln_ZhatW, dln_epsW, dln_epsEta = sp.symbols(
        "dlnchi0 dlndeltaU dlnZhatW dlnepsilonW dlnepsilon_eta", real=True
    )
    sigmas = continuum_monomial_drifts(
        dln_chi0, dln_deltaU, dln_ZhatW, dln_epsW, dln_epsEta,
        chi0_star, deltaU_star, epsW_star, eps_star,
    )
    print("Sigma_tr =")
    sp.pprint(sigmas["Sigma_tr"])
    print("Sigma_nt =")
    sp.pprint(sigmas["Sigma_nt"])
    print("Sigma_eta =")
    sp.pprint(sigmas["Sigma_eta"])

    subbanner("IV. Exact actual-branch tester inputs")
    print("The actual orbit-lock tester now needs only the reduced branch observables")
    print("  chi0, deltaU, ZhatW, epsilonW, epsilon_eta")
    print("or, infinitesimally, their logarithmic drifts.")
    print()
    print("The five-drift packet maps exactly to the three branch-adapted monomial drifts")
    print("  (Sigma_tr, Sigma_nt, Sigma_eta).")
