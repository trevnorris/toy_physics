#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage326_328_common import *
from fivepn_stage319_322_common import microscopic_reduced_drifts

"""
Stage 327 — exact microscopic coherent placement-drift extractor.

This stage extends the older microscopic -> reduced drift extractor to the *actual*
coherent placement variables used by the Stage 323–325 branch compilers.

The new drift variables are
    dln Z_W, dln Lambda, dln zeta,
so the actual coherent branch can now be tested directly in the language of the
physical placement map rather than through the older reduced variable ZhatW.
"""

if __name__ == "__main__":
    banner("STAGE 327 — MICROSCOPIC COHERENT PLACEMENT-DRIFT EXTRACTOR")

    lambda1, c1, gamma1 = sp.symbols("lambda_1 c_1 gamma_1", real=True)
    kappaU, kappaEta, kappaW = sp.symbols("kappa_U kappa_eta kappa_W", real=True)
    mu1, tau1 = sp.symbols("mu_1 tau_1", real=True)
    phi1, kappaPhi = sp.symbols("phi_1 kappa_phi", real=True)

    drifts = microscopic_coherent_placement_drifts(
        lambda1, c1, gamma1, kappaU, kappaEta, kappaW, mu1, tau1, phi1, kappaPhi
    )
    drifts_old = microscopic_reduced_drifts(lambda1, c1, gamma1, kappaU, kappaEta, kappaW, mu1, tau1)

    subbanner("I. Exact microscopic -> actual placement drifts")
    print("dln chi0 =")
    sp.pprint(drifts["dln_chi0"])
    print("dln deltaU =")
    sp.pprint(drifts["dln_deltaU"])
    print("dln Z_W =")
    sp.pprint(drifts["dln_ZW"])
    print("dln epsilon_W =")
    sp.pprint(drifts["dln_epsW"])
    print("dln epsilon_eta =")
    sp.pprint(drifts["dln_epsEta"])
    print("dln Lambda =")
    sp.pprint(drifts["dln_Lambda"])
    print("dln zeta =")
    sp.pprint(drifts["dln_zeta"])

    subbanner("II. Exact bridge back to the older reduced drift packet")
    print("dln ZhatW from the actual placement drifts =")
    sp.pprint(drifts["dln_ZhatW"])
    print("dln ZhatW from the older reduced extractor =")
    sp.pprint(drifts_old["dln_ZhatW"])
    expect_zero("dln ZhatW equivalence", drifts["dln_ZhatW"] - drifts_old["dln_ZhatW"])
    expect_zero("dln chi0 equivalence", drifts["dln_chi0"] - drifts_old["dln_chi0"])
    expect_zero("dln deltaU equivalence", drifts["dln_deltaU"] - drifts_old["dln_deltaU"])
    expect_zero("dln epsilon_W equivalence", drifts["dln_epsW"] - drifts_old["dln_epsW"])
    expect_zero("dln epsilon_eta equivalence", drifts["dln_epsEta"] - drifts_old["dln_epsEta"])

    subbanner("III. Exact support / orbit split at drift level")
    expect_zero("partial_{phi1} dln chi0", sp.diff(drifts["dln_chi0"], phi1))
    expect_zero("partial_{phi1} dln deltaU", sp.diff(drifts["dln_deltaU"], phi1))
    expect_zero("partial_{phi1} dln Z_W", sp.diff(drifts["dln_ZW"], phi1))
    expect_zero("partial_{phi1} dln epsilon_W", sp.diff(drifts["dln_epsW"], phi1))
    expect_zero("partial_{phi1} dln epsilon_eta", sp.diff(drifts["dln_epsEta"], phi1))
    expect_zero("partial_{phi1} dln Lambda", sp.diff(drifts["dln_Lambda"], phi1))

    print("Only dln zeta depends on the explicit support-lane drifts:")
    print("dln zeta =")
    sp.pprint(drifts["dln_zeta"])

    subbanner("IV. Interpretation")
    print("The actual coherent placement drift packet is now completely explicit in the")
    print("microscopic grouped weak-axisymmetric drifts.  The orbit side uses")
    print("  (dln chi0, dln deltaU, dln Z_W, dln epsilon_W, dln epsilon_eta, dln Lambda),")
    print("while the support side adds only the separate support drift dln zeta.")
    print("So support transport cannot hide inside the orbit packet: it is an exactly")
    print("separate logarithmic direction already at the microscopic-kernel level.")
