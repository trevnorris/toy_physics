#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage319_322_common import *

"""
Stage 319 — exact microscopic -> reduced actual-branch drift extractor.

This stage bridges the Stage-168 microscopic grouped weak-axisymmetric drifts
    (lambda1, c1, gamma1, kappaU, kappaEta, kappaW, mu1, tau1)
into the reduced Stage-315/317 continuum drift packet
    (dln chi0, dln deltaU, dln ZhatW, dln epsilonW, dln epsilon_eta).

The payoff is that the Stage-318 actual-branch orbit tester can now be fed
*directly* from microscopic kernel drifts rather than from an abstract reduced
five-drift packet.
"""

if __name__ == "__main__":
    banner("STAGE 319 — MICROSCOPIC REDUCED DRIFT EXTRACTOR")

    lambda1, c1, gamma1, kappaU, kappaEta, kappaW, mu1, tau1 = sp.symbols(
        "lambda_1 c_1 gamma_1 kappa_U kappa_eta kappa_W mu_1 tau_1", real=True
    )
    mic = sp.Matrix([lambda1, c1, gamma1, kappaU, kappaEta, kappaW, mu1, tau1])

    red = microscopic_reduced_drifts(lambda1, c1, gamma1, kappaU, kappaEta, kappaW, mu1, tau1)
    redvec = sp.Matrix([
        red["dln_chi0"],
        red["dln_deltaU"],
        red["dln_ZhatW"],
        red["dln_epsW"],
        red["dln_epsEta"],
    ])

    subbanner("I. Exact microscopic -> reduced drift formulas")
    print("dln chi0 =")
    sp.pprint(red["dln_chi0"])
    print("dln deltaU =")
    sp.pprint(red["dln_deltaU"])
    print("dln ZhatW =")
    sp.pprint(red["dln_ZhatW"])
    print("dln epsilonW =")
    sp.pprint(red["dln_epsW"])
    print("dln epsilon_eta =")
    sp.pprint(red["dln_epsEta"])

    subbanner("II. Exact extractor matrix")
    X = sp.Matrix([
        [sp.diff(redvec[i], mic[j]) for j in range(len(mic))]
        for i in range(len(redvec))
    ])
    print("X_micro->red =")
    sp.pprint(X)
    print("rank(X_micro->red) =", X.rank())

    dlnchi0_expected = gamma1 + c1 - kappaU
    dlndeltaU_expected = tau1 - kappaU
    dlnZhatW_expected = 2 * lambda1 + mu1 - kappaEta - 2 * kappaW
    dlnepsW_expected = 2 * gamma1 + 2 * lambda1 - kappaU - kappaW
    dlnepsEta_expected = 2 * c1 - kappaU - kappaEta

    expect_zero("dln chi0 match", red["dln_chi0"] - dlnchi0_expected)
    expect_zero("dln deltaU match", red["dln_deltaU"] - dlndeltaU_expected)
    expect_zero("dln ZhatW match", red["dln_ZhatW"] - dlnZhatW_expected)
    expect_zero("dln epsilonW match", red["dln_epsW"] - dlnepsW_expected)
    expect_zero("dln epsilon_eta match", red["dln_epsEta"] - dlnepsEta_expected)

    subbanner("III. Interpretation")
    print("The reduced Stage-318 tester inputs are now exact linear images of the")
    print("microscopic grouped weak-axisymmetric drift vector.  No extra reduced")
    print("closure assumptions are needed to pass from the eight microscopic drifts")
    print("to the five actual-branch drifts.")
