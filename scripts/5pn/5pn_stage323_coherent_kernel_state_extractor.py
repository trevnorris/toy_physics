#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage323_325_common import *

"""
Stage 323 — coherent local-kernel state extractor.

This stage pushes the Stage-319/321 microscopic state down to the actual coherent
placement variables used in the moving-throat compact program:

    (chi0, deltaU, Z_W, epsilon_W, epsilon_eta, Lambda, zeta).

The new point is that the reduced Stage-318 state can already be reconstructed from
this placement map without reopening the full microscopic 8-tuple.  The only extra
constant needed is

    Lambda0 = 27*pi^2*G*c_s^5 / (20*a^5*c^5),

because

    ZhatW = Z_W / Omega_W^2 = Z_W * Lambda0 / Lambda.

So the actual coherent branch state needed by the orbit-lock packet is now expressed
in the same language as the support-compensation stages.
"""

if __name__ == "__main__":
    banner("STAGE 323 — COHERENT KERNEL STATE EXTRACTOR")

    chi0, deltaU, ZW, epsW, epsEta, Lambda, zeta = sp.symbols(
        "chi0 deltaU Z_W epsilon_W epsilon_eta Lambda zeta", positive=True, real=True
    )
    G, c, cs, a = sp.symbols("G c c_s a", positive=True, real=True)
    Lambda0 = lambda0_constant(G, cs, a, c)

    red = placement_to_reduced_state(ZW, Lambda, Lambda0)
    obs = coherent_branch_observables(chi0, deltaU, ZW, epsW, epsEta, Lambda)

    subbanner("I. Exact coherent placement variables")
    print("chi0 =")
    sp.pprint(chi0)
    print("deltaU =")
    sp.pprint(deltaU)
    print("Z_W =")
    sp.pprint(ZW)
    print("epsilon_W =")
    sp.pprint(epsW)
    print("epsilon_eta =")
    sp.pprint(epsEta)
    print("Lambda =")
    sp.pprint(Lambda)
    print("zeta =")
    sp.pprint(zeta)

    subbanner("II. Exact reduced Stage-318 state extracted from the coherent placement map")
    print("Lambda0 =")
    sp.pprint(Lambda0)
    print("ZhatW =")
    sp.pprint(red["ZhatW"])
    print("epsilon =")
    sp.pprint(obs["eps"])
    print("R_tr =")
    sp.pprint(obs["Rtr"])
    print("R_target =")
    sp.pprint(obs["Rtarget"])

    subbanner("III. Exact equivalence of the two R_target formulas")
    Rtarget_reduced = sp.simplify(
        Lambda0 * (1 - epsEta) * (1 - obs["eps"]) ** 2 / (red["ZhatW"] * (1 + chi0) ** 2)
    )
    print("R_target via placement map =")
    sp.pprint(obs["Rtarget"])
    print("R_target via reduced-state map =")
    sp.pprint(Rtarget_reduced)
    expect_zero("R_target equivalence", sp.simplify(obs["Rtarget"] - Rtarget_reduced))

    subbanner("IV. Support-blindness of the extracted orbit state")
    expect_zero("partial_{zeta} ZhatW", sp.diff(red["ZhatW"], zeta))
    expect_zero("partial_{zeta} epsilon", sp.diff(obs["eps"], zeta))
    expect_zero("partial_{zeta} R_tr", sp.diff(obs["Rtr"], zeta))
    expect_zero("partial_{zeta} R_target", sp.diff(obs["Rtarget"], zeta))

    subbanner("V. Interpretation")
    print("The actual coherent local-kernel placement map already determines the reduced")
    print("Stage-318 state used by the orbit-lock packet.  The support ratio zeta is absent")
    print("from that extracted state; it survives only on the separate support-compensation")
    print("lane.  So the actual branch data needed by the weak-axisymmetric tester are")
    print("now expressible directly in the coherent-kernel variables")
    print("  (chi0, deltaU, Z_W, epsilon_W, epsilon_eta, Lambda)")
    print("without reopening the full microscopic kernel tuple.")
