#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage329_332_common import *

"""
Stage 332 — actual coherent local D/N branch theorem gate.

What this stage does
--------------------
1. Packages the actual coherent local D/N branch into one exact two-packet theorem gate.
2. Shows that the weak-axisymmetric orbit-lock test and the coherent support theorem
   live on the same physical branch but remain rigorously separated.
3. Specializes the support side to the first physical support lane, the lowest symmetric
   twin branch, and records the exact regime split.
"""

if __name__ == "__main__":
    banner("STAGE 332 — ACTUAL COHERENT LOCAL D/N BRANCH THEOREM GATE")

    # Actual branch state
    chi0, deltaU, ZW, epsW, epsEta, Lambda = sp.symbols(
        "chi0 deltaU Z_W epsilon_W epsilon_eta Lambda", positive=True, real=True
    )
    zeta = sp.symbols("zeta", nonnegative=True, real=True)
    xi, delta, R = sp.symbols("xi delta R", positive=True, real=True)

    # Reference data for the finite orbit packet
    Rtr_ref, Rtarget_ref, epsEta_ref = sp.symbols(
        "R_tr_ref R_target_ref epsilon_eta_ref", positive=True, real=True
    )
    chi0_star, deltaU_star = sp.symbols("chi0_star deltaU_star", positive=True, real=True)

    branch = coherent_tracking_state(chi0, deltaU, ZW, epsW, epsEta, Lambda, zeta)
    q = direct_branch_packet(
        branch["Rtr"], branch["Rtarget"], epsEta,
        Rtr_ref, Rtarget_ref, epsEta_ref,
        chi0_star, deltaU_star,
    )

    Pi = Pi_tr(xi, delta, R)
    Cmix = branch["Ctr"]
    regime = support_regime_labels(Pi, Cmix)
    zetareq = zeta_req_from_product(Pi, Cmix, branch["eps"])
    twin0 = coherent_tracking_state(chi0, deltaU, ZW, epsW, epsEta, Lambda, sp.Integer(1))

    subbanner("I. Exact actual-branch orbit packet")
    print("Finite quotient packet =")
    sp.pprint(sp.Matrix([q["qtr"], q["qnt"], q["qeta"]]))
    print("Weak-axisymmetric orbit lock on the actual coherent branch is exactly")
    print("  q_tr = q_nt = q_eta = 0,")
    print("or equivalently")
    print("  d ln R_tr = 0,  d ln R_target = 0,  d ln epsilon_eta = 0.")

    subbanner("II. Exact actual-branch support packet")
    print("epsilon =")
    sp.pprint(branch["eps"])
    print("C_mix = 8 Lambda (1-epsilon)/pi^2 =")
    sp.pprint(Cmix)
    print("Pi_tr =")
    sp.pprint(Pi)
    print("zeta_req =")
    sp.pprint(zetareq)

    subbanner("III. Exact regime split on the same physical branch")
    print("Mixed-only enough iff")
    print("  Pi_tr <= C_mix.")
    print("Lowest symmetric twin enough iff")
    print("  C_mix < Pi_tr <= 2 C_mix.")
    print("Non-twin asymmetry required iff")
    print("  Pi_tr > 2 C_mix.")
    print("Boundary vector =")
    sp.pprint(sp.Matrix([Pi - regime["mixed_only_boundary"], Pi - regime["lowest_twin_boundary"]]))

    subbanner("IV. Lowest-twin specialization")
    print("On the lowest symmetric twin branch zeta = 1, so")
    print("M_tr(zeta=1) =")
    sp.pprint(twin0["Mtr"])
    print("R_target * M_tr(zeta=1) =")
    sp.pprint(sp.simplify(twin0["Rtarget"] * twin0["Mtr"]))
    expect_zero(
        "lowest-twin product law",
        sp.simplify(twin0["Rtarget"] * twin0["Mtr"] - 16 * Lambda * (1 - twin0["eps"]) / sp.pi**2),
    )

    subbanner("V. Exact actual-branch theorem gate")
    print("The actual coherent local D/N branch now ends at one exact theorem gate:")
    print()
    print("Orbit packet:")
    print("  q_tr = q_nt = q_eta = 0.")
    print()
    print("Support packet:")
    print("  if Pi_tr <= C_mix, the mixed lane already suffices;")
    print("  if C_mix < Pi_tr <= 2 C_mix, the lowest symmetric twin lane suffices;")
    print("  if Pi_tr > 2 C_mix, non-twin asymmetry is required.")
    print()
    print("So the next remaining PDE-side job is no longer algebraic compression.")
    print("It is to determine the actual branch point Pi_tr and then decide which of the")
    print("three exact support regimes the physical operator lands in.")
