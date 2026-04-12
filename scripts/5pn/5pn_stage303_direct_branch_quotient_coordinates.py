#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage300_302_common import banner, subbanner, expect_zero

"""
Stage 303 — exact finite quotient coordinates in direct branch-observable language.

What this script does
---------------------
1. Replaces the older reduced-ratio monomials by the direct coherent-branch
   observables
      R_tr, R_target, epsilon_eta.
2. Uses the exact branch identities
      T_* = R_tr^(-C_*),
      N_* = Lambda_0 (1-epsilon_eta) R_tr^(B_*) / R_target
   to write the finite quotient coordinates directly in terms of those three
   observables and their reference values.
3. Solves the exact inverse map.

Interpretation
--------------
On the coherent local tracking branch, the actual moving-throat weak-axisymmetric
quotient packet no longer needs the five reduced continuum ratios explicitly.
It can be computed directly from the three exact branch observables
    (R_tr, R_target, epsilon_eta).
"""

if __name__ == "__main__":
    banner("STAGE 303 — DIRECT BRANCH QUOTIENT COORDINATES")

    chi0s, deltaUs, epsEtas = sp.symbols(
        "chi0_star deltaU_star epsilon_eta_star", positive=True, real=True
    )
    Rtr, Rtarget, epsEta = sp.symbols(
        "R_tr R_target epsilon_eta", positive=True, real=True
    )
    Rtr_ref, Rtarget_ref, epsEta_ref = sp.symbols(
        "R_tr_ref R_target_ref epsilon_eta_ref", positive=True, real=True
    )
    Lambda0 = sp.symbols("Lambda_0", positive=True, real=True)

    Cstar = sp.simplify(
        (1 + chi0s) * (1 + deltaUs) * (1 + chi0s + deltaUs) / (chi0s * deltaUs)
    )
    Bstar = sp.simplify(2 * (1 + chi0s + deltaUs) / deltaUs)

    Tstar = sp.simplify(Rtr ** (-Cstar))
    Tstar_ref = sp.simplify(Rtr_ref ** (-Cstar))
    Nstar = sp.simplify(Lambda0 * (1 - epsEta) * Rtr ** Bstar / Rtarget)
    Nstar_ref = sp.simplify(Lambda0 * (1 - epsEta_ref) * Rtr_ref ** Bstar / Rtarget_ref)

    qtr = sp.simplify(sp.log(Tstar / Tstar_ref))
    qnt = sp.simplify(sp.log(Nstar / Nstar_ref))
    qeta = sp.simplify(sp.log(epsEta / epsEta_ref))

    subbanner("I. Exact finite quotient coordinates")
    print("C_* =")
    sp.pprint(Cstar)
    print("B_* =")
    sp.pprint(Bstar)
    print("q_tr =")
    sp.pprint(qtr)
    print("q_nt =")
    sp.pprint(qnt)
    print("q_eta =")
    sp.pprint(qeta)

    qtr_expected = sp.simplify(-Cstar * sp.log(Rtr / Rtr_ref))
    qnt_expected = sp.simplify(
        Bstar * sp.log(Rtr / Rtr_ref)
        + sp.log((1 - epsEta) / (1 - epsEta_ref))
        - sp.log(Rtarget / Rtarget_ref)
    )
    qeta_expected = sp.simplify(sp.log(epsEta / epsEta_ref))

    expect_zero("q_tr direct-observable form", sp.simplify(qtr - qtr_expected))
    expect_zero("q_nt direct-observable form", sp.simplify(qnt - qnt_expected))
    expect_zero("q_eta direct-observable form", sp.simplify(qeta - qeta_expected))

    subbanner("II. Exact inverse map")
    Rtr_from_q = sp.simplify(Rtr_ref * sp.exp(-qtr / Cstar))
    epsEta_from_q = sp.simplify(epsEta_ref * sp.exp(qeta))
    Rtarget_from_q = sp.simplify(
        Rtarget_ref
        * sp.exp(-qnt - Bstar * qtr / Cstar)
        * (1 - epsEta_from_q)
        / (1 - epsEta_ref)
    )

    print("R_tr =")
    sp.pprint(Rtr_from_q)
    print("epsilon_eta =")
    sp.pprint(epsEta_from_q)
    print("R_target =")
    sp.pprint(Rtarget_from_q)

    expect_zero("R_tr inverse", sp.simplify(Rtr_from_q - Rtr))
    expect_zero("epsilon_eta inverse", sp.simplify(epsEta_from_q - epsEta))
    expect_zero("R_target inverse", sp.simplify(Rtarget_from_q - Rtarget))

    subbanner("III. Exact three-observable packet")
    print("The finite quotient packet on the coherent branch depends only on")
    print("  R_tr, R_target, epsilon_eta")
    print("and not on the remaining reduced continuum ratios explicitly.")

    banner("STAGE 303 LEDGER")
    print("1. The direct branch-invariant coordinates T_* and N_* are exact functions of")
    print("   the three actual branch observables R_tr, R_target, epsilon_eta.")
    print("2. Therefore the finite quotient packet (q_tr,q_nt,q_eta) is already available")
    print("   directly from those three observables without reopening the full five-ratio")
    print("   continuum state space.")
    print("3. The inverse map is exact, so the actual coherent branch may be charted")
    print("   equivalently by either")
    print("      (R_tr, R_target, epsilon_eta) or (q_tr, q_nt, q_eta).")
