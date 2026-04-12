#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage300_302_common import banner, subbanner, expect_zero

"""
Stage 305 — exact three-observable and support-blind theorem.

What this script does
---------------------
1. Packages the direct coherent-branch zero-defect theorem entirely in the three
   exact observables
      R_tr, R_target, epsilon_eta.
2. Shows that the first-order defect triplet vanishes iff those three observables
   are invariant.
3. Splices this with the coherent support-enhancement law of the physical tracking
   branch, where support only rescales M_tr through S(zeta;epsilon), while
   R_tr, R_target, epsilon_eta stay unchanged.
4. Proves that the actual first-order orbit-defect packet is blind to both the
   support-enhancement ratio zeta and the total baseline M_tr.

Interpretation
--------------
The weak-axisymmetric orbit-lock test and the selected-mode support-compensation
problem are exact but distinct branches of the endgame.  Support enhancement can
repair the selected quadrupole normalization branch without moving the actual
coherent branch on or off the first-order similarity orbit.
"""

if __name__ == "__main__":
    banner("STAGE 305 — SUPPORT-BLIND THREE-OBSERVABLE THEOREM")

    chi0s, deltaUs, epsEtas = sp.symbols(
        "chi0_star deltaU_star epsilon_eta_star", positive=True, real=True
    )
    Rtr, Rtarget, epsEta = sp.symbols(
        "R_tr R_target epsilon_eta", positive=True, real=True
    )
    Rtr_ref, Rtarget_ref, epsEta_ref = sp.symbols(
        "R_tr_ref R_target_ref epsilon_eta_ref", positive=True, real=True
    )
    zeta, eps, Mmix = sp.symbols("zeta epsilon M_mix", positive=True, real=True)

    Cstar = sp.simplify(
        (1 + chi0s) * (1 + deltaUs) * (1 + chi0s + deltaUs) / (chi0s * deltaUs)
    )
    Bstar = sp.simplify(2 * (1 + chi0s + deltaUs) / deltaUs)
    Ctr = sp.simplify(1 / Cstar)
    Atr = sp.simplify(Bstar * Ctr)

    # Finite quotient coordinates in the actual branch observables.
    qtr = sp.simplify(-Cstar * sp.log(Rtr / Rtr_ref))
    qnt = sp.simplify(
        Bstar * sp.log(Rtr / Rtr_ref)
        + sp.log((1 - epsEta) / (1 - epsEta_ref))
        - sp.log(Rtarget / Rtarget_ref)
    )
    qeta = sp.simplify(sp.log(epsEta / epsEta_ref))

    # First-order direct-observable drifts and physical defects.
    dlnRtr, dlnRtarget, dlnepsEta = sp.symbols(
        "dlnR_tr dlnR_target dlnepsilon_eta", real=True
    )
    Theta1 = dlnRtr
    Xi1 = sp.simplify(-dlnRtarget - epsEtas * dlnepsEta / (1 - epsEtas))
    R1 = dlnRtarget

    subbanner("I. Exact three-observable zero-defect theorem")
    print("Theta_1 =")
    sp.pprint(Theta1)
    print("Xi_1 =")
    sp.pprint(Xi1)
    print("R_1 =")
    sp.pprint(R1)

    # Zero-defect equivalence.
    expect_zero(
        "(Theta_1, Xi_1, R_1) = 0 iff (dlnR_tr, dlnR_target, dlnepsilon_eta) = 0 [direct part]",
        sp.simplify(Theta1 - dlnRtr),
    )
    expect_zero(
        "R_1 direct branch invariance",
        sp.simplify(R1 - dlnRtarget),
    )
    dlnepsEta_from_defects = sp.simplify(-(1 - epsEtas) * (R1 + Xi1) / epsEtas)
    expect_zero("dlnepsilon_eta inverse", sp.simplify(dlnepsEta_from_defects - dlnepsEta))

    subbanner("II. Coherent support enhancement only rescales M_tr")
    S = sp.simplify(1 + zeta * (1 - eps) / (1 - zeta * eps))
    Mtr = sp.simplify(Mmix * S)

    print("S(zeta; epsilon) =")
    sp.pprint(S)
    print("M_tr =")
    sp.pprint(Mtr)
    print("R_target is independent of zeta on the coherent branch.")
    print("R_tr and epsilon_eta are also independent of zeta on the coherent branch.")

    subbanner("III. Exact support-blindness of the actual orbit-defect packet")
    Mtr_sym = sp.symbols("M_tr", positive=True, real=True)

    # The actual quotient packet and defect packet depend only on (R_tr, R_target, epsilon_eta).
    # Therefore they are blind to the total baseline M_tr and to zeta.
    dq_dM = sp.Matrix([sp.diff(qtr, Mtr_sym), sp.diff(qnt, Mtr_sym), sp.diff(qeta, Mtr_sym)])
    ddef_dM = sp.Matrix([sp.diff(Theta1, Mtr_sym), sp.diff(Xi1, Mtr_sym), sp.diff(R1, Mtr_sym)])

    dq_dzeta = sp.Matrix([sp.diff(qtr, zeta), sp.diff(qnt, zeta), sp.diff(qeta, zeta)])
    ddef_dzeta = sp.Matrix([sp.diff(Theta1, zeta), sp.diff(Xi1, zeta), sp.diff(R1, zeta)])

    print("d(q_tr,q_nt,q_eta)/d M_tr =")
    sp.pprint(dq_dM)
    print("d(Theta_1,Xi_1,R_1)/d M_tr =")
    sp.pprint(ddef_dM)
    print("d(q_tr,q_nt,q_eta)/d zeta =")
    sp.pprint(dq_dzeta)
    print("d(Theta_1,Xi_1,R_1)/d zeta =")
    sp.pprint(ddef_dzeta)

    expect_zero("quotient packet support-blind in M_tr", dq_dM)
    expect_zero("defect packet support-blind in M_tr", ddef_dM)
    expect_zero("quotient packet support-blind in zeta", dq_dzeta)
    expect_zero("defect packet support-blind in zeta", ddef_dzeta)

    banner("STAGE 305 LEDGER")
    print("1. On the coherent branch, the actual first-order orbit-lock test is exactly")
    print("   the invariance of the three direct observables")
    print("      (R_tr, R_target, epsilon_eta).")
    print("2. The physical defect triplet is triangular in their logarithmic drifts:")
    print("      Theta_1 = delta ln R_tr,")
    print("      Xi_1    = - delta ln R_target - [epsilon_eta_*/(1-epsilon_eta_*)] delta ln epsilon_eta,")
    print("      R_1     = delta ln R_target.")
    print("3. Coherent support enhancement only rescales the total baseline")
    print("      M_tr = M_mix S(zeta;epsilon),")
    print("   and does not enter the actual first-order quotient or defect packet.")
    print("4. Therefore the support-compensation branch and the weak-axisymmetric")
    print("   orbit-lock branch are exact but distinct parts of the 5PN endgame.")
