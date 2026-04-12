#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage300_302_common import banner, subbanner, expect_zero

"""
Stage 306 — exact coherent-branch observable map in physical variables.

This stage identifies the actual coherent branch observables used by the
Stage-303/304 projector directly in the physical coherent-kernel variables:
    chi_0, delta_U, epsilon_eta, Z_W, Omega_W^2, epsilon.

The main point is that the direct branch observables are not abstract anymore.
They are explicit branch functions:
    R_tr     = (1 + chi_0/(1+delta_U)) / (1 + chi_0),
    R_target = Lambda (1-epsilon_eta) (1-epsilon)^2 / [ Z_W (1+chi_0)^2 ],
    epsilon_eta = epsilon_eta.

So the Stage-303 finite quotient packet can now be evaluated directly on the
physical coherent branch.
"""

if __name__ == "__main__":
    banner("STAGE 306 — COHERENT BRANCH OBSERVABLE MAP")

    chi0, deltaU = sp.symbols("chi_0 delta_U", positive=True, real=True)
    epsEta, eps = sp.symbols("epsilon_eta epsilon", positive=True, real=True)
    ZW, OmW2, Lambda0 = sp.symbols("Z_W Omega_W2 Lambda_0", positive=True, real=True)

    Rtr = sp.simplify((1 + chi0 / (1 + deltaU)) / (1 + chi0))
    T2 = sp.simplify(ZW * (1 + chi0) ** 2 / (OmW2 * (1 - eps) ** 2))
    Rtarget = sp.simplify(Lambda0 * OmW2 * (1 - epsEta) * (1 - eps) ** 2 / (ZW * (1 + chi0) ** 2))

    subbanner("I. Direct coherent-branch observables")
    print("R_tr =")
    sp.pprint(Rtr)
    print("T^2 =")
    sp.pprint(T2)
    print("R_target =")
    sp.pprint(Rtarget)
    print("epsilon_eta =")
    sp.pprint(epsEta)

    expect_zero(
        "R_target * T^2 - Lambda_0 * (1 - epsilon_eta)",
        sp.simplify(Rtarget * T2 - Lambda0 * (1 - epsEta)),
    )

    subbanner("II. Finite quotient packet in the actual branch observables")
    chi0s, deltaUs, epsEtas = sp.symbols(
        "chi0_star deltaU_star epsilon_eta_star", positive=True, real=True
    )
    Rtr_ref, Rtarget_ref, epsEta_ref = sp.symbols(
        "R_tr_ref R_target_ref epsilon_eta_ref", positive=True, real=True
    )

    Cstar = sp.simplify(
        (1 + chi0s) * (1 + deltaUs) * (1 + chi0s + deltaUs) / (chi0s * deltaUs)
    )
    Bstar = sp.simplify(2 * (1 + chi0s + deltaUs) / deltaUs)

    qtr = sp.simplify(-Cstar * sp.log(Rtr / Rtr_ref))
    qnt = sp.simplify(
        Bstar * sp.log(Rtr / Rtr_ref)
        + sp.log((1 - epsEta) / (1 - epsEta_ref))
        - sp.log(Rtarget / Rtarget_ref)
    )
    qeta = sp.simplify(sp.log(epsEta / epsEta_ref))

    print("q_tr =")
    sp.pprint(qtr)
    print("q_nt =")
    sp.pprint(qnt)
    print("q_eta =")
    sp.pprint(qeta)

    subbanner("III. Exact physical meaning")
    print("The Stage-303 direct quotient packet now has a concrete coherent-branch realization:")
    print("  R_tr     = (1 + chi_0/(1+delta_U)) / (1+chi_0)")
    print("  R_target = Lambda_0 Omega_W^2 (1-epsilon_eta)(1-epsilon)^2 / [Z_W (1+chi_0)^2]")
    print("  epsilon_eta = epsilon_eta")
    print()
    print("So the actual physical branch already supplies the finite observable packet")
    print("  (R_tr, R_target, epsilon_eta),")
    print("from which the Stage-303 quotient coordinates follow exactly.")
