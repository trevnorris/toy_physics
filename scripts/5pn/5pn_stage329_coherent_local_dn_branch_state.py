#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage329_332_common import *
from fivepn_stage326_328_common import microscopic_coherent_placement_state, lambda0_constant

"""
Stage 329 — concrete coherent local D/N branch state from the support-compensation family.

What this stage does
--------------------
1. Takes the first explicit coherent local D/N source family of Stages 30–33,
   where the mixed lane W and the support lane phi_n are sourced by the same local
   throat density.
2. Extracts the actual coherent branch state
      (chi0, deltaU, Z_W, epsilon_W, epsilon_eta, Lambda, zeta_n^(phys))
   directly from the microscopic kernel quotients.
3. Specializes to the first uniform-source and same-operator twin families.
4. Shows that the lowest symmetric twin lane has zeta_0 = 1 and S_0 = 2 exactly.
"""

if __name__ == "__main__":
    banner("STAGE 329 — CONCRETE COHERENT LOCAL D/N BRANCH STATE")

    # Microscopic kernel variables
    lam_star = sp.symbols("lambda_star", positive=True, real=True)
    I_0, I_n = sp.symbols("I_0 I_n", positive=True, real=True)
    cetaU, gamma = sp.symbols("c_etaU gamma", positive=True, real=True)
    KU, Ketaeff, KWeff, Kphi_n_eff, muW = sp.symbols(
        "K_U K_etaeff K_W_eff K_phi_n_eff mu_W", positive=True, real=True
    )
    TU, L, sigma = sp.symbols("T_U L sigma", positive=True, real=True)
    G, cs, a, c = sp.symbols("G c_s a c", positive=True, real=True)
    n = sp.symbols("n", integer=True, nonnegative=True)
    x = sp.symbols("x", positive=True, real=True)

    lambdaW = sp.simplify(lam_star * I_0)
    lambdaPhi = sp.simplify(lam_star * I_n)
    Lambda0 = lambda0_constant(G, cs, a, c)

    state = microscopic_coherent_placement_state(
        lambdaW, cetaU, gamma, KU, Ketaeff, KWeff, muW, TU, L, sigma, Lambda0, lambdaPhi, Kphi_n_eff
    )
    zeta_phys = zeta_phys_general(KWeff, Kphi_n_eff, I_n, I_0)
    ratio_uniform = I_uniform_ratio(n)
    zeta_uniform = sp.simplify(zeta_phys.subs({I_n / I_0: ratio_uniform}))
    zeta_twin = zeta_twin_same_operator(n, x)
    twin_state = coherent_tracking_state(state["chi0"], state["deltaU"], state["ZW"], state["epsW"], state["epsEta"], state["Lambda"], state["zeta"])
    S_twin = S_of_zeta(zeta_twin, twin_state["eps"])

    subbanner("I. Exact microscopic extraction of the actual coherent branch state")
    print("lambda_W =")
    sp.pprint(lambdaW)
    print("lambda_phi,n =")
    sp.pprint(lambdaPhi)
    print("chi0 =")
    sp.pprint(state["chi0"])
    print("deltaU =")
    sp.pprint(state["deltaU"])
    print("Z_W =")
    sp.pprint(state["ZW"])
    print("epsilon_W =")
    sp.pprint(state["epsW"])
    print("epsilon_eta =")
    sp.pprint(state["epsEta"])
    print("Lambda =")
    sp.pprint(state["Lambda"])
    print("zeta =")
    sp.pprint(state["zeta"])

    subbanner("II. Exact physical coherent support ratio")
    print("zeta_n^(phys) =")
    sp.pprint(zeta_phys)
    expect_zero("zeta extracted from microscopic state - zeta_n^(phys)", sp.simplify(state["zeta"] - zeta_phys))
    print("Uniform-source overlap ratio I_n/I_0 =")
    sp.pprint(ratio_uniform)
    print("zeta_n^(uniform local source) =")
    sp.pprint(zeta_uniform)
    print("zeta_n^(twin same-operator family) =")
    sp.pprint(zeta_twin)

    subbanner("III. Exact branch observables on the same physical family")
    print("epsilon =")
    sp.pprint(twin_state["eps"])
    print("R_tr =")
    sp.pprint(twin_state["Rtr"])
    print("M_mix =")
    sp.pprint(twin_state["Mmix"])
    print("R_target =")
    sp.pprint(twin_state["Rtarget"])

    subbanner("IV. Exact lowest symmetric twin lane")
    zeta0 = sp.simplify(zeta_twin.subs({n: 0}))
    S0 = sp.simplify(S_twin.subs({n: 0}))
    print("zeta_0^(twin) =")
    sp.pprint(zeta0)
    print("S_0 =")
    sp.pprint(S0)
    expect_zero("zeta_0^(twin) - 1", zeta0 - 1)
    expect_zero("S_0 - 2", S0 - 2)

    subbanner("V. Actual packet on the lowest twin branch")
    twin0 = coherent_tracking_state(state["chi0"], state["deltaU"], state["ZW"], state["epsW"], state["epsEta"], state["Lambda"], sp.Integer(1))
    print("M_mix =")
    sp.pprint(twin0["Mmix"])
    print("M_supp(zeta=1) =")
    sp.pprint(twin0["Msupp"])
    print("M_tr(zeta=1) =")
    sp.pprint(twin0["Mtr"])
    print("R_target * M_tr(zeta=1) =")
    sp.pprint(sp.simplify(twin0["Rtarget"] * twin0["Mtr"]))
    expect_zero(
        "lowest-twin product law",
        sp.simplify(twin0["Rtarget"] * twin0["Mtr"] - 16 * state["Lambda"] * (1 - twin0["eps"]) / sp.pi**2),
    )

    subbanner("VI. Interpretation")
    print("The actual coherent local D/N branch is now extracted directly from the")
    print("microscopic kernel family.  The mixed lane and the support lane are sourced")
    print("by the same local throat density, so the physical support ratio zeta is fixed")
    print("once the D/N support harmonic is chosen.")
    print()
    print("On the same-operator twin family, the lowest symmetric lane has")
    print("  zeta_0 = 1  and  S_0 = 2")
    print("exactly, so it is the universal doubling branch of the concrete coherent")
    print("support sector.")
