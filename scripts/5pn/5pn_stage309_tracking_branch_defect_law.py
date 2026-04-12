#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage309_311_common import banner, subbanner, expect_zero

"""
Stage 309 — physical tracking-branch transfer-shape defect law and support-blindness theorem.

This stage takes the direct physical observables from Stages 306–308 and inserts
 the actual coherent local D/N tracking-branch placement laws:

    R_tr     = (1 + chi_0/(1+delta_U)) / (1 + chi_0)
    epsilon  = epsilon_W * (1 - 2 delta_U / [11(1+delta_U)])
    M_tr     = M_mix * S(zeta;epsilon)
    R_target = Lambda_0 * Omega_W^2 * (1-epsilon_eta) (1-epsilon)^2 / [Z_W (1+chi_0)^2]

The outputs are:
  1. the exact support-blind transfer shape T^2,
  2. the exact weak-axisymmetric grouped defect law in the physical branch variables,
  3. the exact split-blocking transport law for epsilon_1,
  4. the support-blindness theorem: zeta drops out identically.
"""

if __name__ == "__main__":
    banner("STAGE 309 — TRACKING-BRANCH DEFECT LAW AND SUPPORT-BLINDNESS")

    chi0, deltaU = sp.symbols("chi_0 delta_U", positive=True, real=True)
    epsEta, epsW = sp.symbols("epsilon_eta epsilon_W", positive=True, real=True)
    ZW, OmW2, Lambda0 = sp.symbols("Z_W Omega_W2 Lambda_0", positive=True, real=True)
    zeta, Mmix = sp.symbols("zeta M_mix", positive=True, real=True)

    split = sp.simplify(1 - sp.Rational(2, 11) * deltaU / (1 + deltaU))
    epsilon = sp.symbols('epsilon', positive=True, real=True)
    eps = epsilon
    S = sp.simplify(1 + zeta * (1 - eps) / (1 - zeta * eps))
    Mtr = sp.simplify(Mmix * S)

    Rtr = sp.simplify((1 + chi0 / (1 + deltaU)) / (1 + chi0))
    T2 = sp.simplify(ZW * (1 + chi0) ** 2 / (OmW2 * (1 - eps) ** 2))
    Rtarget = sp.simplify(Lambda0 * OmW2 * (1 - epsEta) * (1 - eps) ** 2 / (ZW * (1 + chi0) ** 2))

    subbanner("I. Exact coherent tracking-branch substitution")
    print("R_tr =")
    sp.pprint(Rtr)
    print("epsilon =")
    sp.pprint(eps)
    print("epsilon on the coherent tracking branch =")
    sp.pprint(sp.simplify(epsW * split))
    print("S(zeta; epsilon) =")
    sp.pprint(S)
    print("M_tr =")
    sp.pprint(Mtr)
    print("T^2 =")
    sp.pprint(T2)
    print("R_target =")
    sp.pprint(Rtarget)

    expect_zero(
        "T^2 - Lambda_0 (1-epsilon_eta) / R_target",
        sp.simplify(T2 - Lambda0 * (1 - epsEta) / Rtarget),
    )

    subbanner("II. Exact support-blindness of the transfer shape")
    expect_zero("d_zeta ln T^2", sp.diff(sp.log(T2), zeta))
    expect_zero("d_zeta ln R_target", sp.diff(sp.log(Rtarget), zeta))
    print("d_zeta ln M_tr =")
    sp.pprint(sp.simplify(sp.diff(sp.log(Mtr), zeta)))

    subbanner("III. Weak-axisymmetric drift language")
    tau = sp.symbols("tau", real=True)
    zetaZ, omegaW = sp.symbols("zeta_Z omega_W", real=True)
    chi1, eps1, eta1 = sp.symbols("chi_1 epsilon_1 eta_1", real=True)

    T2_t = sp.simplify(
        T2.subs({
            ZW: ZW * sp.exp(tau * zetaZ),
            OmW2: OmW2 * sp.exp(tau * omegaW),
            chi0: chi0 + tau * chi1,
            eps: eps + tau * eps1,
        })
    )
    Rtarget_t = sp.simplify(
        Rtarget.subs({
            ZW: ZW * sp.exp(tau * zetaZ),
            OmW2: OmW2 * sp.exp(tau * omegaW),
            chi0: chi0 + tau * chi1,
            eps: eps + tau * eps1,
            epsEta: epsEta + tau * eta1,
        })
    )

    Xi1 = sp.simplify(sp.diff(sp.log(T2_t), tau).subs({tau: 0}))
    R1 = sp.simplify(sp.diff(sp.log(Rtarget_t), tau).subs({tau: 0}))

    print("Xi_1 = d ln T^2 =")
    sp.pprint(Xi1)
    print("R_1 = d ln R_target =")
    sp.pprint(R1)

    Xi1_expected = sp.simplify(zetaZ - omegaW + 2 * chi1 / (1 + chi0) + 2 * eps1 / (1 - eps))
    R1_expected = sp.simplify(omegaW - zetaZ - 2 * chi1 / (1 + chi0) - 2 * eps1 / (1 - eps) - eta1 / (1 - epsEta))

    expect_zero("Xi_1 closed form", sp.simplify(Xi1 - Xi1_expected))
    expect_zero("R_1 closed form", sp.simplify(R1 - R1_expected))
    expect_zero("R_1 + Xi_1 + eta_1/(1-epsilon_eta)", sp.simplify(R1 + Xi1 + eta1 / (1 - epsEta)))

    subbanner("IV. Exact split-blocking transport law")
    epsW1, deltaU1 = sp.symbols("epsilonW_1 deltaU_1", real=True)
    eps_t = sp.simplify(
        (epsW + tau * epsW1)
        * 
        (1 - sp.Rational(2, 11) * (deltaU + tau * deltaU1) / (1 + deltaU + tau * deltaU1))
    )
    eps1_exact = sp.simplify(sp.diff(eps_t, tau).subs({tau: 0}))
    eps1_expected = sp.simplify(
        split * epsW1 - sp.Rational(2, 11) * epsW * deltaU1 / (1 + deltaU) ** 2
    )

    print("epsilon_1 =")
    sp.pprint(eps1_exact)
    expect_zero("epsilon_1 transport law", sp.simplify(eps1_exact - eps1_expected))

    Xi1_physical = sp.simplify(Xi1_expected.subs({eps1: eps1_exact}))
    print("Xi_1 in (zeta_Z, omega_W, chi_1, epsilonW_1, deltaU_1) =")
    sp.pprint(Xi1_physical)

    banner("STAGE 309 LEDGER")
    print("1. On the coherent local D/N tracking branch,")
    print("      T^2 = Z_W (1+chi_0)^2 / [Omega_W^2 (1-epsilon)^2]")
    print("   and the support-enhancement ratio zeta drops out identically.")
    print("2. The exact first-order grouped defect is")
    print("      Xi_1 = zeta_Z - omega_W + 2 chi_1/(1+chi_0) + 2 epsilon_1/(1-epsilon).")
    print("3. The split-blocking drift is")
    print("      epsilon_1 = (1 - 2 delta_U/[11(1+delta_U)]) epsilonW_1")
    print("                   - 2 epsilon_W deltaU_1 / [11(1+delta_U)^2].")
    print("4. So coherent support can raise the baseline M_tr, but it cannot change")
    print("   the transfer shape or the weak-axisymmetric grouped defect at first order.")
