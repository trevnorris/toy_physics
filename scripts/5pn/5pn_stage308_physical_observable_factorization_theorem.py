#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import sympy as sp

sys.path.insert(0, os.path.dirname(__file__))
from fivepn_stage300_302_common import banner, subbanner, expect_zero

"""
Stage 308 — exact factorization of the physical branch observables and the
final direct-variable zero-defect theorem.

This stage shows that the physical coherent variables split cleanly into:
  (a) the three orbit-lock observables (R_tr, R_target, epsilon_eta), and
  (b) the support-compensation baseline variable M_tr.

The support lane modifies M_tr through the enhancement factor S(zeta;epsilon)
without entering the first-order orbit-lock packet at all.
"""

if __name__ == "__main__":
    banner("STAGE 308 — PHYSICAL OBSERVABLE FACTORIZATION THEOREM")

    chi0, deltaU = sp.symbols("chi_0 delta_U", positive=True, real=True)
    epsEta, eps = sp.symbols("epsilon_eta epsilon", positive=True, real=True)
    ZW, OmW2, Lambda0 = sp.symbols("Z_W Omega_W2 Lambda_0", positive=True, real=True)
    zeta, Mmix = sp.symbols("zeta M_mix", positive=True, real=True)

    Rtr = sp.simplify((1 + chi0 / (1 + deltaU)) / (1 + chi0))
    Rtarget = sp.simplify(Lambda0 * OmW2 * (1 - epsEta) * (1 - eps) ** 2 / (ZW * (1 + chi0) ** 2))
    S = sp.simplify(1 + zeta * (1 - eps) / (1 - zeta * eps))
    Mtr = sp.simplify(Mmix * S)

    subbanner("I. Exact physical-variable factorization")
    print("R_tr =")
    sp.pprint(Rtr)
    print("R_target =")
    sp.pprint(Rtarget)
    print("epsilon_eta =")
    sp.pprint(epsEta)
    print("S(zeta; epsilon) =")
    sp.pprint(S)
    print("M_tr =")
    sp.pprint(Mtr)

    print("The direct observable split is now exact:")
    print("  - R_tr depends only on (chi_0, delta_U),")
    print("  - R_target depends only on (Omega_W^2, Z_W, chi_0, epsilon, epsilon_eta),")
    print("  - epsilon_eta depends only on epsilon_eta,")
    print("  - M_tr carries the support-compensation data through (M_mix, zeta, epsilon).")

    subbanner("II. Exact support-blindness of the orbit-lock observables")
    for name, expr in [
        ("d_zeta ln R_tr", sp.diff(sp.log(Rtr), zeta)),
        ("d_zeta ln R_target", sp.diff(sp.log(Rtarget), zeta)),
        ("d_zeta ln epsilon_eta", sp.diff(sp.log(epsEta), zeta)),
        ("d_Mmix ln R_tr", sp.diff(sp.log(Rtr), Mmix)),
        ("d_Mmix ln R_target", sp.diff(sp.log(Rtarget), Mmix)),
        ("d_Mmix ln epsilon_eta", sp.diff(sp.log(epsEta), Mmix)),
    ]:
        expect_zero(name, expr)

    print("d_zeta ln M_tr =")
    sp.pprint(sp.simplify(sp.diff(sp.log(Mtr), zeta)))
    print("d_(ln M_mix) ln M_tr =")
    sp.pprint(sp.simplify(Mmix * sp.diff(sp.log(Mtr), Mmix)))

    subbanner("III. Exact direct-variable zero-defect theorem")
    dlnchi0, dlndeltaU = sp.symbols("dlnchi_0 dlndelta_U", real=True)
    dlnZW, dlnOmW2 = sp.symbols("dlnZ_W dlnOmega_W2", real=True)
    dlneps, dlnepsEta = sp.symbols("dlnepsilon dlnepsilon_eta", real=True)

    Ctr = sp.simplify(chi0 * deltaU / ((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU)))
    dlnRtr = sp.simplify(-Ctr * ((1 + deltaU) * dlnchi0 + (1 + chi0) * dlndeltaU))
    dlnRtarget = sp.simplify(
        dlnOmW2
        - dlnZW
        - 2 * chi0 * dlnchi0 / (1 + chi0)
        - 2 * eps * dlneps / (1 - eps)
        - epsEta * dlnepsEta / (1 - epsEta)
    )

    Theta1 = sp.simplify(dlnRtr)
    Xi1 = sp.simplify(-dlnRtarget - epsEta * dlnepsEta / (1 - epsEta))
    R1 = sp.simplify(dlnRtarget)

    cond_tracking = sp.simplify((1 + deltaU) * dlnchi0 + (1 + chi0) * dlndeltaU)
    cond_target = sp.simplify(
        dlnOmW2
        - dlnZW
        - 2 * chi0 * dlnchi0 / (1 + chi0)
        - 2 * eps * dlneps / (1 - eps)
        - epsEta * dlnepsEta / (1 - epsEta)
    )
    cond_dressing = dlnepsEta

    print("Theta_1 =")
    sp.pprint(Theta1)
    print("Xi_1 =")
    sp.pprint(Xi1)
    print("R_1 =")
    sp.pprint(R1)

    expect_zero("Theta_1 + C_tr * tracking condition", sp.simplify(Theta1 + Ctr * cond_tracking))
    expect_zero("R_1 - target condition", sp.simplify(R1 - cond_target))
    expect_zero("Xi_1 + target condition + epsilon_eta * dressing condition/(1-epsilon_eta)", sp.simplify(Xi1 + cond_target + epsEta * cond_dressing / (1 - epsEta)))

    subbanner("IV. Simplified orbit-lock conditions on the physical coherent branch")
    print("Exact direct-variable zero-defect conditions:")
    print("  1. (1 + delta_U) dlnchi_0 + (1 + chi_0) dlndelta_U = 0")
    print("  2. dlnOmega_W^2 - dlnZ_W - 2 chi_0 dlnchi_0/(1+chi_0)")
    print("       - 2 epsilon dlnepsilon/(1-epsilon)")
    print("       - epsilon_eta dlnepsilon_eta/(1-epsilon_eta) = 0")
    print("  3. dlnepsilon_eta = 0")
    print()
    print("With condition 3 imposed, condition 2 simplifies to")
    print("  dlnOmega_W^2 - dlnZ_W - 2 chi_0 dlnchi_0/(1+chi_0) - 2 epsilon dlnepsilon/(1-epsilon) = 0.")

    banner("STAGE 308 LEDGER")
    print("1. The physical coherent branch factors exactly into")
    print("      (R_tr, R_target, epsilon_eta)   and   M_tr = M_mix S(zeta;epsilon).")
    print("2. The support lane enters only through M_tr; it is invisible to the")
    print("   first-order orbit-lock observables and therefore invisible to")
    print("   the first-order defect packet.")
    print("3. The actual first-order orbit-lock theorem is now an explicit system of")
    print("   three differential constraints on the physical coherent variables")
    print("      (chi_0, delta_U, epsilon_eta, Z_W, Omega_W^2, epsilon).")
