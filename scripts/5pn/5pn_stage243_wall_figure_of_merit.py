#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 243 — dimensionless wall figure of merit and explicit fail/succeed thresholds
for the first thin-wall parent confinement family.
"""


def main() -> None:
    banner("STAGE 243 — WALL FIGURE OF MERIT AND THRESHOLDS")

    V0, a, ell, KX, L, TX = sp.symbols("V0 a ell K_X L T_X", positive=True, real=True)
    J1, If, Hw = sp.symbols("J_1 I_f H_w", positive=True, real=True)
    Pe_req = sp.symbols("Pe_req", positive=True, real=True)
    Delta0, DeltaInf = sp.symbols("Delta_0 Delta_inf", positive=True, real=True)
    m, csw, rho_w = sp.symbols("m csw rho_w", positive=True, real=True)

    G_eq_tw = sp.simplify(4 * sp.pi * a**2 * V0**2 * J1 / (KX * ell))
    kappa = sp.simplify(KX * L**2 / TX)
    W_wall = sp.simplify(kappa * G_eq_tw)

    subbanner("I. Dimensionless wall figure of merit")
    print("W_wall =")
    sp.pprint(W_wall)
    expect_zero(
        "W_wall identity",
        sp.simplify(W_wall - 4 * sp.pi * a**2 * L**2 * J1 * V0**2 / (TX * ell)),
    )

    subbanner("II. Exact thin-wall fail/succeed thresholds")
    G_fail = sp.simplify(Pe_req / (kappa * DeltaInf))
    G_suff = sp.simplify(Pe_req / (kappa * Delta0))

    V0_fail_sq = sp.simplify(KX * ell * G_fail / (4 * sp.pi * a**2 * J1))
    V0_suff_sq = sp.simplify(KX * ell * G_suff / (4 * sp.pi * a**2 * J1))

    V0_fail_sq_reduced = sp.simplify(TX * ell * Pe_req / (4 * sp.pi * a**2 * L**2 * J1 * DeltaInf))
    V0_suff_sq_reduced = sp.simplify(TX * ell * Pe_req / (4 * sp.pi * a**2 * L**2 * J1 * Delta0))

    expect_zero("KX-cancelled fail threshold", sp.simplify(V0_fail_sq - V0_fail_sq_reduced))
    expect_zero("KX-cancelled suff threshold", sp.simplify(V0_suff_sq - V0_suff_sq_reduced))

    print("V0_fail^2 =")
    sp.pprint(V0_fail_sq_reduced)
    print()
    print("V0_suff^2 =")
    sp.pprint(V0_suff_sq_reduced)

    W_fail = sp.simplify(kappa * G_fail)
    W_suff = sp.simplify(kappa * G_suff)
    expect_zero("W_fail identity", sp.simplify(W_fail - Pe_req / DeltaInf))
    expect_zero("W_suff identity", sp.simplify(W_suff - Pe_req / Delta0))

    subbanner("III. Constant-compressibility wall layer")
    W_H = sp.simplify(W_wall.subs({J1: If / Hw}))
    print("W_H =")
    sp.pprint(W_H)
    expect_zero(
        "W_H identity",
        sp.simplify(W_H - 4 * sp.pi * a**2 * L**2 * If * V0**2 / (Hw * TX * ell)),
    )

    Hw_id = sp.simplify(m * csw**2 / rho_w)
    V0_fail_sq_H = sp.simplify(V0_fail_sq_reduced.subs({J1: If / Hw, Hw: Hw_id}))
    V0_suff_sq_H = sp.simplify(V0_suff_sq_reduced.subs({J1: If / Hw, Hw: Hw_id}))

    print()
    print("Constant-compressibility fail threshold =")
    sp.pprint(V0_fail_sq_H)
    print()
    print("Constant-compressibility suff threshold =")
    sp.pprint(V0_suff_sq_H)

    banner("STAGE 243 LEDGER")
    print("1. The first explicit thin-wall parent branch is controlled by")
    print("      W_wall = 4 pi a^2 L^2 J1 V0^2 / (T_X ell).")
    print("2. The exact branch test is")
    print("      W_wall <= Pe_req/Delta_inf   -> fail,")
    print("      W_wall >= Pe_req/Delta_0     -> succeed.")
    print("3. The physical wall-amplitude thresholds are")
    print("      V0_fail^2 = T_X ell Pe_req / [4 pi a^2 L^2 J1 Delta_inf],")
    print("      V0_suff^2 = T_X ell Pe_req / [4 pi a^2 L^2 J1 Delta_0].")
    print("4. In the constant-compressibility layer, J1 = I_f/H_w and the same theorem becomes a direct V0 test.")
    print("5. So the support/source gate is now a one-number comparison, not an abstract gain problem.")


if __name__ == "__main__":
    main()
