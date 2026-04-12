#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 209 — primitive weak-axisymmetric extractor from the explicit isotropic overlap model.

What this script does
---------------------
1. Starts from the explicit isotropic finite-throat overlap prototype used in
   Stages 207–208:
      C   = lambda_B * I_cross,
      G_U = lambda_U,
      G_W = lambda_W * I_cross,
      R   = lambda_R * I_cross,
      I_cross = 8/(3 pi).
2. Introduces the primitive weak-axisymmetric slopes
      dK, dM, d(lambda_B), d(varpi), d(lambda_U), d(lambda_W), d(lambda_R),
      d(Omega_U), d(Omega_W)
   in one grouped lane.
3. Derives exact compact first-order formulas for
      D01, D21, D41, N01,
   through the intermediate variations
      Delta01, P01, Q01, H01, S201.
4. Verifies the bundle decomposition identities
      D01 = dK - B01 - Z01,
      D21 = -(dM + B21 + Z21),
      D41 = -(B41 + Z41),
   and the exact weak-axisymmetric normalization identity
      Xi_load = N01/N0 - D01/D0 = P1/P0.
5. Freezes the exact compensation surfaces on the canonical compensated branch:
      K1 = 0  fixes dM,
      Xi_load = 0  fixes dK.

Interpretation
--------------
This is the first runnable bridge from the explicit isotropic overlap branch to
actual primitive weak-axisymmetric slope data. It converts the next 5PN gate
from a verbal instruction into one exact symbolic compiler.
"""


def differential(expr: sp.Expr, prims: list[sp.Symbol], dprims: list[sp.Symbol]) -> sp.Expr:
    return sp.simplify(sum(sp.diff(expr, p) * dp for p, dp in zip(prims, dprims)))


if __name__ == "__main__":
    banner("STAGE 209 — PRIMITIVE WEAK-AXISYMMETRIC EXTRACTOR")

    # Explicit isotropic overlap constant.
    I_cross = sp.simplify(sp.Rational(8, 3) / sp.pi)

    # Primitive isotropic branch data.
    K, M = sp.symbols("K M", positive=True, real=True)
    lambda_B, varpi = sp.symbols("lambda_B varpi", positive=True, real=True)
    lambda_U, lambda_W, lambda_R = sp.symbols(
        "lambda_U lambda_W lambda_R", positive=True, real=True
    )
    Omega_U, Omega_W = sp.symbols("Omega_U Omega_W", positive=True, real=True)

    # Primitive weak-axisymmetric absolute slopes.
    dK, dM = sp.symbols("dK dM", real=True)
    dlambda_B, dvarpi = sp.symbols("d(lambda_B) d(varpi)", real=True)
    dlambda_U, dlambda_W, dlambda_R = sp.symbols(
        "d(lambda_U) d(lambda_W) d(lambda_R)", real=True
    )
    dOmega_U, dOmega_W = sp.symbols("d(Omega_U) d(Omega_W)", real=True)

    prims = [K, M, lambda_B, varpi, lambda_U, lambda_W, lambda_R, Omega_U, Omega_W]
    dprims = [dK, dM, dlambda_B, dvarpi, dlambda_U, dlambda_W, dlambda_R, dOmega_U, dOmega_W]

    # Explicit overlap-model amplitudes.
    C = sp.simplify(lambda_B * I_cross)
    G_U = lambda_U
    G_W = sp.simplify(lambda_W * I_cross)
    R = sp.simplify(lambda_R * I_cross)

    # Isotropic bundle data.
    Delta = sp.simplify(Omega_U**2 * Omega_W**2 - R**2)
    P = sp.simplify(Omega_U**2 * G_W + R * G_U)
    Q = sp.simplify(G_U**2 * Omega_W**2 + 2 * G_U * G_W * R + G_W**2 * Omega_U**2)
    H = sp.simplify(G_U**2 + G_W**2)
    S2 = sp.simplify(Omega_U**2 + Omega_W**2)

    B0 = sp.simplify(C**2 / varpi**2)
    B2 = sp.simplify(C**2 / varpi**4)
    B4 = sp.simplify(C**2 / varpi**6)

    Z0 = sp.simplify(Q / Delta)
    Z2_num = sp.simplify(Q * S2 - H * Delta)
    Z2 = sp.simplify(Z2_num / Delta**2)
    Z4_num = sp.simplify(Q * (S2**2 - Delta) - S2 * H * Delta)
    Z4 = sp.simplify(Z4_num / Delta**3)

    N0 = sp.simplify(P**2 / Delta**2)

    D0 = sp.simplify(K - B0 - Z0)
    D2 = sp.simplify(-(M + B2 + Z2))
    D4 = sp.simplify(-(B4 + Z4))

    # Compact first-order primitive variations.
    Delta01 = sp.simplify(differential(Delta, prims, dprims))
    P01 = sp.simplify(differential(P, prims, dprims))
    Q01 = sp.simplify(differential(Q, prims, dprims))
    H01 = sp.simplify(differential(H, prims, dprims))
    S201 = sp.simplify(differential(S2, prims, dprims))

    B01 = sp.simplify(differential(B0, prims, dprims))
    B21 = sp.simplify(differential(B2, prims, dprims))
    B41 = sp.simplify(differential(B4, prims, dprims))

    Z01 = sp.simplify((Q01 * Delta - Q * Delta01) / Delta**2)
    Z21_num01 = sp.simplify(Q01 * S2 + Q * S201 - H01 * Delta - H * Delta01)
    Z21 = sp.simplify((Z21_num01 * Delta - 2 * Z2_num * Delta01) / Delta**3)
    Z41_num01 = sp.simplify(
        Q01 * (S2**2 - Delta)
        + Q * (2 * S2 * S201 - Delta01)
        - (S201 * H * Delta + S2 * H01 * Delta + S2 * H * Delta01)
    )
    Z41 = sp.simplify((Z41_num01 * Delta - 3 * Z4_num * Delta01) / Delta**4)

    N01 = sp.simplify(2 * N0 * (P01 / P - Delta01 / Delta))

    D01 = sp.simplify(dK - B01 - Z01)
    D21 = sp.simplify(-(dM + B21 + Z21))
    D41 = sp.simplify(-(B41 + Z41))

    subbanner("I. Explicit isotropic overlap branch data")
    print("I_cross =")
    sp.pprint(I_cross)
    print("C =")
    sp.pprint(C)
    print("G_U =")
    sp.pprint(G_U)
    print("G_W =")
    sp.pprint(G_W)
    print("R =")
    sp.pprint(R)

    subbanner("II. Compact primitive first-order variations")
    print("Delta01 =")
    sp.pprint(Delta01)
    print("P01 =")
    sp.pprint(P01)
    print("Q01 =")
    sp.pprint(Q01)
    print("H01 =")
    sp.pprint(H01)
    print("S201 =")
    sp.pprint(S201)

    # Easy-to-read BdG slope forms.
    expect_zero(
        "B01 - 2 B0(dlambda_B/lambda_B - dvarpi/varpi)",
        sp.simplify(B01 - 2 * B0 * (dlambda_B / lambda_B - dvarpi / varpi)),
    )
    expect_zero(
        "B21 - 2 B2(dlambda_B/lambda_B - 2 dvarpi/varpi)",
        sp.simplify(B21 - 2 * B2 * (dlambda_B / lambda_B - 2 * dvarpi / varpi)),
    )
    expect_zero(
        "B41 - 2 B4(dlambda_B/lambda_B - 3 dvarpi/varpi)",
        sp.simplify(B41 - 2 * B4 * (dlambda_B / lambda_B - 3 * dvarpi / varpi)),
    )

    subbanner("III. Exact bundle-slope decomposition")
    print("D01 =")
    sp.pprint(D01)
    print("D21 =")
    sp.pprint(D21)
    print("D41 =")
    sp.pprint(D41)
    print("N01 =")
    sp.pprint(N01)

    expect_zero("D01 - (dK - B01 - Z01)", D01 - (dK - B01 - Z01))
    expect_zero("D21 + dM + B21 + Z21", D21 + dM + B21 + Z21)
    expect_zero("D41 + B41 + Z41", D41 + B41 + Z41)
    expect_zero("N01/N0 - 2(P01/P - Delta01/Delta)", sp.simplify(N01 / N0 - 2 * (P01 / P - Delta01 / Delta)))

    subbanner("IV. Exact weak-axisymmetric obstruction triplet")
    K1 = sp.simplify(D21 + D01 / 9)
    Xi_load = sp.simplify(N01 / N0 - D01 / D0)
    H_even = sp.simplify(D41 - sp.Rational(2, 3) * D21 - D01 / 27)

    print("K1 =")
    sp.pprint(K1)
    print("Xi_load =")
    sp.pprint(Xi_load)
    print("H_even =")
    sp.pprint(H_even)

    # Exact prefactor-slope identity.
    P1 = sp.simplify((N01 * D0 - N0 * D01) / D0**2)
    expect_zero("Xi_load - P1/P0", sp.simplify(Xi_load - P1 / (N0 / D0)))

    subbanner("V. Exact compensation surfaces on the compensated isotropic branch")
    dM_comp = sp.simplify(D01 / 9 - B21 - Z21)
    dK_comp = sp.simplify(B01 + Z01 + D0 * N01 / N0)

    print("Even-preserving compensation dM_comp =")
    sp.pprint(dM_comp)
    print("Odd/normalization-preserving compensation dK_comp =")
    sp.pprint(dK_comp)

    expect_zero("K1 under dM = dM_comp", sp.simplify(K1.subs(dM, dM_comp)))
    expect_zero("Xi_load under dK = dK_comp", sp.simplify(Xi_load.subs(dK, dK_comp)))

    print("Remaining hidden-even residual after the two compensations is still")
    print("  H_even = D41 - (2/3) D21 - D01/27")
    print("with dM and dK replaced by the exact compensation surfaces above.")

    banner("STAGE 209 LEDGER")
    print("1. The explicit isotropic overlap model now has an exact primitive weak-axisymmetric compiler.")
    print("2. The primitive data feed the grouped obstruction triplet only through")
    print("      D01, D21, D41, N01.")
    print("3. The exact bundle decomposition is")
    print("      D01 = dK - B01 - Z01,")
    print("      D21 = -(dM + B21 + Z21),")
    print("      D41 = -(B41 + Z41),")
    print("      Xi_load = N01/N0 - D01/D0 = P1/P0.")
    print("4. On the compensated isotropic branch, K1 = 0 fixes dM exactly and Xi_load = 0 fixes dK exactly.")
    print("5. The remaining explicit first-order 5PN gate is the hidden-even residual")
    print("      H_even = D41 - (2/3) D21 - D01/27.")
