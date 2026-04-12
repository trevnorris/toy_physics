#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 217 — positive-quadrant obstruction for the support-blind packet-null line.

What this script does
---------------------
1. Takes the exact affine packet-null line F_crit(E_*) from Stage 216.
2. Proves that its intercept and slope are both negative.
3. Concludes that the line never intersects the constructive positive quadrant
      E_* >= 0, F_* > 0.
4. Exhibits one mathematical null vector on the line outside that quadrant to show
   the obstruction is physical/branch-based rather than algebraic impossibility.
5. Confirms again that the representative constructive point (E_*,F_*) = (1/4,5/6)
   has only the trivial support-blind packet-null solution.

Interpretation
--------------
The tempting support-blind mixed corridor does not survive the *true* first-order
Packet-A/B endgame on the explicit positive isotropic prototype. To continue the 5PN
program honestly, one must either:

  • activate support-carrying directions again, or
  • move to a different isotropic baseline branch.
"""


def build_packet_matrix(E_star: sp.Symbol, F_star: sp.Symbol) -> sp.Matrix:
    pi = sp.pi
    I_cross = sp.simplify(sp.Rational(8, 3) / pi)

    lambda_B = sp.Integer(1)
    varpi = sp.Integer(2)
    GU = sp.Integer(1)
    GW = sp.Integer(1)
    R = sp.Rational(1, 2)
    OU = sp.Rational(3, 2)
    OW = sp.Integer(2)
    M = sp.Integer(1)

    C = sp.simplify(lambda_B * I_cross)
    Delta = sp.simplify(OU**2 * OW**2 - R**2)
    P = sp.simplify(OU**2 * GW + R * GU)
    Q = sp.simplify(GU**2 * OW**2 + 2 * GU * GW * R + GW**2 * OU**2)
    H = sp.simplify(GU**2 + GW**2)
    S2 = sp.simplify(OU**2 + OW**2)

    B0 = sp.simplify(C**2 / varpi**2)
    B2 = sp.simplify(C**2 / varpi**4)
    B4 = sp.simplify(C**2 / varpi**6)
    Z0 = sp.simplify(Q / Delta)
    Z2 = sp.simplify((Q * S2 - H * Delta) / Delta**2)
    Z4 = sp.simplify((Q * (S2**2 - Delta) - S2 * H * Delta) / Delta**3)
    N0 = sp.simplify(P**2 / Delta**2)

    K = sp.simplify((3 * (M + B2 + Z2) ** 2 + (B0 + Z0) * (B4 + Z4)) / (B4 + Z4))
    D0 = sp.simplify(K - B0 - Z0)
    D2 = sp.simplify(-(M + B2 + Z2))
    D4 = sp.simplify(-(B4 + Z4))

    alpha_K, alpha_GW, alpha_R = sp.symbols("alpha_K alpha_GW alpha_R", real=True)
    chi0 = sp.simplify(R * GU / (OU**2 * GW))
    deltaU = sp.Integer(1)
    alpha_deltaU = sp.simplify(-(1 + deltaU) / (1 + chi0) * (alpha_R - alpha_GW))
    alpha_OW = sp.simplify((alpha_GW + E_star * alpha_R - sp.Rational(1, 2) * F_star * alpha_deltaU) / (E_star + 2))

    GWs, Rs, OWs = sp.symbols("GW_s R_s OW_s", positive=True, real=True)
    prims = [GWs, Rs, OWs]
    dlogs = [alpha_GW, alpha_R, alpha_OW]

    def differential(expr: sp.Expr) -> sp.Expr:
        return sp.simplify(sum(sp.diff(expr, p) * p * dp for p, dp in zip(prims, dlogs)))

    Delta_expr = OU**2 * OWs**2 - Rs**2
    P_expr = OU**2 * GWs + Rs * GU
    Q_expr = GU**2 * OWs**2 + 2 * GU * GWs * Rs + GWs**2 * OU**2
    H_expr = GU**2 + GWs**2
    S2_expr = OU**2 + OWs**2
    subs = {GWs: GW, Rs: R, OWs: OW}

    Delta01 = sp.simplify(differential(Delta_expr).subs(subs))
    P01 = sp.simplify(differential(P_expr).subs(subs))
    Q01 = sp.simplify(differential(Q_expr).subs(subs))
    H01 = sp.simplify(differential(H_expr).subs(subs))
    S201 = sp.simplify(differential(S2_expr).subs(subs))

    Z2_num = sp.simplify(Q * S2 - H * Delta)
    Z4_num = sp.simplify(Q * (S2**2 - Delta) - S2 * H * Delta)

    Z01 = sp.simplify((Q01 * Delta - Q * Delta01) / Delta**2)
    Z21 = sp.simplify(((Q01 * S2 + Q * S201 - H01 * Delta - H * Delta01) * Delta - 2 * Z2_num * Delta01) / Delta**3)
    Z41 = sp.simplify(((Q01 * (S2**2 - Delta) + Q * (2 * S2 * S201 - Delta01) - (S201 * H * Delta + S2 * H01 * Delta + S2 * H * Delta01)) * Delta - 3 * Z4_num * Delta01) / Delta**4)
    N01 = sp.simplify(2 * N0 * (P01 / P - Delta01 / Delta))

    D01 = sp.simplify(K * alpha_K - Z01)
    D21 = sp.simplify(-(M * alpha_K + Z21))
    D41 = sp.simplify(-Z41)

    u2_1 = sp.simplify((-D0 * D21 + D2 * D01) / D0**2)
    u4_1 = sp.simplify((D0 * (-D0 * D41 + 2 * D2 * D21 - D4 * D01) + 2 * D01 * (D0 * D4 - D2**2)) / D0**3)
    Xi_load = sp.simplify(N01 / N0 - D01 / D0)

    free = [alpha_K, alpha_GW, alpha_R]
    return sp.Matrix(
        [
            [sp.simplify(sp.diff(u2_1, v)) for v in free],
            [sp.simplify(sp.diff(u4_1, v)) for v in free],
            [sp.simplify(sp.diff(Xi_load, v)) for v in free],
        ]
    )


if __name__ == "__main__":
    banner("STAGE 217 — POSITIVE-QUADRANT OBSTRUCTION FOR THE SUPPORT-BLIND PACKET-NULL LINE")

    E_star, F_star = sp.symbols("E_star F_star", real=True)
    pi = sp.pi

    coeff_E = sp.Integer(263797293760000) + sp.Integer(1757766806455275) * pi**2 + sp.Integer(3339557838723645) * pi**4 + sp.Integer(1551622258297188) * pi**6
    coeff_F = sp.Integer(48655861632000) + sp.Integer(389171318788980) * pi**2 + sp.Integer(930178880126748) * pi**4 + sp.Integer(694451446430976) * pi**6
    coeff_0 = sp.Integer(102703468791960) * pi**6 - sp.Integer(155911749769062) * pi**4 - sp.Integer(147002028439770) * pi**2 - sp.Integer(25886544768000)

    F_crit = sp.simplify(-(coeff_E * E_star + coeff_0) / coeff_F)
    slope = sp.simplify(sp.diff(F_crit, E_star))
    intercept = sp.simplify(F_crit.subs(E_star, 0))

    subbanner("I. Exact affine line and sign data")
    print("F_crit(E_star) =")
    sp.pprint(F_crit)
    print("dF_crit/dE_star =")
    sp.pprint(slope)
    print("F_crit(0) =")
    sp.pprint(intercept)
    print("Numerically:")
    print("  slope      ≈", sp.N(slope))
    print("  intercept  ≈", sp.N(intercept))
    print("  F_crit(1/4)≈", sp.N(F_crit.subs(E_star, sp.Rational(1, 4))))
    print("  F_crit(1/2)≈", sp.N(F_crit.subs(E_star, sp.Rational(1, 2))))
    print("  F_crit(1)  ≈", sp.N(F_crit.subs(E_star, sp.Integer(1))))

    # Since coeff_E, coeff_F, coeff_0 are all positive, slope and intercept are both negative.
    if not (float(sp.N(coeff_E)) > 0 and float(sp.N(coeff_F)) > 0 and float(sp.N(coeff_0)) > 0):
        raise AssertionError("Expected positive affine-line coefficients.")
    if not (float(sp.N(slope)) < 0 and float(sp.N(intercept)) < 0):
        raise AssertionError("Expected negative slope and intercept for F_crit.")

    subbanner("II. The constructive positive quadrant does not intersect the line")
    print("Because coeff_E, coeff_F, coeff_0 > 0, we have")
    print("  F_crit(E_star) < 0  for every E_star >= 0.")
    print("So the support-blind packet-null line never enters the constructive quadrant")
    print("  E_star >= 0,  F_star > 0.")

    subbanner("III. Mathematical null vector exists on the line outside the quadrant")
    M_packet = build_packet_matrix(E_star, F_star)
    E_sample = sp.Rational(1, 4)
    F_sample = sp.simplify(F_crit.subs(E_star, E_sample))
    M_line = sp.simplify(M_packet.subs({E_star: E_sample, F_star: F_sample}))
    print("rank on the line at E_star = 1/4 =", M_line.rank())
    if M_line.rank() != 2:
        raise AssertionError("Expected a 1-parameter null corridor on the packet-null line.")
    null_vec = M_line.nullspace()[0]
    null_vec = sp.simplify(null_vec / null_vec[2])
    print("Normalized null vector on the line (alpha_K, alpha_GW, alpha_R) with alpha_R = 1 =")
    sp.pprint(null_vec)
    print("Numerically =")
    sp.pprint(sp.Matrix([sp.N(x) for x in null_vec]))

    subbanner("IV. The representative constructive point remains trivial")
    M_constructive = sp.simplify(M_packet.subs({E_star: sp.Rational(1, 4), F_star: sp.Rational(5, 6)}))
    print("rank at (E_*,F_*) = (1/4, 5/6) =", M_constructive.rank())
    if M_constructive.rank() != 3:
        raise AssertionError("Expected the constructive positive sample to remain full rank.")

    banner("STAGE 217 LEDGER")
    print("1. The exact support-blind packet-null line is affine:")
    print("      F_* = F_crit(E_*).")
    print("2. Its slope and intercept are both negative, so it never intersects the")
    print("   constructive quadrant E_* >= 0, F_* > 0.")
    print("3. Therefore the support-blind mixed sector does admit a mathematical first-order")
    print("   packet-null corridor, but only outside the positive constructive-orbit branch.")
    print("4. On the physically motivated constructive sample (E_*,F_*) = (1/4, 5/6),")
    print("   the support-blind packet-null system is full rank and therefore trivial.")
    print("5. The next honest 5PN continuation is to restore support-carrying directions")
    print("   (alpha_GU, alpha_OU, beta_B, beta_varpi), or else move to a different")
    print("   isotropic baseline branch before re-testing the first-order endgame.")
