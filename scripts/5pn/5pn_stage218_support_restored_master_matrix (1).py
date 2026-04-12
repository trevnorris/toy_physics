#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 218 — support-restored first-order packet master matrix.

What this script does
---------------------
1. Reuses the same explicit positive isotropic finite-throat prototype from Stages
   215–217.
2. Restores the support-carrying directions
      (alpha_GU, alpha_OU, beta_B, beta_varpi)
   in addition to the mixed/wall directions
      (alpha_K, alpha_GW, alpha_R).
3. Builds the *true* first-order Packet-A/B master matrix for
      u2^(1), u4^(1), Xi_load,
   with the orbit-lock drifts (alpha_deltaU, alpha_OW, alpha_M) imposed exactly.
4. Proves that a support-only 3x3 minor is already nonzero throughout the positive
   constructive quadrant E_* >= 0, F_* > 0.
5. Concludes that the current isotropic baseline is not dead: the full master matrix
   has rank 3 and nullity 4 there, while the support-only slice already has nullity 2.

Interpretation
--------------
The Stage-217 obstruction killed only the *support-blind* corridor. It did not force
an isotropic-baseline pivot. The next honest move is therefore to restore the support
variables and solve the packet master matrix on the current branch.
"""


def build_master_system(E_star: sp.Symbol | sp.Expr, F_star: sp.Symbol | sp.Expr):
    """Return the explicit isotropic prototype data and the full first-order packet master matrix.

    Free directions are ordered as
      (alpha_K, alpha_GW, alpha_GU, alpha_R, alpha_OU, beta_B, beta_varpi).
    """
    pi = sp.pi

    # Explicit positive isotropic finite-throat prototype.
    I_cross = sp.simplify(sp.Rational(8, 3) / pi)
    lambda_B = sp.Integer(1)
    varpi = sp.Integer(2)
    G_U = sp.Integer(1)
    G_W = sp.Integer(1)
    R = sp.Rational(1, 2)
    Omega_U = sp.Rational(3, 2)
    Omega_W = sp.Integer(2)
    M = sp.Integer(1)

    C = sp.simplify(lambda_B * I_cross)
    Delta = sp.simplify(Omega_U**2 * Omega_W**2 - R**2)
    P = sp.simplify(Omega_U**2 * G_W + R * G_U)
    Q = sp.simplify(G_U**2 * Omega_W**2 + 2 * G_U * G_W * R + G_W**2 * Omega_U**2)
    H = sp.simplify(G_U**2 + G_W**2)
    S2 = sp.simplify(Omega_U**2 + Omega_W**2)

    B0 = sp.simplify(C**2 / varpi**2)
    B2 = sp.simplify(C**2 / varpi**4)
    B4 = sp.simplify(C**2 / varpi**6)
    Z0 = sp.simplify(Q / Delta)
    Z2 = sp.simplify((Q * S2 - H * Delta) / Delta**2)
    Z4 = sp.simplify((Q * (S2**2 - Delta) - S2 * H * Delta) / Delta**3)
    N0 = sp.simplify(P**2 / Delta**2)

    # Fix K on the exact isotropic one-pole surface.
    K = sp.simplify((3 * (M + B2 + Z2) ** 2 + (B0 + Z0) * (B4 + Z4)) / (B4 + Z4))
    D0 = sp.simplify(K - B0 - Z0)
    D2 = sp.simplify(-(M + B2 + Z2))
    D4 = sp.simplify(-(B4 + Z4))

    # Orbit-lock-consistent logarithmic drifts.
    alpha_K, alpha_GW, alpha_GU, alpha_R, alpha_OU = sp.symbols(
        "alpha_K alpha_GW alpha_GU alpha_R alpha_OmegaU", real=True
    )
    beta_B, beta_varpi = sp.symbols("beta_B beta_varpi", real=True)

    chi0 = sp.simplify(R * G_U / (Omega_U**2 * G_W))
    deltaU = sp.Integer(1)
    alpha_deltaU = sp.simplify(
        -(1 + deltaU) / (1 + chi0) * (alpha_R + alpha_GU - alpha_GW - 2 * alpha_OU)
    )
    alpha_OW = sp.simplify(
        (
            alpha_GW
            - alpha_GU
            + (1 - E_star) * alpha_OU
            + E_star * alpha_R
            - sp.Rational(1, 2) * F_star * alpha_deltaU
        )
        / (E_star + 2)
    )
    alpha_M = sp.simplify(alpha_K - 2 * alpha_GU + 2 * alpha_OU)

    # Logarithmic primitive variations of the mixed block.
    GWs, GUs, Rs, OUs, OWs = sp.symbols("GW_s GU_s R_s OU_s OW_s", positive=True, real=True)
    prims = [GWs, GUs, Rs, OUs, OWs]
    dlogs = [alpha_GW, alpha_GU, alpha_R, alpha_OU, alpha_OW]

    def differential_log(expr: sp.Expr) -> sp.Expr:
        return sp.simplify(sum(sp.diff(expr, p) * p * dp for p, dp in zip(prims, dlogs)))

    Delta_expr = OUs**2 * OWs**2 - Rs**2
    P_expr = OUs**2 * GWs + Rs * GUs
    Q_expr = GUs**2 * OWs**2 + 2 * GUs * GWs * Rs + GWs**2 * OUs**2
    H_expr = GUs**2 + GWs**2
    S2_expr = OUs**2 + OWs**2
    subs = {GWs: G_W, GUs: G_U, Rs: R, OUs: Omega_U, OWs: Omega_W}

    Delta01 = sp.simplify(differential_log(Delta_expr).subs(subs))
    P01 = sp.simplify(differential_log(P_expr).subs(subs))
    Q01 = sp.simplify(differential_log(Q_expr).subs(subs))
    H01 = sp.simplify(differential_log(H_expr).subs(subs))
    S201 = sp.simplify(differential_log(S2_expr).subs(subs))

    B01 = sp.simplify(2 * B0 * (beta_B - beta_varpi))
    B21 = sp.simplify(2 * B2 * (beta_B - 2 * beta_varpi))
    B41 = sp.simplify(2 * B4 * (beta_B - 3 * beta_varpi))

    Z2_num = sp.simplify(Q * S2 - H * Delta)
    Z4_num = sp.simplify(Q * (S2**2 - Delta) - S2 * H * Delta)

    Z01 = sp.simplify((Q01 * Delta - Q * Delta01) / Delta**2)
    Z21 = sp.simplify(
        ((Q01 * S2 + Q * S201 - H01 * Delta - H * Delta01) * Delta - 2 * Z2_num * Delta01)
        / Delta**3
    )
    Z41 = sp.simplify(
        (
            (
                Q01 * (S2**2 - Delta)
                + Q * (2 * S2 * S201 - Delta01)
                - (S201 * H * Delta + S2 * H01 * Delta + S2 * H * Delta01)
            )
            * Delta
            - 3 * Z4_num * Delta01
        )
        / Delta**4
    )
    N01 = sp.simplify(2 * N0 * (P01 / P - Delta01 / Delta))

    D01 = sp.simplify(K * alpha_K - B01 - Z01)
    D21 = sp.simplify(-(M * alpha_M + B21 + Z21))
    D41 = sp.simplify(-(B41 + Z41))

    u2_1 = sp.simplify((-D0 * D21 + D2 * D01) / D0**2)
    u4_1 = sp.simplify(
        (D0 * (-D0 * D41 + 2 * D2 * D21 - D4 * D01) + 2 * D01 * (D0 * D4 - D2**2)) / D0**3
    )
    Xi_load = sp.simplify(N01 / N0 - D01 / D0)

    free = [alpha_K, alpha_GW, alpha_GU, alpha_R, alpha_OU, beta_B, beta_varpi]
    M_master = sp.Matrix(
        [
            [sp.simplify(sp.diff(u2_1, v)) for v in free],
            [sp.simplify(sp.diff(u4_1, v)) for v in free],
            [sp.simplify(sp.diff(Xi_load, v)) for v in free],
        ]
    )

    return {
        "free": free,
        "M_master": M_master,
        "chi0": chi0,
        "u2": sp.simplify(-D2 / D0),
        "D0": D0,
        "D2": D2,
        "D4": D4,
        "u2_1": u2_1,
        "u4_1": u4_1,
        "Xi_load": Xi_load,
    }


if __name__ == "__main__":
    banner("STAGE 218 — SUPPORT-RESTORED FIRST-ORDER PACKET MASTER MATRIX")

    E_star, F_star = sp.symbols("E_star F_star", real=True)
    data = build_master_system(E_star, F_star)
    M_master = data["M_master"]

    subbanner("I. Full first-order Packet-A/B master matrix")
    print("Free directions are ordered as")
    print("  (alpha_K, alpha_GW, alpha_GU, alpha_R, alpha_OmegaU, beta_B, beta_varpi)")
    print("M_master =")
    sp.pprint(M_master)

    subbanner("II. Support-only minor in (alpha_K, alpha_GU, alpha_OmegaU)")
    M_support = M_master[:, [0, 2, 4, 5, 6]]
    minor_support = sp.factor(sp.together(sp.Matrix(M_support[:, [0, 1, 2]]).det()))
    print("det M_support[:, (alpha_K,alpha_GU,alpha_OmegaU)] =")
    sp.pprint(minor_support)

    factor_core = sp.expand(
        sp.simplify(
            -minor_support * 3961650000 * (490 + 1503 * sp.pi**2) ** 6 * (E_star + 2)
            / (sp.pi**2 * (8575 + 12717 * sp.pi**2) ** 2)
        )
    )
    print("Core positive-affine factor =")
    sp.pprint(factor_core)

    coeff_E = sp.simplify(sp.diff(factor_core, E_star))
    coeff_F = sp.simplify(sp.diff(factor_core, F_star))
    coeff_0 = sp.simplify(factor_core.subs({E_star: 0, F_star: 0}))
    print("coeff_E =")
    sp.pprint(coeff_E)
    print("coeff_F =")
    sp.pprint(coeff_F)
    print("coeff_0 =")
    sp.pprint(coeff_0)

    # Positive-quadrant sign check: all three affine coefficients are positive.
    if not (float(sp.N(coeff_E)) > 0 and float(sp.N(coeff_F)) > 0 and float(sp.N(coeff_0)) > 0):
        raise AssertionError("Expected positive affine coefficients in the support-only minor.")

    print("Numerically:")
    print("  coeff_E ≈", sp.N(coeff_E))
    print("  coeff_F ≈", sp.N(coeff_F))
    print("  coeff_0 ≈", sp.N(coeff_0))

    subbanner("III. Rank / nullity verdict")
    print("Because the support-only 3x3 minor is nonzero for every E_* >= 0, F_* >= 0, we have:")
    print("  rank(M_support) = 3")
    print("  rank(M_master)  = 3")
    print("  nullity(M_support) = 2")
    print("  nullity(M_master)  = 4")

    # Concrete constructive check.
    constructive = {E_star: sp.Rational(1, 4), F_star: sp.Rational(5, 6)}
    M_master_c = sp.simplify(M_master.subs(constructive))
    M_support_c = sp.simplify(M_support.subs(constructive))
    print("rank(M_support) at (E_*,F_*)=(1/4,5/6) =", M_support_c.rank())
    print("rank(M_master)  at (E_*,F_*)=(1/4,5/6) =", M_master_c.rank())
    if M_support_c.rank() != 3 or M_master_c.rank() != 3:
        raise AssertionError("Unexpected rank drop on the constructive branch.")

    banner("STAGE 218 LEDGER")
    print("1. Restoring the support carriers gives a 3x7 first-order Packet-A/B master matrix on the")
    print("   same explicit isotropic baseline that had killed the support-blind corridor.")
    print("2. A support-only 3x3 minor in (alpha_K, alpha_GU, alpha_OmegaU) is already nonzero")
    print("   throughout the positive constructive quadrant E_* >= 0, F_* >= 0.")
    print("3. Therefore the support-only slice has nullity 2 and the full support-restored master")
    print("   system has nullity 4 on the current isotropic baseline.")
    print("4. So the Stage-217 obstruction does not force an isotropic-baseline pivot. It only kills")
    print("   the support-blind corridor.")
    print("5. The next honest move is to solve the support-restored master matrix, not to pivot the")
    print("   isotropic baseline prematurely.")
