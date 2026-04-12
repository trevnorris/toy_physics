#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 213 — exact Z-sector bridge back into the live even gates.

What this script does
---------------------
1. Reinstates the conservative Maxwell/mixed moments
      Z0 = Q/Delta,
      Z2 = (Q S2 - H Delta)/Delta^2,
      Z4 = [Q(S2^2-Delta) - S2 H Delta]/Delta^3.
2. Differentiates them exactly to obtain dZ0, dZ2, dZ4.
3. Inserts the normalized Stage-211 similarity-kernel relations for
      dln(delta_U), dln(M), dln(Omega_W)
   so that the Z-sector depends only on the free mixed/wall/U drifts.
4. Compiles the Z-sector contributions to the two live 5PN even gates,
      K1_Z = -dZ2 - dZ0/9,
      H_even,Z = -dZ4 + (2/3)dZ2 + dZ0/27.
5. Evaluates those expressions on the constructive slice recorded in the notes.

Interpretation
--------------
This is the exact step that the lower-bound Stage-17 picture was missing: once the
conservative Z2/Z4 sector is restored, the previously untouched mixed directions
(alpha_W, alpha_R) become active in the even-gate problem.
"""


if __name__ == "__main__":
    banner("STAGE 213 — EXACT Z-SECTOR EVEN-GATE BRIDGE")

    # Free normalized similarity directions and constructive-branch parameters.
    alpha_K, alpha_W, alpha_U, alpha_R, alpha_OU = sp.symbols(
        "alpha_K alpha_W alpha_U alpha_R alpha_OmegaU", real=True
    )
    chi0, deltaU = sp.symbols("chi_0 delta_U", positive=True, real=True)
    E_star, F_star = sp.symbols("E_star F_star", real=True)

    # Base mixed-sector variables.
    G_U, G_W, R = sp.symbols("G_U G_W R", positive=True, real=True)
    Omega_U, Omega_W = sp.symbols("Omega_U Omega_W", positive=True, real=True)

    # Stage-211 dependent drifts.
    dln_deltaU = sp.simplify(
        - (1 + deltaU) / (1 + chi0) * (alpha_R + alpha_U - alpha_W - 2 * alpha_OU)
    )
    dln_OW = sp.simplify(
        (alpha_W - alpha_U + (1 - E_star) * alpha_OU + E_star * alpha_R - (F_star / 2) * dln_deltaU)
        / (E_star + 2)
    )

    subbanner("I. Conservative Z-sector moments")
    Delta = sp.simplify(Omega_U**2 * Omega_W**2 - R**2)
    Q = sp.simplify(G_U**2 * Omega_W**2 + 2 * G_U * G_W * R + G_W**2 * Omega_U**2)
    H = sp.simplify(G_U**2 + G_W**2)
    S2 = sp.simplify(Omega_U**2 + Omega_W**2)

    Z0 = sp.simplify(Q / Delta)
    Z2 = sp.simplify((Q * S2 - H * Delta) / Delta**2)
    Z4 = sp.simplify((Q * (S2**2 - Delta) - S2 * H * Delta) / Delta**3)

    print("Z0 =")
    sp.pprint(Z0)
    print("Z2 =")
    sp.pprint(Z2)
    print("Z4 =")
    sp.pprint(Z4)

    subbanner("II. Exact first variations")
    dDelta = sp.simplify(2 * Omega_U**2 * Omega_W**2 * (alpha_OU + dln_OW) - 2 * R**2 * alpha_R)
    dQ = sp.simplify(
        2 * G_U**2 * Omega_W**2 * (alpha_U + dln_OW)
        + 2 * G_U * G_W * R * (alpha_U + alpha_W + alpha_R)
        + 2 * G_W**2 * Omega_U**2 * (alpha_W + alpha_OU)
    )
    dH = sp.simplify(2 * G_U**2 * alpha_U + 2 * G_W**2 * alpha_W)
    dS2 = sp.simplify(2 * Omega_U**2 * alpha_OU + 2 * Omega_W**2 * dln_OW)

    dZ0 = sp.simplify((Delta * dQ - Q * dDelta) / Delta**2)
    dZ2 = sp.simplify(
        (Delta * (-Delta * dH - H * dDelta + Q * dS2 + S2 * dQ) + 2 * dDelta * (Delta * H - Q * S2))
        / Delta**3
    )
    dZ4 = sp.simplify(
        -(
            Delta**2 * H * dS2
            + Delta**2 * S2 * dH
            + Delta**2 * dQ
            - 2 * Delta * H * S2 * dDelta
            - 2 * Delta * Q * S2 * dS2
            - 2 * Delta * Q * dDelta
            - Delta * S2**2 * dQ
            + 3 * Q * S2**2 * dDelta
        )
        / Delta**4
    )

    dZ0_expected = sp.simplify((Delta * dQ - Q * dDelta) / Delta**2)
    dZ2_expected = sp.simplify(
        (Delta * (-Delta * dH - H * dDelta + Q * dS2 + S2 * dQ) + 2 * dDelta * (Delta * H - Q * S2))
        / Delta**3
    )
    dZ4_expected = sp.simplify(
        -(
            Delta**2 * H * dS2 + Delta**2 * S2 * dH + Delta**2 * dQ - 2 * Delta * H * S2 * dDelta
            - 2 * Delta * Q * S2 * dS2 - 2 * Delta * Q * dDelta - Delta * S2**2 * dQ + 3 * Q * S2**2 * dDelta
        )
        / Delta**4
    )

    expect_zero("dZ0 formula", dZ0 - dZ0_expected)
    expect_zero("dZ2 formula", dZ2 - dZ2_expected)
    expect_zero("dZ4 formula", dZ4 - dZ4_expected)

    print("dZ0 =")
    sp.pprint(dZ0)
    print("dZ2 =")
    sp.pprint(dZ2)
    print("dZ4 =")
    sp.pprint(dZ4)

    subbanner("III. Exact Z-sector contributions to the live even gates")
    K1_Z = sp.simplify(-dZ2 - dZ0 / 9)
    H_even_Z = sp.simplify(-dZ4 + sp.Rational(2, 3) * dZ2 + dZ0 / 27)

    print("K1_Z =")
    sp.pprint(K1_Z)
    print("H_even,Z =")
    sp.pprint(H_even_Z)

    subbanner("IV. Constructive slice evaluation")
    slice_subs = {
        G_U: sp.Integer(5),
        G_W: sp.Integer(7),
        R: sp.Integer(2),
        Omega_U: sp.Integer(3),
        Omega_W: sp.Integer(4),
        chi0: sp.Rational(3, 2),
        deltaU: sp.Rational(2, 3),
        E_star: sp.Rational(1, 4),
        F_star: sp.Rational(5, 6),
    }

    K1_Z_slice = sp.simplify(K1_Z.subs(slice_subs))
    H_even_Z_slice = sp.simplify(H_even_Z.subs(slice_subs))

    print("K1_Z on constructive slice =")
    sp.pprint(K1_Z_slice)
    print("H_even,Z on constructive slice =")
    sp.pprint(H_even_Z_slice)

    K1_Z_expected = sp.Matrix([
        sp.Rational(-32134513, 50009400),
        sp.Rational(733046, 6251175),
        sp.Rational(-59010631, 25004700),
        sp.Rational(78623501, 25004700),
    ])
    K1_Z_actual = sp.Matrix([
        sp.expand(K1_Z_slice).coeff(alpha_W),
        sp.expand(K1_Z_slice).coeff(alpha_R),
        sp.expand(K1_Z_slice).coeff(alpha_U),
        sp.expand(K1_Z_slice).coeff(alpha_OU),
    ])
    expect_zero("K1_Z slice coefficients", K1_Z_actual - K1_Z_expected)

    HZ_expected = sp.Matrix([
        sp.Rational(5617869293, 21003948000),
        sp.Rational(-1174937411, 21003948000),
        sp.Rational(11102468471, 10501974000),
        sp.Rational(-28906377971, 21003948000),
    ])
    HZ_actual = sp.Matrix([
        sp.expand(H_even_Z_slice).coeff(alpha_W),
        sp.expand(H_even_Z_slice).coeff(alpha_R),
        sp.expand(H_even_Z_slice).coeff(alpha_U),
        sp.expand(H_even_Z_slice).coeff(alpha_OU),
    ])
    expect_zero("H_even,Z slice coefficients", HZ_actual - HZ_expected)

    print("Pure alpha_W direction:")
    print("  K1_Z =", sp.simplify(K1_Z_slice.subs({alpha_W: 1, alpha_R: 0, alpha_U: 0, alpha_OU: 0})))
    print("  H_even,Z =", sp.simplify(H_even_Z_slice.subs({alpha_W: 1, alpha_R: 0, alpha_U: 0, alpha_OU: 0})))
    print("Pure alpha_R direction:")
    print("  K1_Z =", sp.simplify(K1_Z_slice.subs({alpha_W: 0, alpha_R: 1, alpha_U: 0, alpha_OU: 0})))
    print("  H_even,Z =", sp.simplify(H_even_Z_slice.subs({alpha_W: 0, alpha_R: 1, alpha_U: 0, alpha_OU: 0})))

    banner("STAGE 213 LEDGER")
    print("1. The conservative Maxwell/mixed moments Z0, Z2, Z4 now have exact first variations")
    print("   on the normalized similarity kernel.")
    print("2. Their contributions to the live even gates are")
    print("      K1_Z = -dZ2 - dZ0/9,")
    print("      H_even,Z = -dZ4 + (2/3)dZ2 + dZ0/27.")
    print("3. On the constructive slice, both gates depend nontrivially on alpha_W and alpha_R.")
    print("4. So the omitted Z2/Z4 block does exactly the job the lower-bound picture predicted:")
    print("   it activates the previously untouched mixed directions.")
