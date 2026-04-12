#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 214 — full even-gate solve on the constructive branch.

What this script does
---------------------
1. Builds the lower-bound even-gate matrix from the matched wall sector plus the
   explicit one-mode BdG sector.
2. Shows that in that lower-bound picture the mixed directions (alpha_W, alpha_R)
   remain untouched.
3. Reinstates the exact Z-sector contributions from Stage 213.
4. Forms the full constructive-slice even-gate matrix in the seven directions
      (alpha_K, alpha_W, alpha_U, alpha_R, alpha_OmegaU, beta_B, beta_varpi).
5. Proves that the mixed-sector 2x2 minor in (alpha_W, alpha_R) is nonzero, so the
   fake lower-bound freedom is removed.
6. Solves the full even-gate system exactly for (alpha_W, alpha_R).

Interpretation
--------------
This is the first executable proof that the mixed-sector freedom left by the Stage-17
lower-bound solve was never genuine. Once Z2/Z4 are reinstated, the constructive branch
solves directly for the previously untouched mixed directions.
"""


if __name__ == "__main__":
    banner("STAGE 214 — FULL EVEN-GATE SOLVE ON THE CONSTRUCTIVE BRANCH")

    alpha_K, alpha_W, alpha_U, alpha_R, alpha_OU = sp.symbols(
        "alpha_K alpha_W alpha_U alpha_R alpha_OmegaU", real=True
    )
    beta_B, beta_varpi = sp.symbols("beta_B beta_varpi", real=True)

    # ------------------------------------------------------------------
    # Lower-bound wall + BdG solve
    # ------------------------------------------------------------------
    subbanner("I. Lower-bound even-gate matrix before reinstating Z2/Z4")
    K = sp.Integer(2)
    M = sp.Integer(3)
    B0 = sp.Integer(2)
    varpi = sp.Integer(3)

    dln_M = sp.simplify(alpha_K - 2 * alpha_U + 2 * alpha_OU)

    K1_wall = sp.simplify(K * alpha_K / 9 - M * dln_M)
    H_even_wall = sp.simplify(-K * alpha_K / 27 + sp.Rational(2, 3) * M * dln_M)

    B2 = sp.simplify(B0 / varpi**2)
    B4 = sp.simplify(B0 / varpi**4)
    K1_B = sp.simplify(-2 * B2 * (beta_B - 2 * beta_varpi) - 2 * B0 * (beta_B - beta_varpi) / 9)
    H_even_B = sp.simplify(
        -2 * B4 * (beta_B - 3 * beta_varpi)
        + sp.Rational(4, 3) * B2 * (beta_B - 2 * beta_varpi)
        + 2 * B0 * (beta_B - beta_varpi) / 27
    )

    K1_lower = sp.simplify(K1_wall + K1_B)
    H_lower = sp.simplify(H_even_wall + H_even_B)

    vars7 = sp.Matrix([alpha_K, alpha_W, alpha_U, alpha_R, alpha_OU, beta_B, beta_varpi])
    Gate_lower = sp.Matrix(
        [
            [sp.diff(K1_lower, v) for v in vars7],
            [sp.diff(H_lower, v) for v in vars7],
        ]
    )

    print("Gate_lower =")
    sp.pprint(Gate_lower)
    print("rank(Gate_lower) =", Gate_lower.rank())
    expect_zero("alpha_W lower-bound column", Gate_lower[:, 1])
    expect_zero("alpha_R lower-bound column", Gate_lower[:, 3])

    # ------------------------------------------------------------------
    # Exact Z-sector reinstatement on the constructive slice
    # ------------------------------------------------------------------
    subbanner("II. Reinstating the exact Z-sector contributions")
    # Stage-211 constructive-branch parameters used in the notes.
    chi0 = sp.Rational(3, 2)
    deltaU = sp.Rational(2, 3)
    E_star = sp.Rational(1, 4)
    F_star = sp.Rational(5, 6)
    G_U = sp.Integer(5)
    G_W = sp.Integer(7)
    R = sp.Integer(2)
    Omega_U = sp.Integer(3)
    Omega_W = sp.Integer(4)

    dln_deltaU = sp.simplify(
        - (1 + deltaU) / (1 + chi0) * (alpha_R + alpha_U - alpha_W - 2 * alpha_OU)
    )
    dln_OW = sp.simplify(
        (alpha_W - alpha_U + (1 - E_star) * alpha_OU + E_star * alpha_R - (F_star / 2) * dln_deltaU)
        / (E_star + 2)
    )

    Delta = sp.simplify(Omega_U**2 * Omega_W**2 - R**2)
    Q = sp.simplify(G_U**2 * Omega_W**2 + 2 * G_U * G_W * R + G_W**2 * Omega_U**2)
    H = sp.simplify(G_U**2 + G_W**2)
    S2 = sp.simplify(Omega_U**2 + Omega_W**2)

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
            Delta**2 * H * dS2 + Delta**2 * S2 * dH + Delta**2 * dQ - 2 * Delta * H * S2 * dDelta
            - 2 * Delta * Q * S2 * dS2 - 2 * Delta * Q * dDelta - Delta * S2**2 * dQ + 3 * Q * S2**2 * dDelta
        )
        / Delta**4
    )

    K1_Z = sp.simplify(-dZ2 - dZ0 / 9)
    H_even_Z = sp.simplify(-dZ4 + sp.Rational(2, 3) * dZ2 + dZ0 / 27)

    print("K1_Z =")
    sp.pprint(K1_Z)
    print("H_even,Z =")
    sp.pprint(H_even_Z)

    # ------------------------------------------------------------------
    # Full constructive-slice solve
    # ------------------------------------------------------------------
    subbanner("III. Full constructive-slice even-gate matrix")
    K1_full = sp.simplify(K1_lower + K1_Z)
    H_full = sp.simplify(H_lower + H_even_Z)

    Gate_full = sp.Matrix(
        [
            [sp.diff(K1_full, v) for v in vars7],
            [sp.diff(H_full, v) for v in vars7],
        ]
    )

    print("Gate_full =")
    sp.pprint(Gate_full)
    print("rank(Gate_full) =", Gate_full.rank())
    if Gate_full.rank() != 2:
        raise AssertionError("Full constructive-slice gate matrix did not have rank 2.")
    print("nullity(Gate_full) =", Gate_full.shape[1] - Gate_full.rank())

    Gate_expected = sp.Matrix(
        [
            [
                -sp.Rational(25, 9),
                -sp.Rational(32134513, 50009400),
                sp.Rational(91017569, 25004700),
                sp.Rational(733046, 6251175),
                -sp.Rational(71404699, 25004700),
                -sp.Rational(8, 9),
                sp.Rational(4, 3),
            ],
            [
                sp.Rational(52, 27),
                sp.Rational(5617869293, 21003948000),
                -sp.Rational(30905427529, 10501974000),
                -sp.Rational(1174937411, 21003948000),
                sp.Rational(55109414029, 21003948000),
                sp.Rational(32, 81),
                -sp.Rational(16, 27),
            ],
        ]
    )
    expect_zero("full constructive-slice gate matrix", Gate_full - Gate_expected)

    subbanner("IV. Mixed-sector minor and exact solve for (alpha_W, alpha_R)")
    mixed_minor = sp.simplify(Gate_full[:, [1, 3]].det())
    print("det Gate_(alpha_W, alpha_R) =", mixed_minor)
    if mixed_minor == 0:
        raise AssertionError("Mixed-sector minor unexpectedly vanished.")

    solve_pair = sp.solve([sp.Eq(K1_full, 0), sp.Eq(H_full, 0)], [alpha_W, alpha_R], dict=True)
    if len(solve_pair) != 1:
        raise AssertionError("Expected a unique exact solve for (alpha_W, alpha_R).")
    sol = solve_pair[0]

    print("alpha_W =")
    sp.pprint(sp.simplify(sol[alpha_W]))
    print("alpha_R =")
    sp.pprint(sp.simplify(sol[alpha_R]))

    alphaW_expected = (
        sp.Rational(14503089433000, 942737330573) * alpha_K
        + sp.Rational(30450672110098, 942737330573) * alpha_OU
        - sp.Rational(29120459867142, 942737330573) * alpha_U
        - sp.Rational(18876066395200, 25453907925471) * beta_B
        + sp.Rational(9438033197600, 8484635975157) * beta_varpi
    )
    alphaR_expected = (
        sp.Rational(101802968743000, 942737330573) * alpha_K
        + sp.Rational(189815725996721, 942737330573) * alpha_OU
        - sp.Rational(188832473718440, 942737330573) * alpha_U
        + sp.Rational(89510801038400, 25453907925471) * beta_B
        - sp.Rational(44755400519200, 8484635975157) * beta_varpi
    )
    expect_zero("alpha_W exact solve", sp.simplify(sol[alpha_W] - alphaW_expected))
    expect_zero("alpha_R exact solve", sp.simplify(sol[alpha_R] - alphaR_expected))

    # Verify back-substitution.
    expect_zero("K1_full after exact solve", sp.simplify(K1_full.subs(sol)))
    expect_zero("H_full after exact solve", sp.simplify(H_full.subs(sol)))

    banner("STAGE 214 LEDGER")
    print("1. Before reinstating Z2/Z4, the lower-bound even-gate matrix leaves alpha_W and alpha_R")
    print("   untouched exactly.")
    print("2. After reinstating the exact Z-sector contributions, the full constructive-slice gate")
    print("   matrix still has rank 2 but its mixed-sector 2x2 minor is nonzero.")
    print("3. Therefore the previously untouched mixed directions are no longer genuine null directions.")
    print("4. On the constructive branch, the full even-gate system solves directly for alpha_W and alpha_R.")
    print("5. So the fake mixed-sector freedom in the Stage-17 lower-bound picture is removed exactly")
    print("   once the conservative Z2/Z4 block is restored.")
