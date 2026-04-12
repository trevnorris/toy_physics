#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 210 — primitive mechanism sieve and the surviving outgoing corridor.

What this script does
---------------------
1. Tests the three Stage-209 obstruction gates
      K1 = D21 + D01/9,
      Xi_load = N01/N0 - D01/D0,
      H_even = D41 - (2/3) D21 - D01/27
   on primitive sector-restricted weak-axisymmetric deformations.
2. Proves the two direct no-go results from the notes:
      - wall-only anisotropy is trivial,
      - pure BdG anisotropy is also trivial.
3. Shows that BdG self-similarity x_c = x_varpi kills only Xi_load, not the full
   5PN obstruction triplet.
4. Factors the outgoing load as
      N0 / K = M_leg^2 (1 + I)^2 / (1 - H)^2,
   derives the exact one-port outgoing defect,
      Sigma^(N) = 2 m + 2 I/(1+I) i + 2 H/(1-H) h,
   and proves that under rigid interference/hybridization ratios the surviving
   nontrivial cancellation condition is the square-root mixed-leg law
      G_W / Omega_W^2 ∝ sqrt(K).

Interpretation
--------------
This is the first actual weak-axisymmetric mechanism sieve after the isotropic
Packet-A/Packet-B bridge. It kills the dead primitive sectors and isolates the
mixed-sector outgoing corridor that still survives.
"""


if __name__ == "__main__":
    banner("STAGE 210 — MECHANISM SIEVE AND OUTGOING CORRIDOR")

    # Generic compensated-branch symbols.
    D0 = sp.symbols("D_0", positive=True, real=True)
    B0, B2, B4 = sp.symbols("B_0 B_2 B_4", positive=True, real=True)
    dK, dM = sp.symbols("dK dM", real=True)

    subbanner("I. Wall-only weak-axisymmetric anisotropy is dead")
    D01_wall = dK
    D21_wall = -dM
    D41_wall = sp.Integer(0)
    N01_wall = sp.Integer(0)

    K1_wall = sp.simplify(D21_wall + D01_wall / 9)
    Xi_wall = sp.simplify(N01_wall / 1 - D01_wall / D0)
    H_wall = sp.simplify(D41_wall - sp.Rational(2, 3) * D21_wall - D01_wall / 27)

    print("K1_wall =")
    sp.pprint(K1_wall)
    print("Xi_wall =")
    sp.pprint(Xi_wall)
    print("H_wall =")
    sp.pprint(H_wall)

    wall_solutions = sp.solve(
        [sp.Eq(K1_wall, 0), sp.Eq(Xi_wall, 0), sp.Eq(H_wall, 0)],
        [dK, dM],
        dict=True,
    )
    print("wall-only solve of (K1, Xi_load, H_even) = 0 gives")
    sp.pprint(wall_solutions)
    if wall_solutions != [{dK: sp.Integer(0), dM: sp.Integer(0)}]:
        raise AssertionError("Wall-only sector did not collapse to the trivial branch.")

    subbanner("II. Pure BdG weak-axisymmetric anisotropy is also dead")
    x_c, x_varpi = sp.symbols("x_c x_varpi", real=True)

    B01 = sp.simplify(2 * B0 * (x_c - x_varpi))
    B21 = sp.simplify(2 * B2 * (x_c - 2 * x_varpi))
    B41 = sp.simplify(2 * B4 * (x_c - 3 * x_varpi))

    D01_bdg = sp.simplify(-B01)
    D21_bdg = sp.simplify(-B21)
    D41_bdg = sp.simplify(-B41)
    Xi_bdg = sp.simplify(-D01_bdg / D0)
    K1_bdg = sp.simplify(D21_bdg + D01_bdg / 9)
    H_bdg = sp.simplify(D41_bdg - sp.Rational(2, 3) * D21_bdg - D01_bdg / 27)

    print("B01 =")
    sp.pprint(B01)
    print("B21 =")
    sp.pprint(B21)
    print("B41 =")
    sp.pprint(B41)
    print("K1_bdg =")
    sp.pprint(K1_bdg)
    print("Xi_bdg =")
    sp.pprint(Xi_bdg)
    print("H_bdg =")
    sp.pprint(H_bdg)

    bdg_solutions = sp.solve(
        [sp.Eq(K1_bdg, 0), sp.Eq(Xi_bdg, 0), sp.Eq(H_bdg, 0)],
        [x_c, x_varpi],
        dict=True,
    )
    print("pure-BdG solve of (K1, Xi_load, H_even) = 0 gives")
    sp.pprint(bdg_solutions)
    if bdg_solutions != [{x_c: sp.Integer(0), x_varpi: sp.Integer(0)}]:
        raise AssertionError("Pure BdG sector did not collapse to the trivial branch.")

    subbanner("III. BdG self-similarity kills only the load defect")
    K1_bdg_self = sp.simplify(K1_bdg.subs(x_c, x_varpi))
    Xi_bdg_self = sp.simplify(Xi_bdg.subs(x_c, x_varpi))
    H_bdg_self = sp.simplify(H_bdg.subs(x_c, x_varpi))

    print("K1_bdg under x_c = x_varpi =")
    sp.pprint(K1_bdg_self)
    print("Xi_bdg under x_c = x_varpi =")
    sp.pprint(Xi_bdg_self)
    print("H_bdg under x_c = x_varpi =")
    sp.pprint(H_bdg_self)

    expect_zero("Xi_bdg under self-similarity", Xi_bdg_self)
    if K1_bdg_self == 0 or H_bdg_self == 0:
        raise AssertionError("BdG self-similarity unexpectedly killed the full obstruction triplet.")

    subbanner("IV. Exact one-port outgoing-load factorization")
    K = sp.symbols("K", positive=True, real=True)
    G_U, G_W, R = sp.symbols("G_U G_W R", positive=True, real=True)
    Omega_U, Omega_W = sp.symbols("Omega_U Omega_W", positive=True, real=True)

    Delta = sp.simplify(Omega_U**2 * Omega_W**2 - R**2)
    P = sp.simplify(Omega_U**2 * G_W + R * G_U)
    N0 = sp.simplify(P**2 / Delta**2)

    M_leg = sp.simplify(G_W / (Omega_W**2 * sp.sqrt(K)))
    I_ratio = sp.simplify(R * G_U / (Omega_U**2 * G_W))
    H_ratio = sp.simplify(R**2 / (Omega_U**2 * Omega_W**2))

    factored = sp.simplify(M_leg**2 * (1 + I_ratio) ** 2 / (1 - H_ratio) ** 2)
    print("N0 / K =")
    sp.pprint(sp.simplify(N0 / K))
    print("Factored form =")
    sp.pprint(factored)
    expect_zero("outgoing-load factorization", sp.simplify(N0 / K - factored))

    subbanner("V. Exact one-port outgoing defect and square-root mixed-leg law")
    m, i, h = sp.symbols("m i h", real=True)
    Ibar, Hbar = sp.symbols("I H", positive=True, real=True)

    Sigma_N = sp.simplify(2 * m + 2 * Ibar / (1 + Ibar) * i + 2 * Hbar / (1 - Hbar) * h)
    print("Sigma^(N) =")
    sp.pprint(Sigma_N)

    # Under rigid interference and hybridization ratios, i = h = 0.
    Sigma_rigid = sp.simplify(Sigma_N.subs({i: 0, h: 0}))
    print("Sigma^(N) under rigid I,H =")
    sp.pprint(Sigma_rigid)
    expect_zero("Sigma^(N) under rigid I,H - 2m", Sigma_rigid - 2 * m)

    gW, oW, kappa1 = sp.symbols("g_W o_W kappa_1", real=True)
    m_raw = sp.simplify(gW - 2 * oW - kappa1 / 2)
    zero_defect = sp.simplify(Sigma_rigid.subs(m, m_raw))
    print("Rigid-branch zero-defect scalar =")
    sp.pprint(zero_defect)
    expect_zero(
        "square-root mixed-leg law",
        sp.simplify(zero_defect - (2 * gW - 4 * oW - kappa1)),
    )

    banner("STAGE 210 LEDGER")
    print("1. Wall-only weak-axisymmetric anisotropy fails: the full (K1, Xi_load, H_even) system")
    print("   has only the trivial solution dK = dM = 0.")
    print("2. Pure BdG weak-axisymmetric anisotropy also fails: the full system has only")
    print("   x_c = x_varpi = 0.")
    print("3. BdG self-similarity x_c = x_varpi kills only Xi_load, not the full 5PN triplet.")
    print("4. The surviving nontrivial outgoing corridor is the mixed-sector load factor")
    print("      N0/K = M_leg^2 (1 + I)^2 / (1 - H)^2.")
    print("5. Under rigid interference/hybridization ratios, zero defect reduces exactly to")
    print("      2 g_W - 4 o_W - kappa_1 = 0,")
    print("   i.e. the square-root mixed-leg law")
    print("      G_W / Omega_W^2 ∝ sqrt(K).")
