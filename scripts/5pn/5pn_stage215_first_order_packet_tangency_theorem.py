#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage199_201_common import banner, subbanner, expect_zero, grouped_trace_anomaly

"""
Stage 215 — exact first-order Packet-A tangency theorem on the weak-axisymmetric branch.

What this script does
---------------------
1. Starts from the exact grouped weak-axisymmetric signature
      lambda_(20)=1, lambda_(21)=1/2, lambda_(22)=-1.
2. Perturbs the grouped-lane conservative/outgoing data
      (D_A0, D_A2, D_A4, N_A0)
   around an isotropic baseline and computes the first-order grouped response
   moments
      u2^(A), u4^(A), P0^(A).
3. Proves the exact grouped-anomaly formulas
      a = (epsilon/4) s,   b = (3 epsilon/4) s
   for any weak-axisymmetric scalar slope s.
4. Derives the exact true packet-tangency conditions in operator form:
      u2^(1) = 0,
      u4^(1) = 0,
      P1     = 0.
5. Shows how the Stage-209 surrogate gates (K1, H_even, Xi_load) arise as an
   exact reparameterization of those true packet conditions on the canonical
   compensated branch u2 = 1/9.

Interpretation
--------------
This script cleanly separates two notions that had been easy to blur together:

  • true first-order Packet-A tangency of the grouped response,
  • the surrogate even-gate combinations (K1, H_even).

They agree exactly on the canonical compensated branch, but not on a generic
isotropic prototype with a different baseline u2.
"""


if __name__ == "__main__":
    banner("STAGE 215 — EXACT FIRST-ORDER PACKET-A TANGENCY THEOREM")

    eps = sp.symbols("epsilon", real=True)
    D0, D2, D4, N0 = sp.symbols("D0 D2 D4 N0", nonzero=True, real=True)
    D01, D21, D41, N01 = sp.symbols("D01 D21 D41 N01", real=True)

    # Exact isotropic baseline response/prefactor moments.
    u2 = sp.simplify(-D2 / D0)
    u4 = sp.simplify((D2**2 - D0 * D4) / D0**2)
    P0 = sp.simplify(N0 / D0)

    lam = {
        "20": sp.Integer(1),
        "21": sp.Rational(1, 2),
        "22": -sp.Integer(1),
    }

    def lane_response(lmbd: sp.Expr) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
        D0A = sp.simplify(D0 + eps * lmbd * D01)
        D2A = sp.simplify(D2 + eps * lmbd * D21)
        D4A = sp.simplify(D4 + eps * lmbd * D41)
        N0A = sp.simplify(N0 + eps * lmbd * N01)

        u2A = sp.simplify(-D2A / D0A)
        u4A = sp.simplify((D2A**2 - D0A * D4A) / D0A**2)
        P0A = sp.simplify(N0A / D0A)
        return u2A, u4A, P0A

    u2_20, u4_20, P0_20 = lane_response(lam["20"])
    u2_21, u4_21, P0_21 = lane_response(lam["21"])
    u2_22, u4_22, P0_22 = lane_response(lam["22"])

    # First-order weak-axisymmetric scalar slopes.
    u2_1 = sp.simplify(sp.diff(u2_20, eps).subs(eps, 0))
    u4_1 = sp.simplify(sp.diff(u4_20, eps).subs(eps, 0))
    P1 = sp.simplify(sp.diff(P0_20, eps).subs(eps, 0))
    Xi_load = sp.simplify(N01 / N0 - D01 / D0)

    subbanner("I. Exact weak-axisymmetric grouped anomaly formulas")
    u2_an = grouped_trace_anomaly(u2_20, u2_21, u2_22)
    u4_an = grouped_trace_anomaly(u4_20, u4_21, u4_22)
    P0_an = grouped_trace_anomaly(P0_20, P0_21, P0_22)

    # We care about first-order packet tangency, so extract O(epsilon) coefficients.
    u2_bar_1 = sp.simplify(sp.diff(u2_an["bar"], eps).subs(eps, 0))
    u4_bar_1 = sp.simplify(sp.diff(u4_an["bar"], eps).subs(eps, 0))
    P0_bar_1 = sp.simplify(sp.diff(P0_an["bar"], eps).subs(eps, 0))

    u2_a_1 = sp.simplify(sp.diff(u2_an["a"], eps).subs(eps, 0))
    u2_b_1 = sp.simplify(sp.diff(u2_an["b"], eps).subs(eps, 0))
    u4_a_1 = sp.simplify(sp.diff(u4_an["a"], eps).subs(eps, 0))
    u4_b_1 = sp.simplify(sp.diff(u4_an["b"], eps).subs(eps, 0))
    P0_a_1 = sp.simplify(sp.diff(P0_an["a"], eps).subs(eps, 0))
    P0_b_1 = sp.simplify(sp.diff(P0_an["b"], eps).subs(eps, 0))

    print("d/d epsilon ubar(u2)|0 =")
    sp.pprint(u2_bar_1)
    print("d/d epsilon ubar(u4)|0 =")
    sp.pprint(u4_bar_1)
    print("d/d epsilon ubar(P0)|0 =")
    sp.pprint(P0_bar_1)

    expect_zero("first-order ubar(u2) shift", u2_bar_1)
    expect_zero("first-order ubar(u4) shift", u4_bar_1)
    expect_zero("first-order ubar(P0) shift", P0_bar_1)

    expect_zero("first-order a(u2) - u2^(1)/4", sp.simplify(u2_a_1 - u2_1 / 4))
    expect_zero("first-order b(u2) - 3 u2^(1)/4", sp.simplify(u2_b_1 - 3 * u2_1 / 4))
    expect_zero("first-order a(u4) - u4^(1)/4", sp.simplify(u4_a_1 - u4_1 / 4))
    expect_zero("first-order b(u4) - 3 u4^(1)/4", sp.simplify(u4_b_1 - 3 * u4_1 / 4))
    expect_zero("first-order a(P0) - P1/4", sp.simplify(P0_a_1 - P1 / 4))
    expect_zero("first-order b(P0) - 3 P1/4", sp.simplify(P0_b_1 - 3 * P1 / 4))

    subbanner("II. Exact true packet-tangency slopes")
    u2_1_expected = sp.simplify((-D0 * D21 + D2 * D01) / D0**2)
    u4_1_expected = sp.simplify((D0 * (-D0 * D41 + 2 * D2 * D21 - D4 * D01) + 2 * D01 * (D0 * D4 - D2**2)) / D0**3)
    P1_expected = sp.simplify((N01 * D0 - N0 * D01) / D0**2)

    expect_zero("u2^(1) formula", sp.simplify(u2_1 - u2_1_expected))
    expect_zero("u4^(1) formula", sp.simplify(u4_1 - u4_1_expected))
    expect_zero("P1 formula", sp.simplify(P1 - P1_expected))
    expect_zero("P1 - P0*Xi_load", sp.simplify(P1 - P0 * Xi_load))

    print("u2^(1) =")
    sp.pprint(sp.simplify(u2_1))
    print("u4^(1) =")
    sp.pprint(sp.simplify(u4_1))
    print("P1/P0 =")
    sp.pprint(Xi_load)

    subbanner("III. One-pole specialization")
    u2sym = sp.symbols("u2", real=True)
    one_pole_sub = {D2: -u2sym * D0, D4: D0 * (u2sym**2 - 4 * u2sym**2)}
    # D4 = D0(u2^2 - u4) with u4 = 4u2^2 -> D4 = -3u2^2 D0.
    one_pole_sub[D4] = -3 * u2sym**2 * D0

    u2_1_pole = sp.simplify(u2_1.subs(one_pole_sub))
    u4_1_pole = sp.simplify(u4_1.subs(one_pole_sub))
    print("u2^(1) on one-pole baseline =")
    sp.pprint(u2_1_pole)
    print("u4^(1) on one-pole baseline =")
    sp.pprint(u4_1_pole)

    # If u2^(1)=0, then D21 = -u2 D01, and u4^(1)=0 reduces to D41 = -3 u2^2 D01.
    D41_req = sp.simplify(sp.solve(sp.Eq(u4_1_pole.subs(D21, -u2sym * D01), 0), D41)[0])
    print("D41 required by u2^(1)=u4^(1)=0 on one-pole baseline =")
    sp.pprint(D41_req)
    expect_zero("one-pole D41 relation", sp.simplify(D41_req + 3 * u2sym**2 * D01))

    subbanner("IV. Exact canonical compensated-branch translation to (K1, H_even, Xi)")
    # Canonical compensated branch used in the notes: u2 = 1/9, u4 = 4/81.
    canonical_sub = {D2: -D0 / 9, D4: -D0 / 27}
    K1 = sp.simplify(D21 + D01 / 9)
    H_even = sp.simplify(D41 - sp.Rational(2, 3) * D21 - D01 / 27)

    u2_1_canonical = sp.simplify(u2_1.subs(canonical_sub))
    u4_1_canonical = sp.simplify(u4_1.subs(canonical_sub))

    print("u2^(1) on canonical branch =")
    sp.pprint(u2_1_canonical)
    print("u4^(1) on canonical branch =")
    sp.pprint(u4_1_canonical)
    print("K1 =")
    sp.pprint(K1)
    print("H_even =")
    sp.pprint(H_even)

    expect_zero("u2^(1) + K1/D0 on canonical branch", sp.simplify(u2_1_canonical + K1 / D0))
    expect_zero("u4^(1) + (H_even + 8 K1/9)/D0 on canonical branch", sp.simplify(u4_1_canonical + (H_even + sp.Rational(8, 9) * K1) / D0))

    banner("STAGE 215 LEDGER")
    print("1. On a genuine weak-axisymmetric grouped branch, the weighted trace of the lane data")
    print("   is unchanged at first order, so the scalar trace defects (Delta_pole, Delta_norm)")
    print("   are automatically invisible at O(epsilon).")
    print("2. The live first-order Packet-A slots are therefore exactly the anisotropy amplitudes")
    print("      u2^(1),  u4^(1),  P1.")
    print("3. Their grouped defects are universal:")
    print("      a = epsilon s/4,   b = 3 epsilon s/4")
    print("   for the grouped signature (20,21,22) ~ (1,1/2,-1).")
    print("4. True first-order packet tangency is")
    print("      u2^(1) = 0,   u4^(1) = 0,   P1 = 0.")
    print("5. On the canonical compensated branch u2 = 1/9, this is exactly equivalent to")
    print("      K1 = 0,   H_even = 0,   Xi_load = 0.")
    print("6. Off that canonical branch, the surrogate gates (K1, H_even) need not coincide")
    print("   with the true packet-tangency conditions.")
