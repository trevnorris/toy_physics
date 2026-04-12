#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 275 — exact no-shift family and first-order outgoing-defect sensitivities.

What this script does
---------------------
1. Starts from the exact isotropic deformation law

       chi_Q = 3 (S beta^5 + 9 Sigma5) / (3 S - Sigma0).

2. Solves exactly for the full no-shift family chi_Q = 1.
3. Linearizes the outgoing defect Delta_Q = chi_Q - 1 around the canonical compact
   branch and shows

       Delta_Q = 5 delta_beta + Sigma0/3 + 9 Sigma5 + O(2),

   while the overall amplitude deformation in S drops out to first order.
4. Records the exact implication that the canonical outgoing branch is not unique:
   there is a two-parameter isotropic deformation family that leaves chi_Q = 1
   unchanged once the even moments are preserved.
"""


def main() -> None:
    banner("STAGE 275 — OUTGOING DEFECT LINEARIZATION")

    S, beta = sp.symbols("S beta", positive=True, real=True)
    Sigma0, Sigma5 = sp.symbols("Sigma_0 Sigma_5", real=True)
    chiQ = sp.simplify(3 * (S * beta**5 + 9 * Sigma5) / (3 * S - Sigma0))

    subbanner("I. Exact no-shift family chi_Q = 1")
    Sigma0_noshift = sp.solve(sp.Eq(chiQ, 1), Sigma0)[0]
    print("Exact chi_Q = 1 surface solved for Sigma0 =")
    sp.pprint(sp.factor(Sigma0_noshift))

    expect_zero(
        "chi_Q(no-shift family) - 1",
        sp.simplify(chiQ.subs(Sigma0, Sigma0_noshift) - 1),
    )

    subbanner("II. First-order outgoing-defect sensitivity")
    s, db, s0, s5, eps = sp.symbols("s d_beta s0 s5 eps", real=True)
    chi_lin = sp.series(
        chiQ.subs({S: 1 + eps * s, beta: 1 + eps * db, Sigma0: eps * s0, Sigma5: eps * s5}),
        eps,
        0,
        2,
    ).removeO()
    Delta_lin = sp.expand(chi_lin - 1)
    Delta1 = sp.expand(sp.collect(Delta_lin, eps).coeff(eps, 1))

    print("Delta_Q^(1) =")
    sp.pprint(Delta1)
    expect_zero(
        "Delta_Q^(1) - (5 d_beta + Sigma0/3 + 9 Sigma5)",
        sp.simplify(Delta1 - (5 * db + s0 / 3 + 9 * s5)),
    )

    print("So the overall amplitude deformation in S cancels out at first order.")

    subbanner("III. Useful special slice: conservative normalization already fixed")
    NQ = sp.simplify(3 / (3 * S - Sigma0))
    Sigma0_unitN = sp.solve(sp.Eq(NQ, 1), Sigma0)[0]
    chi_out_only = sp.simplify(chiQ.subs(Sigma0, Sigma0_unitN))

    print("If the conservative branch is already normalized (N_Q = 1), then Sigma0 =")
    sp.pprint(sp.factor(Sigma0_unitN))
    print("and the remaining outgoing factor reduces to")
    sp.pprint(sp.factor(chi_out_only))

    subbanner("IV. Best current reading")
    print("1. The exact isotropic no-shift family is")
    print("      Sigma0 = 3 S (1 - beta^5) - 27 Sigma5.")
    print("2. Therefore the canonical outgoing branch is not unique once the even moments are fixed;")
    print("   there is an exact two-parameter isotropic deformation family that preserves chi_Q = 1.")
    print("3. Near the canonical branch, the first-order outgoing defect is")
    print("      Delta_Q = 5 delta_beta + Sigma0/3 + 9 Sigma5 + O(2),")
    print("   so the overall amplitude scaling S is not load-bearing at linear order.")
    print("4. On the slice where the conservative normalization is already fixed (N_Q = 1),")
    print("   the remaining outgoing defect is carried only by the retarded DtN deformation data.")


if __name__ == "__main__":
    main()
