#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp

from fivepn_stage258_261_common import banner, subbanner, expect_zero, family1_refreshed_numbers

"""
Stage 264 — exact second-order weak-anisotropy tolerance theorem on the actual
selected grouped-P2 branch.

Stage 78 in the notes says that the first nonzero geometry contamination enters
only at O(chi^2).  This script turns that qualitative statement into exact
algebra.

Take
    epsilon_2 = a2 * chi^2,
    epsilon_4 = a4 * chi^2.

Then the exact twin-safe and Family-1-admissible conditions become explicit
quadratics in y := chi^2.  This yields three useful theorem-level facts:

1. every weak-anisotropy branch has a finite safe neighborhood around y=0,
2. the initial drift of the actual branch is governed only by a4 - 2 a2,
3. exit from the twin-safe or Family-1 strip can occur only through an exact
   quadratic root structure, not arbitrarily close to the isotropic point.
"""


def main() -> None:
    banner("STAGE 264 — EXACT SECOND-ORDER WEAK-ANISOTROPY TOLERANCE")

    a2, a4, y, eps_blk, zmax = sp.symbols("a_2 a_4 y eps_blk zeta_max", real=True)
    eps2 = a2 * y
    eps4 = a4 * y

    c_pole = sp.simplify((1 + eps4) / (4 * (1 + eps2) ** 2))
    zeta_req = sp.simplify(c_pole / (1 - eps_blk - (1 - 2 * eps_blk) * c_pole))
    cmax = sp.simplify(zmax * (1 - eps_blk) / (1 + (1 - 2 * eps_blk) * zmax))

    M_twin = sp.expand(2 * (1 + eps2) ** 2 - (1 + eps4))
    M_F1 = sp.expand(4 * cmax * (1 + eps2) ** 2 - (1 + eps4))

    subbanner("I. Exact O(chi^2) deformation of the actual branch")
    print("epsilon_2 = a_2 * chi^2,   epsilon_4 = a_4 * chi^2,   y := chi^2")
    print("c_pole(y) =")
    sp.pprint(sp.factor(c_pole))
    print("zeta_req(y) =")
    sp.pprint(sp.factor(zeta_req))
    print("M_twin(y) =")
    sp.pprint(M_twin)
    print("M_F1(y) =")
    sp.pprint(sp.factor(M_F1))

    expect_zero("M_twin(0)-1", sp.simplify(M_twin.subs(y, 0) - 1))
    expect_zero("M_F1(0)-(4 cmax - 1)", sp.simplify(M_F1.subs(y, 0) - (4 * cmax - 1)))

    subbanner("II. Exact initial drift at the isotropic point")
    dc0 = sp.simplify(sp.diff(c_pole, y).subs(y, 0))
    dz0 = sp.simplify(sp.diff(zeta_req, y).subs(y, 0))
    print("d c_pole / d(chi^2) at chi=0 =")
    sp.pprint(dc0)
    print("d zeta_req / d(chi^2) at chi=0 =")
    sp.pprint(sp.factor(dz0))

    expect_zero(
        "dc/dy|0 - (a4-2 a2)/4",
        sp.simplify(dc0 - (a4 - 2 * a2) / 4),
    )
    expect_zero(
        "dz/dy|0 - 4(1-eps_blk)(a4-2 a2)/(3-2 eps_blk)^2",
        sp.simplify(dz0 - 4 * (1 - eps_blk) * (a4 - 2 * a2) / (3 - 2 * eps_blk) ** 2),
    )
    print()
    print("So the weak-anisotropy branch initially moves")
    print("  toward larger pole weight / support demand iff  a_4 - 2 a_2 > 0,")
    print("  toward smaller pole weight / support demand iff a_4 - 2 a_2 < 0.")

    subbanner("III. Exact twin-safe quadratic and its roots")
    Delta_twin = sp.factor(sp.discriminant(M_twin, y))
    print("Discriminant Delta_twin =")
    sp.pprint(Delta_twin)

    y_twin_roots = sp.solve(sp.Eq(M_twin, 0), y)
    print("Twin-boundary roots in y = chi^2 are")
    for root in y_twin_roots:
        sp.pprint(sp.factor(root))

    if len(y_twin_roots) == 2:
        y_twin_minus, y_twin_plus = y_twin_roots
        expect_zero("M_twin(y_twin_minus)", sp.simplify(M_twin.subs(y, y_twin_minus)))
        expect_zero("M_twin(y_twin_plus)", sp.simplify(M_twin.subs(y, y_twin_plus)))

    print()
    print("Therefore:")
    print("  - if Delta_twin <= 0 (or the positive roots are absent), the O(chi^2) branch stays")
    print("    twin-safe for all y >= 0;")
    print("  - if Delta_twin > 0 and positive roots exist, the first loss of twin-safety occurs at")
    print("    the smaller positive root in y = chi^2.")

    subbanner("IV. Exact Family-1 admissibility quadratic and its roots")
    Delta_F1 = sp.factor(sp.discriminant(M_F1, y))
    print("Discriminant Delta_F1 =")
    sp.pprint(Delta_F1)

    y_F1_roots = sp.solve(sp.Eq(M_F1, 0), y)
    print("Family-1-boundary roots in y = chi^2 are")
    for root in y_F1_roots:
        sp.pprint(sp.factor(root))

    if len(y_F1_roots) == 2:
        y_F1_minus, y_F1_plus = y_F1_roots
        expect_zero("M_F1(y_F1_minus)", sp.simplify(M_F1.subs(y, y_F1_minus)))
        expect_zero("M_F1(y_F1_plus)", sp.simplify(M_F1.subs(y, y_F1_plus)))

    print()
    print("Therefore the exact blocked Family-1 corridor can fail only through the corresponding")
    print("quadratic root structure in y = chi^2.  In particular, because M_F1(0) = 4 cmax - 1 > 0")
    print("on every admissible Lambda_EM-refreshed branch, there is always a finite safe neighborhood")
    print("around the isotropic point before any Family-1 failure can occur.")

    subbanner("V. Linear special cases")
    print("If a_2 = 0, the twin and Family-1 tests collapse to exact linear thresholds:")
    print("  M_twin = 1 - a_4 y,")
    print("  M_F1   = (4 cmax - 1) - a_4 y.")
    print("So then")
    print("  y_twin,crit = 1/a_4           (if a_4 > 0),")
    print("  y_F1,crit   = (4 cmax - 1)/a_4 (if a_4 > 0).")

    subbanner("VI. Exact refreshed unblocked Family-1 margin")
    nums = family1_refreshed_numbers()
    zmax_num = sp.nsimplify(str(nums["zeta_max"]))
    cmax0_num = sp.nsimplify(str(nums["c_pole_max_0"]))
    print("On the exact Lambda_EM-refreshed unblocked branch:")
    print("zeta_max^(F1) ≈")
    sp.pprint(zmax_num)
    print("c_pole,max^(F1)(0) ≈")
    sp.pprint(cmax0_num)
    print("So M_F1(0) = 4 cmax^(F1)(0) - 1 is strictly positive numerically as well.")

    banner("STAGE 264 LEDGER")
    print("1. Stage-78 weak anisotropy enters the actual grouped-P2 branch only through")
    print("      epsilon_2 = a_2 chi^2,   epsilon_4 = a_4 chi^2,")
    print("   so the exact twin-safe and Family-1 corridor tests become explicit quadratics in y = chi^2.")
    print("2. Every such branch has a finite safe neighborhood around y = 0 because")
    print("      M_twin(0) = 1,   M_F1(0) = 4 cmax^(F1) - 1 > 0.")
    print("3. The initial drift is governed only by a_4 - 2 a_2:")
    print("      a_4 - 2 a_2 > 0  ->  the branch initially moves toward larger support demand,")
    print("      a_4 - 2 a_2 < 0  ->  the branch initially moves toward smaller support demand.")
    print("4. Any eventual exit from the twin-safe strip or from the exact Family-1 corridor can occur")
    print("   only through the corresponding exact quadratic roots in y = chi^2, not arbitrarily close")
    print("   to the actual isotropic point.")


if __name__ == "__main__":
    main()
