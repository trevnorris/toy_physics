#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp
import sympy as sp

from fivepn_stage258_261_common import banner, subbanner, expect_zero

"""
Stage 260 — exact blocking monotonicity and non-twin asymmetry demand.

This script takes the Stage-35 support-demand quantity

    zeta_req = c_pole / [1 - eps_blk - (1 - 2 eps_blk) c_pole]

and proves three exact facts that matter for the post-257 continuation:

1. zeta_req is strictly increasing in c_pole on the admissible branch;
2. blocking raises zeta_req in the twin-safe regime c_pole < 1/2,
   leaves it unchanged at c_pole = 1/2, and lowers it in the non-twin regime
   c_pole > 1/2;
3. the exact non-twin asymmetry excess is

       Delta_zeta = zeta_req - 1
                  = (1 - eps_blk)(2 c_pole - 1)
                    / [1 - eps_blk - (1 - 2 eps_blk) c_pole].

So blocking hurts the symmetric-twin side but helps the non-twin side.  This is
exactly the kind of non-obvious branch behavior that gets lost if one freezes too
many coefficients too early.
"""


def main() -> None:
    banner("STAGE 260 — EXACT BLOCKING MONOTONICITY AND ASYMMETRY DEMAND")

    c_pole, eps_blk = sp.symbols("c_pole eps_blk", real=True)
    zeta_req = sp.simplify(c_pole / (1 - eps_blk - (1 - 2 * eps_blk) * c_pole))
    dz_dc = sp.factor(sp.diff(zeta_req, c_pole))
    dz_de = sp.factor(sp.diff(zeta_req, eps_blk))
    Delta_zeta = sp.factor(sp.simplify(zeta_req - 1))

    subbanner("I. Exact demand and its monotonic derivatives")
    print("zeta_req(c_pole; eps_blk) =")
    sp.pprint(zeta_req)
    print("d zeta_req / d c_pole =")
    sp.pprint(dz_dc)
    print("d zeta_req / d eps_blk =")
    sp.pprint(dz_de)
    print("Delta_zeta := zeta_req - 1 =")
    sp.pprint(Delta_zeta)

    expect_zero(
        "Delta_zeta - explicit formula",
        sp.simplify(
            Delta_zeta
            - (1 - eps_blk) * (2 * c_pole - 1) / (1 - eps_blk - (1 - 2 * eps_blk) * c_pole)
        ),
    )

    subbanner("II. Exact sign theorems")
    print("For 0 <= eps_blk < 1, the derivative with respect to c_pole is strictly positive:")
    print("  d zeta_req / d c_pole = (1 - eps_blk) / [1 - eps_blk - (1 - 2 eps_blk)c_pole]^2 > 0.")
    print()
    print("The blocking derivative is controlled by the universal sign factor (1 - 2 c_pole):")
    print("  d zeta_req / d eps_blk = c_pole (1 - 2 c_pole) / [1 - eps_blk - (1 - 2 eps_blk)c_pole]^2.")
    print()
    print("So the exact branch classification is:")
    print("  c_pole < 1/2   -> blocking increases zeta_req (hurts the twin-safe side),")
    print("  c_pole = 1/2   -> blocking leaves zeta_req unchanged,")
    print("  c_pole > 1/2   -> blocking decreases zeta_req (helps the non-twin side).")

    subbanner("III. Exact checkpoints")
    zeta_half = sp.simplify(zeta_req.subs(c_pole, sp.Rational(1, 2)))
    zeta_actual = sp.simplify(zeta_req.subs(c_pole, sp.Rational(1, 4)))
    print("zeta_req(c_pole = 1/2) =")
    sp.pprint(zeta_half)
    print("zeta_req(c_pole = 1/4) =")
    sp.pprint(zeta_actual)
    expect_zero("zeta_req(c_pole=1/2)-1", sp.simplify(zeta_half - 1))
    expect_zero("actual isotropic branch", sp.simplify(zeta_actual - 1 / (3 - 2 * eps_blk)))

    subbanner("IV. Numerical illustration of the sign change")
    mp.mp.dps = 60
    eps_num = mp.mpf("0.2")
    for c_num in [mp.mpf("0.25"), mp.mpf("0.5"), mp.mpf("0.6666666666666667")]:
        z = c_num / (1 - eps_num - (1 - 2 * eps_num) * c_num)
        dz = c_num * (1 - 2 * c_num) / (1 - eps_num - (1 - 2 * eps_num) * c_num) ** 2
        print(f"c_pole = {c_num}: zeta_req ≈ {z}, d zeta_req / d eps_blk ≈ {dz}")

    banner("STAGE 260 LEDGER")
    print("1. The required coherent support ratio zeta_req is strictly increasing in c_pole on the")
    print("   admissible branch, so the selected-branch demand ordering in c_pole language is exact.")
    print("2. Blocking changes zeta_req with the sign of (1 - 2 c_pole):")
    print("      c_pole < 1/2  -> blocking hurts the twin-safe regime,")
    print("      c_pole > 1/2  -> blocking helps the non-twin regime.")
    print("3. The exact excess beyond the symmetric-twin threshold is")
    print("      Delta_zeta = (1-eps_blk)(2 c_pole - 1) / [1-eps_blk-(1-2 eps_blk)c_pole].")
    print("   So the actual isotropic c_pole = 1/4 branch remains twin-safe for every admissible")
    print("   blocking value, while non-twin branches become easier to support as blocking grows.")


if __name__ == "__main__":
    main()
