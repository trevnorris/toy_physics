#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp
import sympy as sp

from fivepn_stage258_261_common import banner, subbanner, expect_zero, family1_refreshed_numbers

"""
Stage 258 — exact selected-branch regime split in generic isotropic one-pole
conservative language.

This script continues directly from the exact Lambda_EM refresh without
re-assuming the canonical c_pole = 1/4 split until the end.  It proves that the
Stage-35 mixed/twin/non-twin regime split becomes exceptionally sharp in the
one-pole carrier language:

    Pi_tr / C_mix = rho_alpha = 1 / (1 - c_pole),

so the three support regimes are determined only by the sign and size of the
pole weight c_pole.  This is the first exact branch-level statement that the
mixed-only region is impossible on the physical positive-pole branch, while the
actual isotropic c_pole = 1/4 branch sits comfortably inside the universal
lowest-twin-safe strip.
"""


def main() -> None:
    banner("STAGE 258 — EXACT ONE-POLE SELECTED-BRANCH REGIME SPLIT")

    c_pole, eps_blk = sp.symbols("c_pole eps_blk", real=True)
    rho_alpha = sp.simplify(1 / (1 - c_pole))
    zeta_req = sp.simplify(c_pole / (1 - eps_blk - (1 - 2 * eps_blk) * c_pole))
    Pi_over_Cmix = sp.simplify((1 + (1 - 2 * eps_blk) * zeta_req) / (1 - eps_blk * zeta_req))

    subbanner("I. Exact one-pole map into the selected-branch demand variables")
    print("rho_alpha(c_pole) =")
    sp.pprint(rho_alpha)
    print("zeta_req(c_pole; eps_blk) =")
    sp.pprint(sp.factor(zeta_req))
    print("Pi_tr / C_mix =")
    sp.pprint(Pi_over_Cmix)
    expect_zero("Pi/C_mix - rho_alpha", sp.simplify(Pi_over_Cmix - rho_alpha))

    subbanner("II. Exact support-regime boundaries in c_pole language")
    c_mixed = sp.solve(sp.Eq(rho_alpha, 1), c_pole)[0]
    c_twin = sp.solve(sp.Eq(rho_alpha, 2), c_pole)[0]
    print("Boundary rho_alpha = 1  ->  c_pole =")
    sp.pprint(c_mixed)
    print("Boundary rho_alpha = 2  ->  c_pole =")
    sp.pprint(c_twin)
    expect_zero("mixed-only boundary", sp.simplify(c_mixed))
    expect_zero("lowest-twin boundary", sp.simplify(c_twin - sp.Rational(1, 2)))
    print()
    print("So the exact regime split is:")
    print("  c_pole <= 0         : mixed-only enough,")
    print("  0 < c_pole <= 1/2   : symmetric lowest twin enough,")
    print("  c_pole > 1/2        : non-twin asymmetry required.")

    subbanner("III. Exact translation into zeta_req language")
    zeta_at_zero = sp.simplify(zeta_req.subs(c_pole, 0))
    zeta_at_half = sp.simplify(zeta_req.subs(c_pole, sp.Rational(1, 2)))
    print("zeta_req(c_pole = 0) =")
    sp.pprint(zeta_at_zero)
    print("zeta_req(c_pole = 1/2) =")
    sp.pprint(zeta_at_half)
    expect_zero("zeta_req(c_pole=0)", zeta_at_zero)
    expect_zero("zeta_req(c_pole=1/2)-1", sp.simplify(zeta_at_half - 1))
    print()
    print("So the same exact regime split is:")
    print("  zeta_req <= 0   : mixed-only enough,")
    print("  0 < zeta_req <= 1 : symmetric lowest twin enough,")
    print("  zeta_req > 1    : non-twin asymmetry required.")

    subbanner("IV. Physical positive-pole branch and the actual isotropic point")
    nums = family1_refreshed_numbers()
    cmax0 = nums["c_pole_max_0"]
    ccan = mp.mpf("1") / 4
    rho_can = 1 / (1 - ccan)
    zeta_can = ccan / (1 - ccan)
    print(f"Family-1 hard ceiling at eps_blk = 0: c_pole < {cmax0}")
    print(f"Actual isotropic point: c_pole = {ccan}, rho_alpha = {rho_can}, zeta_req = {zeta_can}")
    print(f"Distance to the universal lowest-twin boundary c_pole = 1/2: {mp.mpf('1')/2 - ccan}")
    print(f"Residual non-twin Family-1 corridor at eps_blk = 0: {cmax0 - mp.mpf('1')/2}")
    print()
    print("So on the physical positive-pole branch:")
    print("  - mixed-only sufficiency is impossible unless c_pole <= 0,")
    print("  - the actual isotropic c_pole = 1/4 branch is exactly in the lowest-twin-safe regime,")
    print("  - and the explicit Family-1 side still leaves a finite non-twin corridor above c_pole = 1/2.")

    banner("STAGE 258 LEDGER")
    print("1. The selected-branch demand variable Pi_tr / C_mix is exactly rho_alpha = 1/(1-c_pole),")
    print("   independent of blocking.")
    print("2. Therefore the support-regime split is exact in one-pole language:")
    print("      c_pole <= 0         : mixed-only enough,")
    print("      0 < c_pole <= 1/2   : symmetric lowest twin enough,")
    print("      c_pole > 1/2        : non-twin asymmetry required.")
    print("3. So on every physical positive-pole branch the mixed-only regime is already excluded,")
    print("   while the actual isotropic c_pole = 1/4 branch stays safely inside the universal")
    print("   lowest-twin strip.")


if __name__ == "__main__":
    main()
