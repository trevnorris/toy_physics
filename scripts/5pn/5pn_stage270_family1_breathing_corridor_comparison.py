#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp

from fivepn_stage268_271_common import banner, subbanner, expect_zero, family1_refreshed_numbers

"""
Stage 270 — exact corridor comparison for the actual Family-1 two-pole breathing lane.

This script takes the exact contamination-equivalent one-pole diagnostic of
Stage 269 and plugs it into the Stage-266 corridor thresholds.

Because the actual constructive scalar lane is two-pole, the comparison is
performed in the *exact* coefficient-equivalent sense: the resulting thresholds
are the exact thresholds for the contamination pair (a_2,a_4), not a heuristic
one-pole guess.

The outcome is a single danger variable

    u_eff = r_eff (1 - c_eff),

where `r_eff = Omega_Q^2 / lambda_coeff` is the grouped-quadrupole / scalar-lane
pole ratio in the same reduced spectral variable.

Then:

- support demand grows initially iff      r_eff > 2 / (1 - c_eff),
- the twin margin shrinks initially iff   r_eff > 4 / (1 - c_eff),
- actual twin failure needs               r_eff >= (4+2*sqrt(2)) / (1 - c_eff),
- actual Family-1 failure needs           r_eff >= u_F1^fail / (1 - c_eff).

So the only remaining scalar-lane uncertainty is the actual pole ratio
`r_eff`; the scalar-lane composition is already fixed by the concrete Family-1
breathing channel.
"""


def main() -> None:
    banner("STAGE 270 — FAMILY-1 BREATHING LANE VS THE CORRIDOR THRESHOLDS")

    # Rebuild the Stage-269 coefficient-equivalent diagnostic.
    G0 = sp.Rational(4, 45)
    lam_minus = sp.Float("6.405572392138922", 80)
    lam_plus = sp.Float("254.444968136936126", 80)
    R_minus = sp.Float("0.002552474771738", 80)
    R_plus = sp.Float("0.386733239513976", 80)

    G2 = sp.simplify(R_minus / lam_minus + R_plus / lam_plus)
    G4 = sp.simplify(R_minus / lam_minus**2 + R_plus / lam_plus**2)

    lam_coeff = sp.simplify(G2 / G4)
    c_eff = sp.simplify(G2**2 / (G0 * G4))
    one_minus_c = sp.simplify(1 - c_eff)

    r_eff = sp.symbols("r_eff", positive=True, real=True)
    u_eff = sp.simplify(r_eff * one_minus_c)

    subbanner("I. Exact coefficient-equivalent danger variable")
    print("lambda_coeff =")
    sp.pprint(sp.N(lam_coeff, 30))
    print("c_eff =")
    sp.pprint(sp.N(c_eff, 30))
    print("1 - c_eff =")
    sp.pprint(sp.N(one_minus_c, 30))
    print("u_eff := r_eff (1 - c_eff) =")
    sp.pprint(u_eff)

    nums = family1_refreshed_numbers()
    cstar = sp.Float(str(nums["c_pole_max_0"]), 80)

    u_demand = sp.Integer(2)
    u_twin_shrink = sp.Integer(4)
    u_twin_fail = sp.simplify(4 + 2 * sp.sqrt(2))
    u_F1_shrink = sp.simplify(8 * cstar)
    u_F1_fail = sp.simplify(8 * cstar + 4 * sp.sqrt(cstar * (4 * cstar - 1)))

    r_demand = sp.simplify(u_demand / one_minus_c)
    r_twin_shrink = sp.simplify(u_twin_shrink / one_minus_c)
    r_twin_fail = sp.simplify(u_twin_fail / one_minus_c)
    r_F1_shrink = sp.simplify(u_F1_shrink / one_minus_c)
    r_F1_fail = sp.simplify(u_F1_fail / one_minus_c)

    subbanner("II. Exact threshold ladder for the actual Family-1 breathing lane")
    print("support demand grows initially iff      r_eff >")
    sp.pprint(sp.N(r_demand, 20))
    print("twin margin shrinks initially iff       r_eff >")
    sp.pprint(sp.N(r_twin_shrink, 20))
    print("actual twin failure requires            r_eff >=")
    sp.pprint(sp.N(r_twin_fail, 20))
    print("Family-1 margin shrinks initially iff   r_eff >")
    sp.pprint(sp.N(r_F1_shrink, 20))
    print("actual Family-1 failure requires        r_eff >=")
    sp.pprint(sp.N(r_F1_fail, 20))

    expect_zero(
        "u_eff(twin failure threshold) - (4 + 2 sqrt(2))",
        sp.simplify(u_eff.subs(r_eff, r_twin_fail) - u_twin_fail),
    )
    expect_zero(
        "u_eff(Family-1 failure threshold) - (8 c_* + 4 sqrt(c_*(4 c_* - 1)))",
        sp.simplify(u_eff.subs(r_eff, r_F1_fail) - u_F1_fail),
    )

    subbanner("III. Interpretation")
    print("1. The explicit Family-1 breathing lane fixes the contamination-equivalent pole fraction to")
    print("      c_eff ~= ", sp.N(c_eff, 16))
    print("   so the only remaining scalar-lane ambiguity is the actual pole ratio r_eff = Omega_Q^2/lambda_coeff.")
    print("2. Numerically this pushes the exact danger thresholds to")
    print("      r_eff > ", sp.N(r_demand, 10), "   (support demand starts growing),")
    print("      r_eff > ", sp.N(r_twin_shrink, 10), "   (twin margin starts shrinking),")
    print("      r_eff >= ", sp.N(r_twin_fail, 10), " (actual twin failure possible),")
    print("      r_eff >= ", sp.N(r_F1_fail, 10), " (actual Family-1 failure possible).")
    print("3. So the first known constructive scalar lane is not dangerous by composition alone.")
    print("   It becomes dangerous only if the grouped quadrupole pole is much faster than the")
    print("   coefficient-equivalent breathing-lane pole in the same reduced spectral variable.")

    banner("STAGE 270 LEDGER")
    print("1. The actual Family-1 two-pole monopole breathing lane fixes the exact coefficient-equivalent")
    print("   contact/pole split seen by Stage-78 contamination.")
    print("2. On that branch the scalar-lane danger variable is")
    print("      u_eff = r_eff (1 - c_eff),")
    print("   with c_eff ~= ", sp.N(c_eff, 16), ".")
    print("3. Therefore the exact Stage-266 thresholds translate into the explicit pole-ratio ladder")
    print("      r_eff > ", sp.N(r_demand, 12), "   -> support demand grows initially,")
    print("      r_eff > ", sp.N(r_twin_shrink, 12), "   -> twin margin shrinks,")
    print("      r_eff >= ", sp.N(r_twin_fail, 12), " -> actual twin failure can occur,")
    print("      r_eff >= ", sp.N(r_F1_fail, 12), " -> actual Family-1 failure can occur.")
    print("4. So the known constructive scalar lane is compatible with the corridor unless the physical")
    print("   grouped-quadrupole pole is very fast relative to the contamination-equivalent breathing pole.")


if __name__ == "__main__":
    main()
