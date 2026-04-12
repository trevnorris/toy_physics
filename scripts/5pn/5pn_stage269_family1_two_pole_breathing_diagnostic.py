#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp

from fivepn_stage268_271_common import banner, subbanner, expect_zero

"""
Stage 269 — the actual 2PN Family-1 monopole breathing channel as an exact
two-pole scalar-lane diagnostic for Stage-78 contamination.

This script uses the explicit Family-1 scalar monopole channel recorded in the
2PN constructive appendix,

    K_00(s) = -757/2520 + R_-/(1-s/lambda_-) + R_+/(1-s/lambda_+),

with
    K_00(0) = 4/45,
and the quoted pole / residue values.

It then computes the first three low-frequency moments (G_0,G_2,G_4) and turns
them into the exact contamination-equivalent one-pole diagnostic

    lambda_coeff = G_2 / G_4,
    c_eff        = G_2^2 / (G_0 G_4).

This is the right one-pole reduction for Stage-78 contamination because it
preserves (a_2,a_4) exactly.

A crucial bookkeeping warning is also recorded:
the coefficient-equivalent one-pole scale `lambda_coeff` is not the same object
as the interval-accuracy Padé pole `lambda_eff` quoted in the 2PN note.
They answer different questions and should not be conflated.
"""


def main() -> None:
    banner("STAGE 269 — FAMILY-1 TWO-POLE BREATHING DIAGNOSTIC")

    # Exact scalar static value fixed in the 2PN note.
    G0 = sp.Rational(4, 45)

    # Exact static contact term from the same Family-1 channel.
    Gc = -sp.Rational(757, 2520)

    # Pole / residue data quoted in the constructive 2PN appendix.
    lam_minus = sp.Float("6.405572392138922", 50)
    lam_plus = sp.Float("254.444968136936126", 50)
    R_minus = sp.Float("0.002552474771738", 50)
    R_plus = sp.Float("0.386733239513976", 50)

    # The note also quotes a Padé-on-interval effective pole.
    lam_pade = sp.Float("202.923516367519028", 50)

    G2 = sp.simplify(R_minus / lam_minus + R_plus / lam_plus)
    G4 = sp.simplify(R_minus / lam_minus**2 + R_plus / lam_plus**2)

    lam_coeff = sp.simplify(G2 / G4)
    c_eff = sp.simplify(G2**2 / (G0 * G4))
    Gpole_eff = sp.simplify(G0 * c_eff)
    Gcontact_eff = sp.simplify(G0 * (1 - c_eff))

    subbanner("I. Actual Family-1 two-pole scalar lane")
    print("K_00(s) = -757/2520 + R_-/(1-s/lambda_-) + R_+/(1-s/lambda_+)")
    print("with")
    print("  lambda_- =", lam_minus)
    print("  lambda_+ =", lam_plus)
    print("  R_-      =", R_minus)
    print("  R_+      =", R_plus)
    print("and exact static value K_00(0) =")
    sp.pprint(G0)

    subbanner("II. First three low-frequency moments")
    print("G_0 =")
    sp.pprint(G0)
    print("G_2 =")
    sp.pprint(G2)
    print("G_4 =")
    sp.pprint(G4)

    expect_zero("G_0 - 4/45", sp.simplify(G0 - sp.Rational(4, 45)))
    # Because the pole/residue inputs are decimal values from the note, do not over-claim exact
    # recovery of G_0 from G_c + R_- + R_+. Just display the small mismatch.
    static_mismatch = sp.N(Gc + R_minus + R_plus - G0, 30)
    print("Gc + R_- + R_+ - G_0  (numerical residue truncation check) =")
    sp.pprint(static_mismatch)

    subbanner("III. Exact contamination-equivalent one-pole diagnostic")
    print("lambda_coeff := G_2 / G_4 =")
    sp.pprint(sp.N(lam_coeff, 30))
    print("c_eff := G_2^2 / (G_0 G_4) =")
    sp.pprint(sp.N(c_eff, 30))
    print("G_pole,eff = G_0 c_eff =")
    sp.pprint(sp.N(Gpole_eff, 30))
    print("G_contact,eff = G_0 (1-c_eff) =")
    sp.pprint(sp.N(Gcontact_eff, 30))

    print()
    print("Padé-on-interval one-pole scale quoted in the 2PN note:")
    print("lambda_pade =")
    sp.pprint(lam_pade)
    print("lambda_pade / lambda_coeff =")
    sp.pprint(sp.N(lam_pade / lam_coeff, 30))

    subbanner("IV. Interpretation")
    print("1. The actual Family-1 scalar lane is explicitly a coupled two-pole monopole breathing channel.")
    print("2. For Stage-78 contamination, however, only (G_0,G_2,G_4) matter, so the correct exact")
    print("   one-pole diagnostic is the coefficient-equivalent pair (lambda_coeff, c_eff).")
    print("3. On this branch the diagnostic pole fraction is")
    print("      c_eff ~= ", sp.N(c_eff, 16))
    print("   so the breathing lane is neither pure-pole nor pure-contact at the contamination level.")
    print("4. The positive effective contact term G_contact,eff comes partly from the explicit contact")
    print("   slot and partly from pole-separation support; it is not identical to the raw static contact")
    print("   coefficient G_c.")
    print("5. The interval-accuracy Padé pole lambda_pade is a different object. It is useful for fitting")
    print("   K_00(s) on a finite s-window, but it is not the exact contamination-equivalent pole.")

    banner("STAGE 269 LEDGER")
    print("1. The actual constructive scalar lane from the 2PN Family-1 module is a coupled two-pole")
    print("   monopole breathing channel K_00(s), not a literal one-pole lane.")
    print("2. For Stage-78 contamination coefficients, its exact one-pole diagnostic is")
    print("      lambda_coeff = G_2/G_4,")
    print("      c_eff        = G_2^2/(G_0 G_4).")
    print("3. Numerically on the quoted Family-1 branch:")
    print("      lambda_coeff ~= ", sp.N(lam_coeff, 16))
    print("      c_eff        ~= ", sp.N(c_eff, 16))
    print("4. So the actual Family-1 breathing lane behaves, for contamination purposes, like a mixed")
    print("   contact-plus-pole scalar lane with pole fraction around 0.61.")
    print("5. The quoted Padé pole lambda_pade ~= ", sp.N(lam_pade, 16), " is not the same object and")
    print("   should not be used as the Stage-78 contamination pole without an explicit additional argument.")


if __name__ == "__main__":
    main()
