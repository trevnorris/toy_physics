#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp

from fivepn_stage268_271_common import banner, subbanner, expect_zero

"""
Stage 268 — exact effective one-pole diagnostic for arbitrary scalar/geometry lanes.

What this script does
---------------------
1. Starts from the exact Stage-78 contamination coefficients

       a_2 = Ω_Q^2 M_0^2 G_2 / (K_pole G_0^2),
       a_4 = Ω_Q^4 M_0^2 (G_0 G_4 - G_2^2) / (K_pole G_0^3),

   which depend only on the first three low-frequency scalar-lane moments
   (G_0,G_2,G_4).

2. Shows that *any* scalar/geometry lane with nonzero G_2,G_4 can be mapped
   exactly, for contamination purposes, to one effective contact-plus-one-pole
   lane with

       Ω_eff^2 = G_2 / G_4,
       c_eff   = G_2^2 / (G_0 G_4).

   This is an identity for (a_2,a_4), not an approximation.

3. Applies the same reduction to the first concrete contact-plus-two-pole scalar
   lane and isolates the exact decomposition

       G_0 G_4 - G_2^2
       = G_c (R_1/L_1^2 + R_2/L_2^2)
         + R_1 R_2 (1/L_1 - 1/L_2)^2,

   which cleanly separates contact support from pole-separation support.

So after this stage the corridor theorems from Stage 266 can be applied to a
general scalar lane through the exact effective one-pole diagnostic
(Ω_eff^2, c_eff), without pretending the underlying lane is literally one pole.
"""


def main() -> None:
    banner("STAGE 268 — EXACT EFFECTIVE ONE-POLE DIAGNOSTIC")

    M0, Kpole, OmegaQ = sp.symbols("M_0 K_pole Omega_Q", positive=True, real=True)
    G0, G2, G4 = sp.symbols("G_0 G_2 G_4", positive=True, real=True)

    mu_mix = sp.symbols("mu_mix", positive=True, real=True)
    Omega_eff2, c_eff, r_eff = sp.symbols("Omega_eff2 c_eff r_eff", positive=True, real=True)

    a2 = sp.simplify(OmegaQ**2 * M0**2 * G2 / (Kpole * G0**2))
    a4 = sp.simplify(OmegaQ**4 * M0**2 * (G0 * G4 - G2**2) / (Kpole * G0**3))

    subbanner("I. Exact Stage-78 coefficient pair")
    print("a_2 =")
    sp.pprint(a2)
    print("a_4 =")
    sp.pprint(a4)

    subbanner("II. Exact effective one-pole reduction of any scalar lane")
    Omega_eff2_def = sp.simplify(G2 / G4)
    c_eff_def = sp.simplify(G2**2 / (G0 * G4))
    Gpole_eff = sp.simplify(G2**2 / G4)
    Gcontact_eff = sp.simplify(G0 - Gpole_eff)
    r_eff_def = sp.simplify(OmegaQ**2 / Omega_eff2_def)

    a2_eff = sp.simplify((M0**2 / (Kpole * G0)) * r_eff_def * c_eff_def)
    a4_eff = sp.simplify((M0**2 / (Kpole * G0)) * r_eff_def**2 * c_eff_def * (1 - c_eff_def))

    print("Omega_eff^2 := G_2 / G_4 =")
    sp.pprint(Omega_eff2_def)
    print("c_eff := G_2^2 / (G_0 G_4) =")
    sp.pprint(c_eff_def)
    print("G_pole,eff := G_2^2 / G_4 =")
    sp.pprint(Gpole_eff)
    print("G_contact,eff := G_0 - G_pole,eff =")
    sp.pprint(Gcontact_eff)

    expect_zero("a_2 - a_2^(eff one-pole)", sp.simplify(a2 - a2_eff))
    expect_zero("a_4 - a_4^(eff one-pole)", sp.simplify(a4 - a4_eff))
    expect_zero("G_contact,eff/G_0 - (1-c_eff)", sp.simplify(Gcontact_eff / G0 - (1 - c_eff_def)))
    expect_zero("G_pole,eff/G_0 - c_eff", sp.simplify(Gpole_eff / G0 - c_eff_def))

    subbanner("III. First exact consequence")
    print("The contamination pair (a_2,a_4) depends only on (G_0,G_2,G_4).")
    print("Therefore every scalar/geometry lane with those first three moments is exactly")
    print("equivalent, for Stage-78 contamination purposes, to one effective contact-plus-pole lane")
    print("with pole scale Omega_eff^2 = G_2/G_4 and pole fraction c_eff = G_2^2/(G_0 G_4).")

    subbanner("IV. Exact two-pole contact+pole decomposition")
    Gc, R1, R2, L1, L2 = sp.symbols("G_c R_1 R_2 L_1 L_2", positive=True, real=True)

    G0_2p = sp.simplify(Gc + R1 + R2)
    G2_2p = sp.simplify(R1 / L1 + R2 / L2)
    G4_2p = sp.simplify(R1 / L1**2 + R2 / L2**2)

    decomposition = sp.simplify(G0_2p * G4_2p - G2_2p**2)
    decomposition_expected = sp.simplify(
        Gc * (R1 / L1**2 + R2 / L2**2)
        + R1 * R2 * (1 / L1 - 1 / L2) ** 2
    )

    print("For D_g = G_c + R_1/(1-s/L_1) + R_2/(1-s/L_2),")
    print("G_0 G_4 - G_2^2 =")
    sp.pprint(decomposition)
    print("Exact decomposition =")
    sp.pprint(decomposition_expected)

    expect_zero("two-pole decomposition identity", sp.simplify(decomposition - decomposition_expected))

    c_eff_2p = sp.simplify((G2_2p**2) / (G0_2p * G4_2p))
    print("c_eff^(two-pole lane) =")
    sp.pprint(c_eff_2p)

    subbanner("V. Structural interpretation")
    print("1. If G_contact,eff = 0, the contamination-equivalent diagnostic is pure-pole (c_eff = 1).")
    print("2. If G_pole,eff = 0, the contamination-equivalent diagnostic is pure-contact (c_eff = 0).")
    print("3. For a genuine two-pole lane the quantity G_0 G_4 - G_2^2 receives")
    print("   support both from the contact slot and from pole separation.")
    print("4. So a multi-pole scalar lane can be safely compared to the Stage-266 corridor")
    print("   without collapsing it naively to a literal one-pole model first: the exact")
    print("   bridge is the moment-matched pair (Omega_eff^2, c_eff).")

    banner("STAGE 268 LEDGER")
    print("1. The Stage-78 contamination pair (a_2,a_4) depends only on (G_0,G_2,G_4).")
    print("2. Any scalar/geometry lane is therefore exactly equivalent, for contamination purposes,")
    print("   to one effective contact-plus-pole diagnostic with")
    print("      Omega_eff^2 = G_2/G_4,")
    print("      c_eff       = G_2^2/(G_0 G_4).")
    print("3. This is an identity for (a_2,a_4), not a heuristic fit.")
    print("4. For a two-pole scalar lane, the obstruction term G_0 G_4 - G_2^2 splits exactly into")
    print("   a contact contribution and a pole-separation contribution.")
    print("5. So the next honest step is to evaluate (Omega_eff^2, c_eff) on the actual Family-1")
    print("   monopole breathing channel rather than guess whether it is 'pure-pole' or 'mixed'.")


if __name__ == "__main__":
    main()
