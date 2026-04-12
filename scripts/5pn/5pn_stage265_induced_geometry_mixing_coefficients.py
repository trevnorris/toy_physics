#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp

from fivepn_stage265_267_common import banner, subbanner, expect_zero

"""
Stage 265 — exact coefficient extraction for the first l=0 <-> l=2 induced
geometry-mixing mechanism.

This script starts from the exact Stage-78 mixed scalar-geometry / grouped-P2
Schur complement and turns it into the explicit weak-anisotropy coefficients

    epsilon_2 = a_2 chi^2,
    epsilon_4 = a_4 chi^2.

It then specializes to the first concrete scalar-lane model in which the scalar
geometry kernel itself is one effective contact-plus-pole branch.  On that
branch the coefficients collapse to a compact three-parameter form controlled by

    mu_mix := M_0^2 / (K_pole G_0),
    r      := Omega_Q^2 / Omega_g^2,
    c_g    := G_pole / G_0.

The exact outputs are

    a_2 = mu_mix * r * c_g,
    a_4 = mu_mix * r^2 * c_g * (1 - c_g),

and therefore

    a_4 - 2 a_2 = mu_mix * r * c_g * ( r(1-c_g) - 2 ).

So the first induced l=0 <-> l=2 geometry contamination already splits into one
overall amplitude, one pole-ratio, and one scalar-lane pole fraction.
"""


def main() -> None:
    banner("STAGE 265 — EXACT INDUCED GEOMETRY-MIXING COEFFICIENTS")

    chi, M0, Kpole = sp.symbols("chi M_0 K_pole", positive=True, real=True)
    G0, G2, G4 = sp.symbols("G_0 G_2 G_4", real=True)
    OmegaQ, w = sp.symbols("Omega_Q omega", positive=True, real=True)

    Dg = G0 + G2 * w**2 + G4 * w**4
    schur = sp.series(-chi**2 * M0**2 / Dg, w, 0, 6).removeO()

    Kg2_eff = sp.simplify(sp.expand(schur).coeff(w, 2))
    Kg4_eff = sp.simplify(sp.expand(schur).coeff(w, 4))
    eps2 = sp.simplify(OmegaQ**2 * Kg2_eff / Kpole)
    eps4 = sp.simplify(OmegaQ**4 * Kg4_eff / Kpole)

    subbanner("I. Exact Stage-78 Schur-complement contamination coefficients")
    print("D_g(omega) =")
    sp.pprint(Dg)
    print("-chi^2 M_0^2 / D_g(omega) through O(omega^4) =")
    sp.pprint(sp.expand(schur))
    print("K_(g,2)^eff =")
    sp.pprint(Kg2_eff)
    print("K_(g,4)^eff =")
    sp.pprint(Kg4_eff)
    print("epsilon_2 =")
    sp.pprint(eps2)
    print("epsilon_4 =")
    sp.pprint(eps4)

    expect_zero(
        "K_(g,2)^eff - chi^2 M_0^2 G_2 / G_0^2",
        sp.simplify(Kg2_eff - chi**2 * M0**2 * G2 / G0**2),
    )
    expect_zero(
        "K_(g,4)^eff - chi^2 M_0^2 (G_0 G_4 - G_2^2) / G_0^3",
        sp.simplify(Kg4_eff - chi**2 * M0**2 * (G0 * G4 - G2**2) / G0**3),
    )
    expect_zero(
        "epsilon_2 - chi^2 Omega_Q^2 M_0^2 G_2 / (K_pole G_0^2)",
        sp.simplify(eps2 - chi**2 * OmegaQ**2 * M0**2 * G2 / (Kpole * G0**2)),
    )
    expect_zero(
        "epsilon_4 - chi^2 Omega_Q^4 M_0^2 (G_0 G_4 - G_2^2) / (K_pole G_0^3)",
        sp.simplify(
            eps4
            - chi**2 * OmegaQ**4 * M0**2 * (G0 * G4 - G2**2) / (Kpole * G0**3)
        ),
    )

    subbanner("II. First concrete scalar-lane model: one contact plus one scalar pole")
    Gc, Gp, Omegag = sp.symbols("G_c G_p Omega_g", positive=True, real=True)
    Dg_pole = sp.simplify(Gc + Gp / (1 - w**2 / Omegag**2))
    Dg_pole_series = sp.series(Dg_pole, w, 0, 6).removeO()
    G0_pole = sp.simplify(Dg_pole_series.coeff(w, 0))
    G2_pole = sp.simplify(Dg_pole_series.coeff(w, 2))
    G4_pole = sp.simplify(Dg_pole_series.coeff(w, 4))

    print("D_g^(contact+pole)(omega) =")
    sp.pprint(Dg_pole)
    print("Its O(omega^4) series is")
    sp.pprint(sp.expand(Dg_pole_series))
    print("G_0 =")
    sp.pprint(G0_pole)
    print("G_2 =")
    sp.pprint(G2_pole)
    print("G_4 =")
    sp.pprint(G4_pole)

    expect_zero("G_0 - (G_c + G_p)", sp.simplify(G0_pole - (Gc + Gp)))
    expect_zero("G_2 - G_p/Omega_g^2", sp.simplify(G2_pole - Gp / Omegag**2))
    expect_zero("G_4 - G_p/Omega_g^4", sp.simplify(G4_pole - Gp / Omegag**4))

    subbanner("III. Exact weak-anisotropy coefficients a_2 and a_4 on that branch")
    mu_mix, c_g, r = sp.symbols("mu_mix c_g r", positive=True, real=True)
    Gtot = sp.symbols("G_tot", positive=True, real=True)

    raw_a2 = sp.simplify((eps2 / chi**2).subs({G0: G0_pole, G2: G2_pole}))
    raw_a4 = sp.simplify((eps4 / chi**2).subs({G0: G0_pole, G2: G2_pole, G4: G4_pole}))

    subs_branch = {
        Gc + Gp: Gtot,
        Gp: c_g * Gtot,
        Gc: (1 - c_g) * Gtot,
        M0**2: mu_mix * Kpole * Gtot,
        OmegaQ**2: r * Omegag**2,
    }

    a2 = sp.simplify(raw_a2.subs(subs_branch))
    a4 = sp.simplify(raw_a4.subs(subs_branch))
    drift = sp.simplify(a4 - 2 * a2)

    print("a_2 =")
    sp.pprint(a2)
    print("a_4 =")
    sp.pprint(a4)
    print("a_4 - 2 a_2 =")
    sp.pprint(sp.factor(drift))

    expect_zero("a_2 - mu_mix*r*c_g", sp.simplify(a2 - mu_mix * r * c_g))
    expect_zero(
        "a_4 - mu_mix*r^2*c_g*(1-c_g)",
        sp.simplify(a4 - mu_mix * r**2 * c_g * (1 - c_g)),
    )
    expect_zero(
        "a_4 - 2 a_2 - mu_mix*r*c_g*(r*(1-c_g)-2)",
        sp.simplify(drift - mu_mix * r * c_g * (r * (1 - c_g) - 2)),
    )

    subbanner("IV. Immediate exact implications")
    print("Contact fraction c_g := G_pole / G_0 lies in [0,1] on a positive contact-plus-pole scalar lane.")
    print("Therefore:")
    print("  - a_2 is always nonnegative on the positive one-pole scalar branch,")
    print("  - a_4 is also nonnegative and vanishes on the pure-pole limit c_g = 1,")
    print("  - the initial support-demand drift is controlled only by r(1-c_g)-2.")
    print()
    print("Special limits:")
    print("  * pure scalar pole   c_g = 1 -> a_4 = 0, so the first induced mixing only generates epsilon_2;")
    print("  * pure scalar contact c_g = 0 -> a_2 = a_4 = 0, so the l=0 lane does not contaminate the quadrupole carrier;")
    print("  * mixed contact+pole -> both coefficients are nonzero, with their relative size set by r = Omega_Q^2/Omega_g^2.")

    banner("STAGE 265 LEDGER")
    print("1. The exact Stage-78 induced geometry-mixing coefficients are")
    print("      a_2 = Omega_Q^2 M_0^2 G_2 / (K_pole G_0^2),")
    print("      a_4 = Omega_Q^4 M_0^2 (G_0 G_4 - G_2^2) / (K_pole G_0^3).")
    print("2. On the first concrete scalar contact+pole branch they collapse to")
    print("      a_2 = mu_mix r c_g,   a_4 = mu_mix r^2 c_g (1-c_g).")
    print("3. So the initial drift of the actual grouped-P2 support demand is controlled only by")
    print("      r(1-c_g) - 2,")
    print("   where r = Omega_Q^2/Omega_g^2 and c_g is the scalar-lane pole fraction.")
    print("4. The first induced l=0<->l=2 mechanism is therefore already much more rigid than a generic")
    print("   O(chi^2) statement: it has one amplitude, one pole-ratio, and one scalar-pole-fraction.")


if __name__ == "__main__":
    main()
