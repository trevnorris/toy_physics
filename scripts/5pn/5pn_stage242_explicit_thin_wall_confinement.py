#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 242 — explicit thin-wall confinement branch on the parent equilibrium-aligned support/source
sector.

This script:
1. turns the abstract parent loading amplitude g_phi into the wall family parameter V0/ell,
2. records the exact shell integral I1 for the first thin-wall wall family,
3. derives the exact equilibrium gain G_eq,
4. isolates the centered-wall and thin-wall limits.
"""


def main() -> None:
    banner("STAGE 242 — EXPLICIT THIN-WALL CONFINEMENT BRANCH")

    V0, a, ell, KX = sp.symbols("V0 a ell K_X", positive=True, real=True)
    J1, J2, J3 = sp.symbols("J_1 J_2 J_3", real=True)

    subbanner("I. Parent loading amplitude on the wall family")
    gphi = sp.simplify(V0 / ell)
    print("For  V_conf(r;a) = V0 f((r-a)/ell),")
    print("the support-loading amplitude is")
    print("g_phi =")
    sp.pprint(gphi)

    subbanner("II. Exact shell integral I1")
    I1 = sp.simplify(4 * sp.pi * ell * (a**2 * J1 + 2 * a * ell * J2 + ell**2 * J3))
    print("I1 =")
    sp.pprint(I1)

    I1_centered = sp.simplify(I1.subs({J2: 0}))
    print()
    print("For a centered symmetric wall layer (J2 = 0):")
    print("I1_centered =")
    sp.pprint(I1_centered)
    expect_zero(
        "centered-shell identity",
        sp.simplify(I1_centered - 4 * sp.pi * ell * (a**2 * J1 + ell**2 * J3)),
    )

    subbanner("III. Exact equilibrium gain")
    G_eq = sp.simplify(gphi**2 * I1 / KX)
    print("G_eq =")
    sp.pprint(G_eq)

    G_eq_expanded = sp.simplify(4 * sp.pi * V0**2 * (a**2 * J1 / ell + 2 * a * J2 + ell * J3) / KX)
    expect_zero("expanded G_eq identity", sp.simplify(G_eq - G_eq_expanded))

    G_eq_centered = sp.simplify(G_eq.subs({J2: 0}))
    print()
    print("Centered-wall branch:")
    print("G_eq_centered =")
    sp.pprint(G_eq_centered)

    G_eq_tw = sp.simplify(4 * sp.pi * a**2 * V0**2 * J1 / (KX * ell))
    remainder = sp.simplify(G_eq_centered - (G_eq_tw + 4 * sp.pi * V0**2 * ell * J3 / KX))
    expect_zero("thin-wall decomposition", remainder)

    subbanner("IV. Thin-wall leading term")
    print("Leading thin-wall gain G_eq^(tw) =")
    sp.pprint(G_eq_tw)

    banner("STAGE 242 LEDGER")
    print("1. The explicit confinement branch fixes the parent loading amplitude to g_phi = V0/ell.")
    print("2. The exact shell integral is")
    print("      I1 = 4 pi ell [ a^2 J1 + 2 a ell J2 + ell^2 J3 ].")
    print("3. The exact equilibrium gain is")
    print("      G_eq = 4 pi V0^2 [ a^2 J1/ell + 2 a J2 + ell J3 ] / K_X.")
    print("4. On the centered branch J2 = 0, the thin-wall leading term is")
    print("      G_eq^(tw) = 4 pi a^2 V0^2 J1 / (K_X ell).")
    print("5. So the parent support/source theorem is now a direct wall-amplitude problem on the first explicit wall family.")


if __name__ == "__main__":
    main()
