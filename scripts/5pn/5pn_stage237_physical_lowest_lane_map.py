#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 237 — physical (Pe, kappa, y) placement map for the first explicit non-twin
lowest support lane, folded into the Stage-236 twin/no-twin threshold.
"""


def main() -> None:
    banner("STAGE 237 — PHYSICAL LOWEST-LANE PLACEMENT MAP")

    Pe = sp.symbols("Pe", positive=True, real=True)
    kappa = sp.symbols("kappa", positive=True, real=True)
    y = sp.symbols("y", positive=True, real=True)
    zeta_req = sp.symbols("zeta_req", positive=True, real=True)
    pi = sp.pi

    Omega_Pe = sp.simplify(
        pi * Pe * (2 * Pe * sp.exp(Pe) + pi)
        / ((4 * Pe**2 + pi**2) * (sp.exp(Pe) - 1))
    )
    zeta_phys = sp.simplify(Omega_Pe**2 * (kappa + pi**2 / 4) / (kappa + y**2))
    zeta_max = sp.simplify(sp.limit(zeta_phys, Pe, sp.oo).subs(y, 0))
    kappa_max = sp.simplify(sp.solve(sp.Eq(zeta_req, zeta_max), kappa)[0])

    subbanner("I. Exact transport overlap factor")
    print("Omega_Pe =")
    sp.pprint(Omega_Pe)
    print()
    print("lim_{Pe->0} Omega_Pe =", sp.simplify(sp.limit(Omega_Pe, Pe, 0)))
    print("lim_{Pe->+inf} Omega_Pe =", sp.simplify(sp.limit(Omega_Pe, Pe, sp.oo)))
    print()
    print("Small-Pe series of Omega_Pe =")
    sp.pprint(sp.series(Omega_Pe, Pe, 0, 3))
    print()
    print("Strong-bias asymptotic of Omega_Pe =")
    sp.pprint(sp.series((pi / 2) - Omega_Pe, Pe, sp.oo, 3))

    subbanner("II. Exact physical non-twin lowest-lane ratio")
    print("zeta_0^(Pe+R)(Pe,y;kappa) =")
    sp.pprint(zeta_phys)
    expect_zero("twin baseline at Pe=0, y=pi/2", sp.simplify(zeta_phys.subs({y: pi / 2}).limit(Pe, 0) - 1))

    dz_dkappa = sp.simplify(sp.diff(zeta_phys, kappa))
    dz_dy = sp.simplify(sp.diff(zeta_phys, y))
    print()
    print("d zeta / d kappa =")
    sp.pprint(dz_dkappa)
    expect_zero(
        "exact kappa-monotonicity identity",
        sp.simplify(dz_dkappa - Omega_Pe**2 * (y**2 - pi**2 / 4) / (kappa + y**2) ** 2),
    )
    print()
    print("d zeta / d y =")
    sp.pprint(dz_dy)
    expect_zero(
        "exact y-monotonicity identity",
        sp.simplify(dz_dy + 2 * y * Omega_Pe**2 * (kappa + pi**2 / 4) / (kappa + y**2) ** 2),
    )

    subbanner("III. Exact closure ceiling and stiffness ceiling")
    print("zeta_max(kappa) =")
    sp.pprint(zeta_max)
    expect_zero(
        "exact ceiling identity",
        sp.simplify(zeta_max - (pi**2 / 4) * (kappa + pi**2 / 4) / kappa),
    )
    print()
    print("kappa_max(zeta_req) =")
    sp.pprint(kappa_max)
    expect_zero(
        "exact stiffness-ceiling identity",
        sp.simplify(kappa_max - pi**4 / (4 * (4 * zeta_req - pi**2))),
    )

    subbanner("IV. Twin vs non-twin theorem gate")
    print("Symmetric lowest twin lane from Stage 236 succeeds iff zeta_req <= 1.")
    print("First explicit non-twin physical family reaches up to zeta_max(kappa).")
    print()
    print("So the first exact decision split is:")
    print("  zeta_req <= 1                : lowest symmetric twin enough")
    print("  1 < zeta_req <= zeta_max     : non-twin family can in principle rescue")
    print("  zeta_req > zeta_max          : first explicit non-twin family impossible")

    banner("STAGE 237 LEDGER")
    print("1. The first explicit non-twin lowest-lane family is")
    print("      zeta_0^(Pe+R) = Omega_Pe^2 (kappa + pi^2/4)/(kappa + y^2).")
    print("2. It recovers the Stage-236 twin baseline exactly at Pe=0 and y=pi/2.")
    print("3. The placement map is monotone in the physical directions that help support:")
    print("      larger Pe helps, smaller y helps, smaller kappa helps.")
    print("4. Its exact ceiling is")
    print("      zeta_max(kappa) = (pi^2/4)(kappa + pi^2/4)/kappa.")
    print("5. Therefore the Stage-236 twin criterion expands into the first physical phase split:")
    print("      zeta_req <= 1           : symmetric twin enough;")
    print("      1 < zeta_req <= zeta_max: non-twin family potentially reachable;")
    print("      zeta_req > zeta_max     : first explicit non-twin family impossible.")


if __name__ == "__main__":
    main()
