#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp
import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 255 — exact Family-1 ceiling translated into the geometry-lane contamination variables.

Uses the stage-75 pole-fraction map

    c_pole = (1 + eps4) / (4 (1 + eps2)^2)

and the exact Lambda_EM-refreshed Family-1 admissible c_pole window to derive the
precise safe region in the (eps2, eps4) plane.
"""


def family1_refreshed_numbers() -> tuple[mp.mpf, mp.mpf]:
    # Hard and guaranteed-success pole-fraction ceilings from Stage 254 numbers.
    mp.mp.dps = 80
    x01 = mp.besseljzero(0, 1)
    Lambda_EM = mp.sqrt(2) * mp.pi / x01
    Lambda_ell = 20 * Lambda_EM
    eta = Lambda_ell
    kappa = mp.mpf("9") * Lambda_ell**2 / 5
    f = lambda y: y * mp.tan(y) - eta
    lo = mp.mpf("1.4")
    hi = mp.mpf("1.56")
    for _ in range(400):
        mid = (lo + hi) / 2
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    y = (lo + hi) / 2
    A_F1 = (kappa + mp.pi**2 / 4) / (kappa + y**2)
    zeta_max = A_F1 * mp.pi**2 / 4

    Theta_chi = mp.mpf("4.06863235008162")
    Theta_suff = mp.mpf("0.042149534156997728721243988650267664871034170059267")
    Pe_suff_chi = Theta_chi / Theta_suff

    def omega_pe(Pe: mp.mpf) -> mp.mpf:
        if abs(Pe) < mp.mpf("1e-30"):
            return mp.mpf("1")
        return mp.pi * Pe * (2 * Pe * mp.e**Pe + mp.pi) / (
            (4 * Pe**2 + mp.pi**2) * (mp.e**Pe - 1)
        )

    zeta_suff = A_F1 * omega_pe(Pe_suff_chi) ** 2
    c_pole_max = mp.mpf("1") - 1 / (1 + zeta_max)
    c_pole_suff = mp.mpf("1") - 1 / (1 + zeta_suff)
    return c_pole_max, c_pole_suff


def main() -> None:
    banner("STAGE 255 — GEOMETRY-LANE CONTAMINATION CEILING")

    eps2, eps4, cstar = sp.symbols("eps2 eps4 c_star", real=True)
    c_pole = sp.simplify((1 + eps4) / (4 * (1 + eps2) ** 2))
    eps4_max = sp.simplify(sp.solve(sp.Eq(c_pole, cstar), eps4)[0])
    eps2_floor = sp.simplify(sp.solve(sp.Eq(c_pole, cstar), eps2)[1])  # principal + branch

    subbanner("I. Exact map from geometry contamination to pole fraction")
    print("c_pole(eps2, eps4) =")
    sp.pprint(c_pole)
    print("eps4 ceiling at fixed eps2 and c_star =")
    sp.pprint(eps4_max)
    print("principal eps2 floor at fixed eps4 and c_star =")
    sp.pprint(eps2_floor)
    expect_zero(
        "roundtrip c_pole(eps4_max)-c_star",
        sp.simplify(c_pole.subs(eps4, eps4_max) - cstar),
    )

    c_pole_max_num, c_pole_suff_num = family1_refreshed_numbers()

    hard_eps4_at_0 = 4 * c_pole_max_num - 1
    suff_eps4_at_0 = 4 * c_pole_suff_num - 1
    hard_eps2_at_0 = mp.sqrt(1 / (4 * c_pole_max_num)) - 1
    suff_eps2_at_0 = mp.sqrt(1 / (4 * c_pole_suff_num)) - 1

    subbanner("II. Exact Lambda_EM-refreshed safe strip in (eps2, eps4)")
    print(f"Hard pole-fraction ceiling c_pole,max ≈ {c_pole_max_num}")
    print(f"Guaranteed-success ceiling    c_pole,suff ≈ {c_pole_suff_num}")
    print()
    print("At eps2 = 0:")
    print(f"  eps4 < {hard_eps4_at_0}   (hard ceiling)")
    print(f"  eps4 < {suff_eps4_at_0}   (guaranteed success)")
    print()
    print("At eps4 = 0 on the physical principal branch 1+eps2 > 0:")
    print(f"  eps2 > {hard_eps2_at_0}   (hard ceiling)")
    print(f"  eps2 > {suff_eps2_at_0}   (guaranteed success)")

    subbanner("III. Canonical point and margin")
    c_can = mp.mpf("1") / 4
    eps4_can_margin = hard_eps4_at_0
    eps2_can_margin = -hard_eps2_at_0
    print(f"Canonical isotropic point: c_pole = 1/4, eps2 = 0, eps4 = 0.")
    print(f"Headroom in eps4 at eps2=0: {eps4_can_margin}")
    print(f"Headroom toward negative eps2 at eps4=0: {eps2_can_margin}")

    banner("STAGE 255 LEDGER")
    print("1. The exact Family-1 support/source ceiling translates to the geometry-lane region")
    print("      1 + eps4 < 4 c_star (1 + eps2)^2")
    print("   with c_star taken as either the hard ceiling or the guaranteed-success ceiling.")
    print("2. On the Lambda_EM-refreshed branch, the hard unblocked ceiling still allows")
    print(f"      eps4 < {hard_eps4_at_0}   at eps2 = 0,")
    print(f"      eps2 > {hard_eps2_at_0}   at eps4 = 0.")
    print("3. So the explicit support/source side remains tolerant to order-one geometry-lane")
    print("   contamination around the canonical c_pole = 1/4 branch.")


if __name__ == "__main__":
    main()
