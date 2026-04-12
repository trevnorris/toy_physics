#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp
import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 254 — exact admissible one-pole conservative window on the refreshed Family-1 branch.

This script keeps the exact Lambda_EM-refreshed Family-1 ceiling and derives the
exact map from a generic isotropic one-pole conservative carrier

    Y_Q^cons(omega) = c_geom + c_pole / (1 - omega^2/Omega_Q^2),

with c_geom + c_pole = 1, into the demand variables

    rho_alpha,
    zeta_req,
    Pi_tr / C_mix.

It then converts the exact Family-1 support ceiling into an exact admissible
interval for c_pole.
"""


def omega_pe(Pe: mp.mpf) -> mp.mpf:
    if abs(Pe) < mp.mpf("1e-30"):
        return mp.mpf("1")
    return mp.pi * Pe * (2 * Pe * mp.e**Pe + mp.pi) / (
        (4 * Pe**2 + mp.pi**2) * (mp.e**Pe - 1)
    )


def family1_refreshed_numbers() -> dict[str, mp.mpf]:
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
    zeta_suff_chi = A_F1 * omega_pe(Pe_suff_chi) ** 2

    return {
        "Lambda_EM": Lambda_EM,
        "zeta_max": zeta_max,
        "zeta_suff_chi": zeta_suff_chi,
    }


def main() -> None:
    banner("STAGE 254 — EXACT FAMILY-1 ADMISSIBLE ONE-POLE WINDOW")

    c_pole, eps_blk, zeta = sp.symbols("c_pole eps_blk zeta", real=True)
    rho_alpha = sp.simplify(1 / (1 - c_pole))
    zeta_req = sp.simplify((rho_alpha - 1) / (1 - eps_blk * (2 - rho_alpha)))
    Pi_over_Cmix = sp.simplify((1 + (1 - 2 * eps_blk) * zeta_req) / (1 - eps_blk * zeta_req))

    subbanner("I. Exact one-pole map into rho_alpha, zeta_req, and Pi_tr/C_mix")
    print("rho_alpha(c_pole) =")
    sp.pprint(rho_alpha)
    print("zeta_req(c_pole; eps_blk) =")
    sp.pprint(sp.factor(zeta_req))
    print("Pi_tr / C_mix =")
    sp.pprint(Pi_over_Cmix)
    expect_zero("Pi/C_mix - rho_alpha", sp.simplify(Pi_over_Cmix - rho_alpha))

    c_pole_star = sp.symbols("c_pole_star", real=True)
    c_pole_max = sp.simplify(sp.solve(sp.Eq(zeta_req, zeta), c_pole)[0])
    print("c_pole(zeta_req = zeta) =")
    sp.pprint(sp.factor(c_pole_max))

    c_geom_min = sp.simplify(1 - c_pole_max)
    print("c_geom_min(zeta ceiling) =")
    sp.pprint(sp.factor(c_geom_min))
    expect_zero(
        "c_geom_min - 1/Q(zeta;eps_blk)",
        sp.simplify(c_geom_min - (1 - eps_blk * zeta) / (1 + (1 - 2 * eps_blk) * zeta)),
    )

    nums = family1_refreshed_numbers()
    zeta_max_num = nums["zeta_max"]
    zeta_suff_num = nums["zeta_suff_chi"]

    rho_max_num = (1 + zeta_max_num)  # unblocked branch
    rho_suff_num = (1 + zeta_suff_num)
    c_pole_max_num = mp.mpf("1") - 1 / rho_max_num
    c_pole_suff_num = mp.mpf("1") - 1 / rho_suff_num
    c_geom_min_num = 1 / rho_max_num

    subbanner("II. Exact Lambda_EM-refreshed Family-1 ceiling in one-pole language")
    print(f"Lambda_EM ≈ {nums['Lambda_EM']}")
    print(f"zeta_max^(F1) ≈ {zeta_max_num}")
    print(f"rho_max^(F1) = 1 + zeta_max^(F1) ≈ {rho_max_num}")
    print(f"c_pole,max^(F1) ≈ {c_pole_max_num}")
    print(f"c_geom,min^(F1) ≈ {c_geom_min_num}")
    print()
    print(f"Guaranteed-success shell-weighted threshold zeta_suff^(chi) ≈ {zeta_suff_num}")
    print(f"rho_suff^(chi) ≈ {rho_suff_num}")
    print(f"c_pole,suff^(chi) ≈ {c_pole_suff_num}")

    subbanner("III. Canonical minimal isotropic point inside the exact window")
    c_pole_can = mp.mpf("1") / 4
    rho_can = 1 / (1 - c_pole_can)
    zeta_can = c_pole_can / (1 - c_pole_can)
    print(f"c_pole^(can) = {c_pole_can}")
    print(f"rho_alpha^(can) = {rho_can}")
    print(f"zeta_req^(can, eps_blk=0) = {zeta_can}")
    print(f"Hard c_pole margin ≈ {c_pole_max_num - c_pole_can}")
    print(f"Hard rho_alpha margin ≈ {rho_max_num - rho_can}")

    banner("STAGE 254 LEDGER")
    print("1. A generic isotropic one-pole conservative carrier is parameterized exactly by")
    print("      rho_alpha = 1/(1-c_pole),")
    print("      zeta_req  = c_pole / [1-eps_blk-(1-2 eps_blk)c_pole],")
    print("      Pi_tr/C_mix = rho_alpha.")
    print(f"2. On the exact Lambda_EM-refreshed Family-1 branch, the hard unblocked ceiling is")
    print(f"      c_pole < {c_pole_max_num}")
    print("   with guaranteed shell-weighted success already for")
    print(f"      c_pole <= {c_pole_suff_num}.")
    print("3. So the explicit support/source side remains far from the canonical")
    print("   c_pole = 1/4 branch and still leaves a large admissible deformation window.")


if __name__ == "__main__":
    main()
