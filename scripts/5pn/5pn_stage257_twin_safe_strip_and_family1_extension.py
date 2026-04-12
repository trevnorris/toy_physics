#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp
import sympy as sp

from fivepn_stage199_201_common import banner, subbanner, expect_zero

"""
Stage 257 — universal twin-safe strip and the exact extension supplied by the
Lambda_EM-refreshed Family-1 ceiling.
"""


def refreshed_family1_ceiling() -> tuple[mp.mpf, mp.mpf]:
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
    rho_max = 1 + zeta_max
    c_pole_max = 1 - 1 / rho_max
    return rho_max, c_pole_max


def main() -> None:
    banner("STAGE 257 — UNIVERSAL TWIN-SAFE STRIP AND FAMILY-1 EXTENSION")

    c_pole, eps_blk = sp.symbols("c_pole eps_blk", real=True)
    rho_alpha = sp.simplify(1 / (1 - c_pole))
    zeta_req = sp.simplify((rho_alpha - 1) / (1 - eps_blk * (2 - rho_alpha)))
    c_half = sp.simplify(sp.solve(sp.Eq(zeta_req, 1), c_pole)[0])

    subbanner("I. Exact universal twin-safe strip")
    print("zeta_req(c_pole; eps_blk) =")
    sp.pprint(sp.factor(zeta_req))
    print("c_pole at zeta_req = 1 =")
    sp.pprint(c_half)
    expect_zero("c_half - 1/2", sp.simplify(c_half - sp.Rational(1, 2)))
    print("So every isotropic one-pole branch with c_pole <= 1/2 stays in the")
    print("symmetric-lowest-twin-safe regime zeta_req <= 1, independently of blocking.")

    rho_max_num, c_pole_max_num = refreshed_family1_ceiling()
    extension_num = c_pole_max_num - mp.mpf("1") / 2
    rho_extension_num = rho_max_num - mp.mpf("2")

    subbanner("II. Exact Lambda_EM Family-1 extension beyond the universal strip")
    print(f"Universal twin-safe endpoint: c_pole = 1/2, rho_alpha = 2.")
    print(f"Refreshed Family-1 hard ceiling: c_pole < {c_pole_max_num}, rho_alpha < {rho_max_num}.")
    print(f"Additional c_pole headroom beyond the universal twin-safe strip ≈ {extension_num}")
    print(f"Additional rho_alpha headroom beyond the universal twin-safe strip ≈ {rho_extension_num}")

    subbanner("III. Canonical actual isotropic branch placement")
    c_can = mp.mpf("1") / 4
    print(f"Canonical branch: c_pole = {c_can}, rho_alpha = {mp.mpf('4')/3}.")
    print(f"Margin to universal twin-safe endpoint c_pole = 1/2: {mp.mpf('1')/2 - c_can}")
    print(f"Margin to exact Family-1 hard ceiling: {c_pole_max_num - c_can}")

    banner("STAGE 257 LEDGER")
    print("1. The universal support theorem already gives an exact blocking-independent safe strip")
    print("      c_pole <= 1/2   <=>   zeta_req <= 1.")
    print("2. The refreshed Lambda_EM Family-1 branch extends that strip to the much larger exact window")
    print(f"      c_pole < {c_pole_max_num},   equivalently   rho_alpha < {rho_max_num}.")
    print("3. So the actual isotropic c_pole = 1/4 branch is not merely safe; it sits deeply inside")
    print("   both the universal twin-safe strip and the exact Family-1 extension beyond it.")


if __name__ == "__main__":
    main()
