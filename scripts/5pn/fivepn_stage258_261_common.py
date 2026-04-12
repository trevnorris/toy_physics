from __future__ import annotations

import mpmath as mp

from fivepn_stage199_201_common import banner, subbanner, expect_zero


def omega_pe(Pe: mp.mpf) -> mp.mpf:
    if abs(Pe) < mp.mpf("1e-30"):
        return mp.mpf("1")
    return mp.pi * Pe * (2 * Pe * mp.e**Pe + mp.pi) / (
        (4 * Pe**2 + mp.pi**2) * (mp.e**Pe - 1)
    )


def family1_refreshed_numbers() -> dict[str, mp.mpf]:
    """
    Recompute the exact Lambda_EM-refreshed Family-1 support ceiling from the
    same finite-throat branch used in Stage 254, so later stages do not need to
    hard-code any of the numerics.
    """
    mp.mp.dps = 100
    x01 = mp.besseljzero(0, 1)
    Lambda_EM = mp.sqrt(2) * mp.pi / x01
    Lambda_ell = 20 * Lambda_EM
    eta = Lambda_ell
    kappa = mp.mpf("9") * Lambda_ell**2 / 5

    # Principal positive Robin root y*tan(y)=eta.
    f = lambda y: y * mp.tan(y) - eta
    lo = mp.mpf("1.4")
    hi = mp.mpf("1.56")
    for _ in range(500):
        mid = (lo + hi) / 2
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    y = (lo + hi) / 2

    A_F1 = (kappa + mp.pi**2 / 4) / (kappa + y**2)
    zeta_max = A_F1 * mp.pi**2 / 4
    c_pole_max_0 = zeta_max / (1 + zeta_max)
    eps_blk_crit = 1 / zeta_max

    Theta_chi = mp.mpf("4.06863235008162")
    Theta_suff = mp.mpf("0.042149534156997728721243988650267664871034170059267")
    Pe_suff_chi = Theta_chi / Theta_suff
    zeta_suff_chi = A_F1 * omega_pe(Pe_suff_chi) ** 2
    c_pole_suff_0 = zeta_suff_chi / (1 + zeta_suff_chi)

    return {
        "Lambda_EM": Lambda_EM,
        "zeta_max": zeta_max,
        "c_pole_max_0": c_pole_max_0,
        "c_pole_suff_0": c_pole_suff_0,
        "eps_blk_crit": eps_blk_crit,
    }
