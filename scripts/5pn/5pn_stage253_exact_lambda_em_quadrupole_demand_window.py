#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp

from fivepn_stage199_201_common import banner, subbanner

"""
Stage 253 — refreshed Family-1 quadrupole-demand window in zeta_req and Pi_tr
after replacing the reference freeze 37/20 by the exact carried EM-branch ratio.
"""


def omega_pe(Pe: mp.mpf) -> mp.mpf:
    if abs(Pe) < mp.mpf('1e-30'):
        return mp.mpf('1')
    return mp.pi * Pe * (2 * Pe * mp.e**Pe + mp.pi) / (
        (4 * Pe**2 + mp.pi**2) * (mp.e**Pe - 1)
    )


def Q(zeta: mp.mpf, eps_blk: mp.mpf) -> mp.mpf:
    return (1 + (1 - 2 * eps_blk) * zeta) / (1 - eps_blk * zeta)


def main() -> None:
    banner("STAGE 253 — EXACT LAMBDA_EM QUADRUPOLE-DEMAND WINDOW")

    mp.mp.dps = 80
    x01 = mp.besseljzero(0, 1)
    Lambda_EM = mp.sqrt(2) * mp.pi / x01
    Lambda_ell = 20 * Lambda_EM
    eta = Lambda_ell
    kappa = mp.mpf('9') * Lambda_ell**2 / 5

    # Solve y tan y = eta on the principal positive branch
    f = lambda y: y * mp.tan(y) - eta
    lo = mp.mpf('1.4')
    hi = mp.mpf('1.56')
    for _ in range(300):
        mid = (lo + hi) / 2
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    y_F1 = (lo + hi) / 2

    A_F1 = (kappa + mp.pi**2 / 4) / (kappa + y_F1**2)
    zeta_max = A_F1 * mp.pi**2 / 4

    Theta_chi = mp.mpf('4.06863235008162')
    Theta_J = mp.mpf('0.927552032539308')
    Theta_fail = mp.mpf('0.00036309892670268561443639899600628005929446369050336')
    Theta_suff = mp.mpf('0.042149534156997728721243988650267664871034170059267')

    Pe_suff_chi = Theta_chi / Theta_suff
    Pe_fail_chi = Theta_chi / Theta_fail
    Pe_suff_J = Theta_J / Theta_suff
    Pe_fail_J = Theta_J / Theta_fail

    zeta_suff_chi = A_F1 * omega_pe(Pe_suff_chi) ** 2
    zeta_fail_chi = A_F1 * omega_pe(Pe_fail_chi) ** 2
    zeta_suff_J = A_F1 * omega_pe(Pe_suff_J) ** 2
    zeta_fail_J = A_F1 * omega_pe(Pe_fail_J) ** 2

    eps_blk = mp.mpf('0')
    Pi_suff_chi = Q(zeta_suff_chi, eps_blk)
    Pi_fail_chi = Q(zeta_fail_chi, eps_blk)
    Pi_suff_J = Q(zeta_suff_J, eps_blk)
    Pi_fail_J = Q(zeta_fail_J, eps_blk)
    Pi_max = Q(zeta_max, eps_blk)
    eps_crit = 1 / zeta_max

    subbanner("I. Refreshed Family-1 support-ratio map")
    print(f"Lambda_EM ≈ {Lambda_EM}")
    print(f"y_F1 ≈ {y_F1}")
    print(f"A_F1 ≈ {A_F1}")
    print(f"zeta_max^(F1) ≈ {zeta_max}")

    subbanner("II. Updated zeta windows at lambda_mu = 1")
    print(f"zeta_suff^(chi)(1) ≈ {zeta_suff_chi}")
    print(f"zeta_fail^(chi)(1) ≈ {zeta_fail_chi}")
    print(f"zeta_suff^(J)(1)   ≈ {zeta_suff_J}")
    print(f"zeta_fail^(J)(1)   ≈ {zeta_fail_J}")

    subbanner("III. Updated Pi_tr window at eps_blk = 0")
    print(f"Pi_suff^(chi) / C_mix ≈ {Pi_suff_chi}")
    print(f"Pi_fail^(chi) / C_mix ≈ {Pi_fail_chi}")
    print(f"Pi_suff^(J)   / C_mix ≈ {Pi_suff_J}")
    print(f"Pi_fail^(J)   / C_mix ≈ {Pi_fail_J}")
    print(f"Pi_max^(F1)   / C_mix ≈ {Pi_max}")
    print(f"eps_blk critical ceiling ≈ {eps_crit}")

    banner("STAGE 253 LEDGER")
    print("1. The exact Lambda_EM refresh leaves the explicit quadrupole-demand picture")
    print("   qualitatively unchanged.")
    print("2. The natural shell-weighted branch still pushes the guaranteed-success threshold")
    print("   extremely close to the hard Family-1 ceiling in zeta_req and Pi_tr language.")
    print("3. So the explicit Family-1 phase remains bottlenecked by the selected quadrupole")
    print("   demand side, not by wall-depth supply.")


if __name__ == "__main__":
    main()
