#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp

from fivepn_stage199_201_common import banner, subbanner

"""
Stage 251 — explicit Family-1 threshold window using the exact carried geometry
ratio Lambda_EM = sqrt(2) pi / x_01 instead of the shorthand 37/20.
"""


def Delta0_fn(kappa: mp.mpf, eta: mp.mpf) -> mp.mpf:
    alpha = mp.sqrt(kappa)
    return eta * (mp.cosh(alpha) - 1) / (
        kappa * (eta * mp.cosh(alpha) + alpha * mp.sinh(alpha))
    )


def DeltaInf_fn(kappa: mp.mpf, eta: mp.mpf) -> mp.mpf:
    alpha = mp.sqrt(kappa)
    return (eta * mp.sinh(alpha) + alpha * (mp.cosh(alpha) - 1)) / (
        alpha * (eta * mp.cosh(alpha) + alpha * mp.sinh(alpha))
    )


def main() -> None:
    banner("STAGE 251 — EXACT LAMBDA_EM THRESHOLD WINDOW")

    mp.mp.dps = 80
    x01 = mp.besseljzero(0, 1)

    Lambda_EM = mp.sqrt(2) * mp.pi / x01
    Lambda_ell = 20 * Lambda_EM
    eta = Lambda_ell
    kappa = mp.mpf('9') * Lambda_ell**2 / 5

    Delta0 = Delta0_fn(kappa, eta)
    DeltaInf = DeltaInf_fn(kappa, eta)

    Upsilon_fail = 1 / (Lambda_ell**2 * DeltaInf)
    Upsilon_suff = 1 / (Lambda_ell**2 * Delta0)
    Theta_fail = Upsilon_fail / 100
    Theta_suff = Upsilon_suff / 100

    old_Upsilon_fail = mp.mpf('0.0362605617972939')
    old_Upsilon_suff = mp.mpf('4.21495341569977')

    subbanner("I. Exact refreshed branch parameters")
    print(f"Lambda_EM ≈ {Lambda_EM}")
    print(f"Lambda_ell = eta ≈ {Lambda_ell}")
    print(f"kappa ≈ {kappa}")

    subbanner("II. Exact threshold window on the refreshed branch")
    print(f"Delta_0 ≈ {Delta0}")
    print(f"Delta_inf ≈ {DeltaInf}")
    print()
    print(f"Upsilon_fail ≈ {Upsilon_fail} * Pe_req")
    print(f"Upsilon_suff ≈ {Upsilon_suff} * Pe_req")
    print()
    print(f"Theta_fail ≈ {Theta_fail} * Pe_req")
    print(f"Theta_suff ≈ {Theta_suff} * Pe_req")

    subbanner("III. Shift relative to the 37/20 reference freeze")
    print(f"relative shift in Upsilon_fail ≈ {(Upsilon_fail - old_Upsilon_fail) / old_Upsilon_fail}")
    print(f"relative shift in Upsilon_suff ≈ {(Upsilon_suff - old_Upsilon_suff) / old_Upsilon_suff}")

    banner("STAGE 251 LEDGER")
    print("1. Replacing 37/20 by the exact EM-branch ratio changes the explicit Family-1")
    print("   fail threshold slightly but leaves the success threshold essentially unchanged.")
    print("2. The refreshed threshold window is")
    print(f"      Theta_fail ≈ {Theta_fail} Pe_req,")
    print(f"      Theta_suff ≈ {Theta_suff} Pe_req.")
    print("3. So the explicit-branch bottleneck analysis can now be rerun without the")
    print("   reference-rational approximation.")


if __name__ == "__main__":
    main()
