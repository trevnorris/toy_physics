#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp

from fivepn_stage199_201_common import banner, subbanner

"""
Stage 252 — direct wall-depth verdict after replacing the shorthand reference
geometry by the exact carried EM-branch ratio.
"""


def main() -> None:
    banner("STAGE 252 — EXACT LAMBDA_EM WALL-DEPTH VERDICT")

    mp.mp.dps = 80

    # Exact refreshed threshold window from Stage 251
    Theta_fail = mp.mpf('0.00036309892670268561443639899600628005929446369050336')
    Theta_suff = mp.mpf('0.042149534156997728721243988650267664871034170059267')

    # Carried explicit wall-depth extraction from the concrete Family-1 wall profile
    Theta_chi = mp.mpf('4.06863235008162')
    Theta_J = mp.mpf('0.927552032539308')

    Pe_suff_chi = Theta_chi / Theta_suff
    Pe_fail_chi = Theta_chi / Theta_fail
    Pe_suff_J = Theta_J / Theta_suff
    Pe_fail_J = Theta_J / Theta_fail

    subbanner("I. Explicit wall-depth data (unchanged by the geometry refresh)")
    print(f"Theta_w^(chi) ≈ {Theta_chi} * lambda_mu^2")
    print(f"Theta_w^(J)   ≈ {Theta_J} * lambda_mu^2")

    subbanner("II. Updated Pe_req comparison on the exact Lambda_EM branch")
    print(f"Pe_suff^(chi) ≈ {Pe_suff_chi} * lambda_mu^2")
    print(f"Pe_fail^(chi) ≈ {Pe_fail_chi} * lambda_mu^2")
    print(f"Pe_suff^(J)   ≈ {Pe_suff_J} * lambda_mu^2")
    print(f"Pe_fail^(J)   ≈ {Pe_fail_J} * lambda_mu^2")

    banner("STAGE 252 LEDGER")
    print("1. The explicit shell-weighted wall-depth extraction is unchanged by the exact")
    print("   Lambda_EM refresh; only the threshold window moves slightly.")
    print("2. On the natural wall datum the refreshed branch still succeeds automatically")
    print("   for modest Pe_req and fails only for anomalously large Pe_req.")
    print("3. So the qualitative verdict survives intact:")
    print("      wall-depth is still not the natural bottleneck of the explicit Family-1 branch.")
    print("4. The remaining serious open quantity is still the quadrupole-side demand Pe_req.")


if __name__ == "__main__":
    main()
