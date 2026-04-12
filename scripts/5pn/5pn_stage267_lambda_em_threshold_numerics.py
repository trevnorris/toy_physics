#!/usr/bin/env python3
from __future__ import annotations

import mpmath as mp

from fivepn_stage265_267_common import banner, subbanner, family1_refreshed_numbers

"""
Stage 267 — exact numerical thresholds for the first actual l=0 <-> l=2
geometry-mixing mechanism on the Lambda_EM-refreshed Family-1 branch.

This script evaluates the exact Stage-266 thresholds on the physically relevant
Family-1 ceiling imported from the exact Lambda_EM branch.  It makes the first
actual quantitative statement about how large the scalar/geometry pole ratio
must be before the induced mixing can threaten the twin-safe or Family-1 strips.
"""


def main() -> None:
    banner("STAGE 267 — LAMBDA_EM-REFRESHED THRESHOLD NUMERICS")

    nums = family1_refreshed_numbers()
    cstar_max = mp.mpf(nums["c_pole_max_0"])
    cstar_suff = mp.mpf(nums["c_pole_suff_0"])

    u_twin_demand = mp.mpf("2")
    u_twin_shrink = mp.mpf("4")
    u_twin_fail = mp.mpf("4") + 2 * mp.sqrt(2)

    def u_F1_shrink(cstar: mp.mpf) -> mp.mpf:
        return 8 * cstar

    def u_F1_fail(cstar: mp.mpf) -> mp.mpf:
        return 8 * cstar + 4 * mp.sqrt(cstar * (4 * cstar - 1))

    u_F1_shrink_max = u_F1_shrink(cstar_max)
    u_F1_fail_max = u_F1_fail(cstar_max)
    u_F1_shrink_suff = u_F1_shrink(cstar_suff)
    u_F1_fail_suff = u_F1_fail(cstar_suff)

    subbanner("I. Exact threshold values on the refreshed Family-1 branch")
    print(f"c_*^(max)   = c_pole,max^(F1)(0)  ≈ {cstar_max}")
    print(f"c_*^(suff)  = c_pole,suff^(F1)(0) ≈ {cstar_suff}")
    print()
    print(f"u_twin^(demand) = {u_twin_demand}")
    print(f"u_twin^(shrink) = {u_twin_shrink}")
    print(f"u_twin^(fail)   = 4 + 2 sqrt(2) ≈ {u_twin_fail}")
    print()
    print(f"u_F1^(shrink,max) = 8 c_*^(max) ≈ {u_F1_shrink_max}")
    print(f"u_F1^(fail,max)   = 8 c_*^(max) + 4 sqrt[c_*(4 c_* - 1)] ≈ {u_F1_fail_max}")
    print()
    print(f"u_F1^(shrink,suff) = 8 c_*^(suff) ≈ {u_F1_shrink_suff}")
    print(f"u_F1^(fail,suff)   = 8 c_*^(suff) + 4 sqrt[c_*(4 c_* - 1)] ≈ {u_F1_fail_suff}")

    subbanner("II. Required frequency-ratio thresholds for representative scalar-pole fractions")
    cg_samples = [mp.mpf("0.00"), mp.mpf("0.25"), mp.mpf("0.50"), mp.mpf("0.75"), mp.mpf("0.90")]

    print("For u = r (1-c_g), the required r-threshold is r_crit = u_crit / (1-c_g).")
    print()
    print("c_g      r_demand     r_twin_shrink     r_twin_fail      r_F1_shrink(max)   r_F1_fail(max)")
    for cg in cg_samples:
        denom = 1 - cg
        r_demand = u_twin_demand / denom
        r_twin_s = u_twin_shrink / denom
        r_twin_f = u_twin_fail / denom
        r_F1_s = u_F1_shrink_max / denom
        r_F1_f = u_F1_fail_max / denom
        print(
            f"{float(cg):0.2f}   {float(r_demand):12.6f}   {float(r_twin_s):15.6f}   {float(r_twin_f):15.6f}   {float(r_F1_s):18.6f}   {float(r_F1_f):15.6f}"
        )

    subbanner("III. Immediate reading")
    print("The first actual l=0<->l=2 geometry-mixing mechanism becomes dangerous only at very large")
    print("quadrupole/scalar pole ratios once the scalar lane carries any sizable contact fraction.")
    print()
    print("Examples on the exact Lambda_EM-refreshed hard Family-1 ceiling:")
    print(f"  c_g = 1/2  ->  r_twin,fail ≈ {u_twin_fail / mp.mpf('0.5')},   r_F1,fail ≈ {u_F1_fail_max / mp.mpf('0.5')}")
    print(f"  c_g = 3/4  ->  r_twin,fail ≈ {u_twin_fail / mp.mpf('0.25')},  r_F1,fail ≈ {u_F1_fail_max / mp.mpf('0.25')}")
    print()
    print("So unless the scalar geometry pole is *much* faster than the grouped quadrupole pole, the")
    print("first induced mixing mechanism will not by itself drive the actual isotropic branch out of")
    print("the universal twin-safe or exact Family-1 corridors.")

    banner("STAGE 267 LEDGER")
    print("1. On the exact Lambda_EM-refreshed Family-1 branch the hard thresholds are approximately")
    print(f"      u_twin^(fail) ≈ {u_twin_fail},   u_F1^(fail,max) ≈ {u_F1_fail_max}.")
    print("2. Because u = r(1-c_g), the required pole-ratio thresholds scale like 1/(1-c_g).")
    print("3. Therefore any scalar lane with a substantial contact fraction c_g is strongly protected:")
    print("   it would need a very large scalar/quadrupole pole separation before the first induced")
    print("   l=0<->l=2 mixing mechanism could threaten the exact twin-safe or Family-1 strips.")
    print("4. This makes the simplest one-pole scalar-lane mixing branch look much more like a small")
    print("   quantitative correction than like the actual source of 5PN failure.")


if __name__ == "__main__":
    main()
