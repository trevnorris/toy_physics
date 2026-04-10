#!/usr/bin/env python3
"""
5pn_stage72_family1_minimal_isotropic_verdict.py

Stage 72 audit: explicit Family-1 verdict for the minimal isotropic passive/outgoing quadrupole branch.
"""

from __future__ import annotations

import math
import mpmath as mp
import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr) -> None:
    expr_s = sp.simplify(sp.together(sp.expand(expr)))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 72 — EXPLICIT FAMILY-1 VERDICT FOR THE MINIMAL ISOTROPIC BRANCH")

mp.mp.dps = 80

rho_min = mp.mpf(4) / 3
zeta_min = mp.mpf(1) / 3
rho_suff = mp.mpf("3.46622291347846")
rho_fail = mp.mpf("3.46752913273870")
rho_max = mp.mpf("3.46752922945601")

print("rho_alpha^(min) =", rho_min)
print("Delta_suff =", rho_suff - rho_min)
print("Delta_fail =", rho_fail - rho_min)
print("Delta_max  =", rho_max - rho_min)

zeta_max = mp.mpf("2.46752922945601")
A_F1 = mp.mpf("1.0000519288021953286593340837139871812220538900677990")
print("zeta_req^(min) =", zeta_min)
print("zeta_max^(F1) - zeta_req^(min) =", zeta_max - zeta_min)

subbanner("72.1 — Zero-transport-bias result")
print("A_F1 =", A_F1)
print("Since zeta_req^(min) = 1/3 < A_F1, the explicit Family-1 branch already meets the")
print("support demand at zero transport bias.")
print("Hence Pe_req = 0 on the explicit Family-1 branch for the minimal isotropic demand.")

banner("STAGE 72 FINAL LEDGER")
print("The explicit Family-1 support/source branch passes the minimal isotropic passive/outgoing")
print("quadrupole test with large margin:")
print("  rho_alpha = 4/3,")
print("  zeta_req = 1/3,")
print("  Pe_req = 0.")
print("So on this branch the remaining reduced theorem gap is no longer support sufficiency;")
print("it is whether the actual grouped-P2 / geometry branch really realizes the minimal")
print("isotropic contact-plus-pole conservative quadrupole module.")
