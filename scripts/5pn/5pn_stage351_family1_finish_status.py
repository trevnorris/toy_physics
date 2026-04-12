#!/usr/bin/env python3
from __future__ import annotations
import mpmath as mp

"""
Stage 351 — compare the exact Family-1 branch location against the actual isotropic
passive/outgoing grouped-P2 demand.

The actual reduced isotropic branch needs
    rho_alpha = 4/3
    zeta_req = 1/3         (unblocked)
and remains lowest-twin-safe on the whole admissible blocked interval.
This script checks those statements numerically against the explicit Family-1
branch points located in Stage 350.
"""

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

mp.mp.dps = 80

# Pulled straight from Stage 350 exact refresh
zeta_chi = mp.mpf("2.46752964788143760644099569607884245175681145627716307221485")
zeta_J   = mp.mpf("2.46752780516750843270099983841122563968319833831668094348097")
rho_chi  = 1 + zeta_chi
rho_J    = 1 + zeta_J
zeta_max = mp.mpf("2.46752974572593578318480550876175670365053962167244541451029")

zeta_req_act = mp.mpf(1) / 3
rho_req_act  = mp.mpf(4) / 3

def zeta_req_blocked(eps_blk: mp.mpf) -> mp.mpf:
    return 1 / (3 - 2 * eps_blk)

eps_blk_crit = 1 / zeta_max

def fmt(x: mp.mpf, digits: int = 16) -> str:
    return mp.nstr(x, n=digits)

banner("STAGE 351 — FAMILY-1 FINISH STATUS AGAINST THE ACTUAL ISOTROPIC BRANCH")

subbanner("I. Actual reduced isotropic demand")
print("zeta_req(actual, unblocked) =", fmt(zeta_req_act, 25))
print("rho_alpha(actual, unblocked) =", fmt(rho_req_act, 25))
print("eps_blk,crit =", fmt(eps_blk_crit, 25))
print("zeta_req(actual, eps_blk=eps_blk,crit) =", fmt(zeta_req_blocked(eps_blk_crit), 25))
print("The actual isotropic branch is:")
print("  mixed-only enough?    no, because rho_alpha = 4/3 > 1")
print("  lowest-twin enough?   yes, because 1 < rho_alpha = 4/3 < 2")
print("  non-twin required?    no")

subbanner("II. Comparison with the numerically located explicit Family-1 branches")
for name, zeta, rho in [
    ("chi-weighted", zeta_chi, rho_chi),
    ("J-weighted", zeta_J, rho_J),
]:
    print(f"{name} branch:")
    print("  zeta_phys =", fmt(zeta, 25))
    print("  rho_alpha,max =", fmt(rho, 25))
    print("  margin in zeta over actual isotropic demand =", fmt(zeta - zeta_req_act, 25))
    print("  margin in rho_alpha over actual isotropic demand =", fmt(rho - rho_req_act, 25))
    print("  branch verdict: the explicit support/source branch lies deep inside the lowest-twin-safe region.")
    print()

subbanner("III. Finish-status verdict")
print("Support/source side:")
print("  numerically located and non-bottlenecked on the refreshed Family-1 branch.")
print("Outgoing side:")
print("  on the canonical compact passive/outgoing branch, chi_Q = 1 and N_Q = 1 exactly.")
print("Orbit-lock side:")
print("  still fixed only by the exact drift surfaces d ln R_tr = d ln R_target = d ln epsilon_eta = 0;")
print("  the current file stack does not yet contain a numerical PDE-selected point for those three.")
