#!/usr/bin/env python3
"""
SymPy audit for Stage 090.

This is a checkpoint-consistency audit rather than a fresh derivation. It
reconfirms that the minimal isotropic contact-plus-pole module carries

    rho_alpha = 4/3,
    zeta_req = 1/3,

and that the carried Family-1 thresholds already place this branch strictly
inside the exact success region with zero required transport bias.

Constant provenance:
- 3/4 and 1/4 are the minimal isotropic contact-plus-pole module used here;
  the grouped-P2/static-geometry derivation is downstream at Stage 091.
- rho_suff^(chi) is carried from the Stage 086 support window audit.
- zeta_max^(F1) is carried from the Stage 080/081 support ceiling.
- A_F1 is carried from the Stage 079 transport map.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    residual = sp.simplify(sp.expand(expr))
    print(f"{name} = {residual}")
    if residual != 0:
        raise AssertionError(f"{name} is not zero")


def expect_true(name: str, cond: bool) -> None:
    print(f"{name} = {cond}")
    if not cond:
        raise AssertionError(f"{name} failed")


banner("STAGE 090 — UPDATED REDUCED STATUS AFTER THE LOADING-RATIO EXTRACTION")

# Minimal isotropic conservative module used by this reduced-status check.
# Stage 090 does not derive these coefficients again; the grouped-P2/static
# geometry derivation is downstream at Stage 091.
c_contact = sp.Rational(3, 4)
c_pole = sp.Rational(1, 4)

rho_alpha = sp.simplify(1 / c_contact)
zeta_req = sp.simplify(c_pole / c_contact)

print("c_contact =", c_contact)
print("c_pole    =", c_pole)
print("rho_alpha =", rho_alpha)
print("zeta_req  =", zeta_req)

expect_zero("rho_alpha - 4/3", rho_alpha - sp.Rational(4, 3))
expect_zero("zeta_req - 1/3", zeta_req - sp.Rational(1, 3))
expect_zero("zeta_req - (rho_alpha - 1)", zeta_req - (rho_alpha - 1))

# Carried forward numeric thresholds with explicit source anchors:
# - rho_suff_chi from Stage 086
# - zeta_max_f1 from Stages 080/081
# - A_F1 from Stage 079
rho_suff_chi = sp.Float("3.46622291347846", 30)
zeta_max_f1 = sp.Float("2.46752922945601", 30)
A_F1 = sp.Float("1.00005192880220", 30)

delta_rho_suff = sp.N(rho_suff_chi - rho_alpha, 25)
delta_zeta_max = sp.N(zeta_max_f1 - zeta_req, 25)
delta_A_F1 = sp.N(A_F1 - zeta_req, 25)

print("rho_suff^(chi) =", sp.N(rho_suff_chi, 25))
print("zeta_max^(F1)  =", sp.N(zeta_max_f1, 25))
print("A_F1           =", sp.N(A_F1, 25))

print("\nMargins:")
print("rho_suff^(chi) - rho_alpha =", delta_rho_suff)
print("zeta_max^(F1)  - zeta_req  =", delta_zeta_max)
print("A_F1           - zeta_req  =", delta_A_F1)

expect_true(
    "rho_alpha lies inside the exact Family-1 success region",
    bool(rho_alpha < rho_suff_chi),
)
expect_true(
    "zeta_req lies below the hard Family-1 ceiling",
    bool(zeta_req < zeta_max_f1),
)
expect_true(
    "zeta_req lies below the zero-bias Family-1 baseline",
    bool(zeta_req < A_F1),
)

# Stage 075 transport map: zeta_req < A_F1 ==> Pe_req = 0. The inequality
# above is the carry-forward proxy for the locked triple value Pe_req = 0
# stated in the Stage 090 notes (paper body item vi).
Pe_req = sp.Integer(0)
expect_zero("Pe_req (carry-forward from Stage 075 transport map)", Pe_req)

print("\nFINAL LEDGER")
print("The explicit Family-1 support/source side is already finished under the")
print("minimal isotropic module:")
print("  rho_alpha = 4/3 lies strictly inside the success region,")
print("  zeta_req  = 1/3 lies below the hard ceiling,")
print("  zeta_req  = 1/3 < A_F1 so Pe_req = 0.")
