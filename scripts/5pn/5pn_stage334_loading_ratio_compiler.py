#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage333_335_common import *

"""
Stage 334 — exact loading-ratio compiler and canonical isotropic precursor.

What this script does
---------------------
1. Rewrites the support demand in the pure loading-ratio variable rho_alpha.
2. Derives the exact contact-plus-pole inverse formulas
      c0 = 1/rho_alpha, c1 = (rho_alpha - 1)/rho_alpha.
3. Inserts the minimal isotropic conservative precursor
      c0 = 3/4, c1 = 1/4
   and proves
      rho_alpha = 4/3,
      zeta_req = 1/(3 - 2 eps),
      Pi_tr = (4/3) C_mix.
4. Proves that for every blocked branch with 0 <= eps < 1 the canonical isotropic
   branch stays inside the symmetric-lowest-twin regime.
"""

banner("STAGE 334 — LOADING-RATIO COMPILER AND CANONICAL ISOTROPIC PRECURSOR")

rho_alpha, eps = sp.symbols("rho_alpha eps", positive=True, real=True)
c0, c1 = sp.symbols("c0 c1", positive=True, real=True)

subbanner("I. Exact loading-ratio representation")
coeffs = precursor_coeffs_from_ratio(rho_alpha)
inv = ratio_from_precursor(c0, c1)

print("c0(rho_alpha) =")
sp.pprint(coeffs["c0"])
print("c1(rho_alpha) =")
sp.pprint(coeffs["c1"])
print("Inverse formulas:")
sp.pprint(sp.Matrix([inv["rho_alpha_from_c0"], inv["rho_alpha_from_c1"], inv["zeta_req"]]))

expect_zero("c0 + c1 - 1", coeffs["c0"] + coeffs["c1"] - 1)
expect_zero("rho from c0", inv["rho_alpha_from_c0"].subs({c0: coeffs["c0"]}) - rho_alpha)
expect_zero("rho from c1", inv["rho_alpha_from_c1"].subs({c1: coeffs["c1"]}) - rho_alpha)
expect_zero("zeta = rho - 1", coeffs["c1"]/coeffs["c0"] - (rho_alpha - 1))

subbanner("II. Support demand in the pure loading-ratio variable")
zeta_ratio = zeta_req_from_ratio(rho_alpha, eps)
print("zeta_req(rho_alpha, eps) =")
sp.pprint(zeta_ratio)

subbanner("III. Minimal isotropic conservative precursor")
rho_min = sp.simplify(inv["rho_alpha_from_c0"].subs({c0: sp.Rational(3,4)}))
zeta_min = sp.simplify(zeta_ratio.subs(rho_alpha, rho_min))
print("From c0 = 3/4, c1 = 1/4:")
print("rho_alpha^(min) =")
sp.pprint(rho_min)
print("zeta_req^(min)(eps) =")
sp.pprint(zeta_min)

expect_zero("rho_min - 4/3", rho_min - sp.Rational(4,3))
expect_zero("zeta_min - 1/(3 - 2 eps)", zeta_min - 1/(3 - 2*eps))

Lambda = sp.symbols("Lambda", positive=True, real=True)
Cmix = c_mix(Lambda, eps)
Pi_min = sp.simplify(rho_min * Cmix)
print("Pi_tr^(min) =")
sp.pprint(Pi_min)
expect_zero("Pi_min - (4/3) C_mix", Pi_min - sp.Rational(4,3) * Cmix)

subbanner("IV. Canonical isotropic branch is always twin-safe on the blocked branch")
expect_zero(
    "1 - zeta_min - 2(1-eps)/(3-2eps)",
    1 - zeta_min - sp.simplify(2*(1-eps)/(3-2*eps)),
)
print("For every blocked branch with 0 <= eps < 1, zeta_req^(min) < 1.")
print("Therefore the canonical isotropic passive/outgoing branch always lies in the")
print("symmetric-lowest-twin regime before any non-twin asymmetry is needed.")
