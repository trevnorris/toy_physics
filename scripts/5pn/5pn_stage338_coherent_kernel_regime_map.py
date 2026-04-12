#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage336_339_common import *

"""
Stage 338 — coherent-kernel regime map in direct continuum variables.

What this script does
---------------------
1. Inserts the coherent-kernel baseline
      M_mix = 8 Z_W (1+chi_0)^2 / [pi^2 (1-eps_eta)(1-eps)]
   into the selected-branch loading ratio.
2. Derives exact lower bounds on Z_W for mixed-only and lowest-twin sufficiency.
3. Derives exact lower bounds on the radiative demand scale Lambda for the same
   two regimes.
4. Shows that at fixed selected point, increasing Z_W lowers rho_alpha.
"""

banner("STAGE 338 — COHERENT-KERNEL REGIME MAP")

xi, delta, R = sp.symbols("xi delta R", positive=True, real=True)
eps_eta, eps, chi0, ZW, Lambda = sp.symbols("eps_eta eps chi_0 Z_W Lambda", positive=True, real=True)

G = G_tr(xi, delta, R)
Pi = Pi_tr(xi, delta, R)
Mmix = sp.simplify(8 * ZW * (1 + chi0)**2 / (sp.pi**2 * (1 - eps_eta) * (1 - eps)))
rho = sp.simplify(G / Mmix)

subbanner("I. Exact coherent-kernel loading ratio")
print("M_mix =")
sp.pprint(Mmix)
print("rho_alpha^(coh) = G_tr / M_mix =")
sp.pprint(rho)

subbanner("II. Exact Z_W thresholds")
ZW_mix = ZW_threshold(sp.Integer(1), xi, delta, R, eps_eta, eps, chi0)
ZW_twin = ZW_threshold(sp.Integer(2), xi, delta, R, eps_eta, eps, chi0)
print("Z_W^(mix-only req) =")
sp.pprint(ZW_mix)
print("Z_W^(twin req) =")
sp.pprint(ZW_twin)

expect_zero(
    "rho_alpha(ZW_mix) - 1",
    sp.simplify(rho.subs(ZW, ZW_mix) - 1),
)
expect_zero(
    "rho_alpha(ZW_twin) - 2",
    sp.simplify(rho.subs(ZW, ZW_twin) - 2),
)

subbanner("III. Exact Lambda thresholds from the branch product")
Lambda_mix = Lambda_threshold(sp.Integer(1), xi, delta, R, eps)
Lambda_twin = Lambda_threshold(sp.Integer(2), xi, delta, R, eps)
print("Lambda^(mix-only req) =")
sp.pprint(Lambda_mix)
print("Lambda^(twin req) =")
sp.pprint(Lambda_twin)

Cmix = c_mix(Lambda, eps)
expect_zero(
    "Pi_tr - C_mix at Lambda_mix",
    sp.simplify(Pi.subs(Lambda, Lambda_mix) - Cmix.subs(Lambda, Lambda_mix)),
)
expect_zero(
    "Pi_tr - 2 C_mix at Lambda_twin",
    sp.simplify(Pi.subs(Lambda, Lambda_twin) - 2 * Cmix.subs(Lambda, Lambda_twin)),
)

subbanner("IV. Monotonicity in wall-to-mixed overlap")
drho_dZW = sp.simplify(sp.diff(rho, ZW))
print("d rho_alpha / d Z_W =")
sp.pprint(drho_dZW)

expect_zero(
    "d rho_alpha / d Z_W + rho_alpha / Z_W",
    sp.simplify(drho_dZW + rho / ZW),
)

print("So at fixed selected branch point, increasing Z_W lowers rho_alpha and makes")
print("mixed-only / lowest-twin support success easier.")
