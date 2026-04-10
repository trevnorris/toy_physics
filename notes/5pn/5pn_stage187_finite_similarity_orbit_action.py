#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage187_192_common import *

"""
Stage 187 — exact finite similarity-orbit action on the full microscopic state.

What this script does
---------------------
1. Writes the exact five-parameter finite orbit action G_* on all eight microscopic
   coordinates.
2. Proves that the three quotient monomials (C_tr, C_nt, epsilon_eta) are preserved
   exactly by that action, not only infinitesimally.
3. Shows that the dependent triple (T_U, K_eta, mu_W) is slaved algebraically to the
   five free coordinates along any fixed orbit.
"""

banner("STAGE 187 — EXACT FINITE SIMILARITY-ORBIT ACTION")

syms = base_symbols()
pars = orbit_generators()
inv0 = invariants_expr(syms)
syms1 = orbit_transform(syms, pars)
inv1 = invariants_expr(syms1)
exps = orbit_exponents(syms, pars)

subbanner("I. Exact finite orbit action on microscopic coordinates")
for key in ("lamW", "c_etaU", "gamma", "KU", "Keta", "KW", "muW", "TU"):
    print(f"{key}' =")
    sp.pprint(syms1[key])
    print()

subbanner("II. Exact orbit exponents")
for key in ("lamW", "c_etaU", "gamma", "KU", "Keta", "KW", "muW", "TU"):
    print(f"Delta ln {key} =")
    sp.pprint(sp.simplify(exps[key]))
    print()

subbanner("III. Exact monomial preservation")
expect_zero("ln(C_tr'/C_tr)", sp.expand_log(sp.log(inv1["Ctr"] / inv0["Ctr"]), force=True))
expect_zero("ln(C_nt'/C_nt)", sp.expand_log(sp.log(inv1["Cnt"] / inv0["Cnt"]), force=True))
expect_zero("ln(epsilon_eta'/epsilon_eta)", sp.expand_log(sp.log(inv1["eps_eta"] / inv0["eps_eta"]), force=True))

subbanner("IV. Exact dependent-triple coevolution along a fixed orbit")
print("K_eta' / K_eta =")
sp.pprint(sp.simplify(syms1["Keta"] / syms["Keta"]))
print("T_U' / T_U =")
sp.pprint(sp.simplify(syms1["TU"] / syms["TU"]))
print("mu_W' / mu_W =")
sp.pprint(sp.simplify(syms1["muW"] / syms["muW"]))

subbanner("V. Exact orbit graph interpretation")
print("The finite similarity orbit G_* is a graph over the five free coordinates")
print("(lambda_W, c_etaU, gamma, K_U, K_W), with the dependent triple")
print("(K_eta, T_U, mu_W) transported multiplicatively by the exact exponents above.")
print("So once a reference point on an orbit is fixed, every other point on that orbit")
print("is determined by five free multiplicative scales.")

banner("STAGE 187 LEDGER")
print("1. The exact weak-axisymmetric similarity orbit is now an explicit finite action")
print("   on the full eight-dimensional microscopic state.")
print("2. The three quotient monomials are preserved exactly under this action.")
print("3. The dependent triple (T_U, K_eta, mu_W) is slaved exactly to the five free")
print("   microscopic coordinates along any fixed orbit.")
