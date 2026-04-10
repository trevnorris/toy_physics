#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage181_186_common import *

"""
5pn_stage181_finite_orbit_dependent_law.py

Stage 181 — exact finite single-orbit law for the dependent microscopic triple
(T_U, K_eta, mu_W).

What this script does
---------------------
1. Starts from the exact Stage-168 direct microscopic monomials
      C_tr,* , C_nt,* , epsilon_eta.
2. Solves them exactly for the dependent coordinates
      (T_U, K_eta, mu_W)
   in terms of the five free microscopic coordinates
      (lambda_W, c_etaU, gamma, K_U, K_W)
   and the fixed invariant triple.
3. Verifies by exact substitution that the reconstructed dependent triple lies on
   the same finite similarity orbit.
"""

banner("STAGE 181 — EXACT FINITE SINGLE-ORBIT LAW FOR THE DEPENDENT TRIPLE")

syms = base_symbols()
inv = invariants_expr(syms)
orb = finite_orbit_dependent_triple(syms)

lamW = syms["lamW"]
c_etaU = syms["c_etaU"]
gamma = syms["gamma"]
KU = syms["KU"]
Keta = syms["Keta"]
KW = syms["KW"]
muW = syms["muW"]
TU = syms["TU"]
L = syms["L"]
sigma = syms["sigma"]
chi0s = syms["chi0s"]
deltaUs = syms["deltaUs"]
Estar = syms["Estar"]
Fstar = syms["Fstar"]
Ctrs = syms["Ctrs"]
Cnts = syms["Cnts"]
epss = syms["epss"]

subbanner("I. Exact monomials carried from Stage 168")
print("C_tr,* =")
sp.pprint(inv["Ctr"])
print("C_nt,* =")
sp.pprint(inv["Cnt"])
print("epsilon_eta =")
sp.pprint(inv["eps_eta"])

subbanner("II. Exact finite orbit solve for (T_U, K_eta, mu_W)")
print("K_eta^(orbit) =")
sp.pprint(orb["Keta_orbit"])
print("T_U^(orbit) =")
sp.pprint(orb["TU_orbit"])
print("mu_W^(orbit) =")
sp.pprint(orb["muW_orbit"])

subbanner("III. Exact substitution check")
subs_orbit = {
    Keta: orb["Keta_orbit"],
    TU: orb["TU_orbit"],
    muW: orb["muW_orbit"],
}
expect_zero("C_tr orbit / C_tr,* - 1", sp.simplify(inv["Ctr"].subs(subs_orbit) / Ctrs - 1))
expect_zero("C_nt orbit / C_nt,* - 1", sp.simplify(inv["Cnt"].subs(subs_orbit) / Cnts - 1))
expect_zero("epsilon_eta orbit / epsilon_eta,* - 1", sp.simplify(inv["eps_eta"].subs(subs_orbit) / epss - 1))

subbanner("IV. Exact finite-law summary")
print("The Stage-170 finite orbit through a given free microscopic point")
print("(lambda_W, c_etaU, gamma, K_U, K_W)")
print("is obtained by setting the dependent triple to the exact values above.")
print("No linearization is used anywhere in this solve.")

banner("STAGE 181 LEDGER")
print("1. The finite quotient theorem is now explicit on the dependent triple.")
print("2. Fixing the invariant triple (C_tr,* , C_nt,* , epsilon_eta,*) and the five")
print("   free microscopic coordinates uniquely fixes")
print("      T_U^(orbit), K_eta^(orbit), mu_W^(orbit).")
print("3. So the finite similarity orbit is no longer abstract: the dependent")
print("   coordinates are exact algebraic functions of the free coordinates and")
print("   the invariant triple.")
