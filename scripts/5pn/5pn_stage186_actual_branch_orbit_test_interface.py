#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage181_186_common import *

"""
5pn_stage186_actual_branch_orbit_test_interface.py

Stage 186 — exact test interface for future PDE branch data.

What this script does
---------------------
1. Packages the finite-orbit formulas into a single symbolic test interface for an
   actual moving-throat branch.
2. From a candidate microscopic state, computes the exact orbit-compatible
   dependent triple, the finite mismatch ratios, the quotient coordinates, and the
   restoration map.
3. Shows that all previous first-order branch-selection diagnostics are obtained by
   taking logarithms / linearizations of these exact finite objects.
"""

banner("STAGE 186 — EXACT TEST INTERFACE FOR FUTURE PDE BRANCH DATA")

syms = base_symbols()
inv = invariants_expr(syms)
orb = finite_orbit_dependent_triple(syms)

Keta = syms["Keta"]
muW = syms["muW"]
TU = syms["TU"]
Ctrs = syms["Ctrs"]
Cnts = syms["Cnts"]
epss = syms["epss"]
chi0s = syms["chi0s"]
Fstar = syms["Fstar"]

subbanner("I. Orbit-compatible dependent triple for a candidate free-coordinate point")
print("T_U^(orbit) =")
sp.pprint(orb["TU_orbit"])
print("K_eta^(orbit) =")
sp.pprint(orb["Keta_orbit"])
print("mu_W^(orbit) =")
sp.pprint(orb["muW_orbit"])

subbanner("II. Exact finite mismatch ratios for an actual candidate branch")
mT_actual = sp.simplify(TU / orb["TU_orbit"])
mK_actual = sp.simplify(Keta / orb["Keta_orbit"])
mMu_actual = sp.simplify(muW / orb["muW_orbit"])

print("m_T(actual) =")
sp.pprint(mT_actual)
print("m_K(actual) =")
sp.pprint(mK_actual)
print("m_mu(actual) =")
sp.pprint(mMu_actual)

subbanner("III. Exact quotient coordinates extracted from the actual candidate branch")
qtr_actual = sp.simplify(sp.log(inv["Ctr"] / Ctrs))
qnt_actual = sp.simplify(sp.log(inv["Cnt"] / Cnts))
qeta_actual = sp.simplify(sp.log(inv["eps_eta"] / epss))

print("q_tr(actual) =")
sp.pprint(qtr_actual)
print("q_nt(actual) =")
sp.pprint(qnt_actual)
print("q_eta(actual) =")
sp.pprint(qeta_actual)

expect_zero("q_tr(actual) - (1+chi0_*) ln m_T", qtr_actual - (1 + chi0s) * sp.log(mT_actual))
expect_zero("q_eta(actual) + ln m_K", qeta_actual + sp.log(mK_actual))
expect_zero(
    "q_nt(actual) - (ln m_mu - ln m_K - F_* ln m_T)",
    qnt_actual - (sp.log(mMu_actual) - sp.log(mK_actual) - Fstar * sp.log(mT_actual)),
)

subbanner("IV. Exact restoration map from actual branch data")
T_restore = sp.simplify(TU / mT_actual)
Keta_restore = sp.simplify(Keta / mK_actual)
mu_restore = sp.simplify(muW / mMu_actual)

expect_zero("restore T_U to orbit", T_restore / orb["TU_orbit"] - 1)
expect_zero("restore K_eta to orbit", Keta_restore / orb["Keta_orbit"] - 1)
expect_zero("restore mu_W to orbit", mu_restore / orb["muW_orbit"] - 1)

print("T_U^(restore) =")
sp.pprint(T_restore)
print("K_eta^(restore) =")
sp.pprint(Keta_restore)
print("mu_W^(restore) =")
sp.pprint(mu_restore)

subbanner("V. Practical verdict form")
print("Actual branch is orbit-locked iff any one of the following equivalent tests holds:")
print("  (a) m_T = m_K = m_mu = 1")
print("  (b) q_tr = q_nt = q_eta = 0")
print("  (c) (T_U, K_eta, mu_W) = (T_U^(orbit), K_eta^(orbit), mu_W^(orbit))")
print("  (d) restoration map acts trivially")

banner("STAGE 186 LEDGER")
print("1. The abstract orbit-lock theorem has now been turned into an explicit interface")
print("   for future PDE data: candidate microscopic branch values can be tested directly.")
print("2. The exact quantities to compute from the actual branch are")
print("      T_U, K_eta, mu_W")
print("   together with the five free microscopic coordinates and the invariant triple.")
print("3. From those, the script returns the exact orbit-compatible dependent triple, the")
print("   finite mismatch ratios, the quotient coordinates, and the restoration map.")
print("4. So the remaining theorem gap is no longer abstract. It is now a direct data test")
print("   on the actual moving-throat branch once the PDE numerics are available.")
