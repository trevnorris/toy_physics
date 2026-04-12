#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage181_186_common import *

"""
5pn_stage182_finite_candidate_factorization.py

Stage 182 — exact finite factorization of a candidate branch into orbit data and
three dependent-coordinate mismatch ratios.

What this script does
---------------------
1. Uses the exact Stage-181 orbit solve.
2. Factors an arbitrary candidate dependent triple as
      T_U = m_T T_U^(orbit),
      K_eta = m_K K_eta^(orbit),
      mu_W = m_mu mu_W^(orbit).
3. Derives the exact invariant-ratio laws in terms of (m_T, m_K, m_mu).
"""

banner("STAGE 182 — EXACT FINITE FACTORIZATION OF A CANDIDATE BRANCH")

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

mT, mK, mMu = sp.symbols("m_T m_K m_mu", positive=True, real=True)

subbanner("I. Exact dependent-coordinate factorization")
print("T_U   = m_T  * T_U^(orbit)")
print("K_eta = m_K  * K_eta^(orbit)")
print("mu_W  = m_mu * mu_W^(orbit)")

subs_candidate = {
    TU: mT * orb["TU_orbit"],
    Keta: mK * orb["Keta_orbit"],
    muW: mMu * orb["muW_orbit"],
}

subbanner("II. Exact invariant-ratio laws")
Ctr_ratio = sp.simplify(inv["Ctr"].subs(subs_candidate) / Ctrs)
Cnt_ratio = sp.simplify(inv["Cnt"].subs(subs_candidate) / Cnts)
eps_ratio = sp.simplify(inv["eps_eta"].subs(subs_candidate) / epss)

print("C_tr / C_tr,* =", Ctr_ratio)
print("C_nt / C_nt,* =", Cnt_ratio)
print("epsilon_eta / epsilon_eta,* =", eps_ratio)

expect_zero("C_tr ratio - m_T^(1+chi0_*)", Ctr_ratio - mT ** (1 + chi0s))
expect_zero("C_nt ratio - m_mu/(m_K m_T^F_*)", Cnt_ratio - mMu / (mK * mT ** Fstar))
expect_zero("epsilon ratio - 1/m_K", eps_ratio - 1 / mK)

subbanner("III. Orbit-lock criterion in finite mismatch variables")
print("Orbit lock  <=>  m_T = m_K = m_mu = 1")
print("because only then are all three invariant ratios equal to 1 exactly.")

banner("STAGE 182 LEDGER")
print("1. Any candidate branch with the same five free microscopic coordinates as a given")
print("   orbit point is described completely by three positive mismatch ratios")
print("      (m_T, m_K, m_mu).")
print("2. The exact quotient motion is then")
print("      C_tr / C_tr,* = m_T^(1+chi_0,*),")
print("      C_nt / C_nt,* = m_mu / (m_K m_T^{F_*}),")
print("      epsilon_eta / epsilon_eta,* = 1 / m_K.")
print("3. So the finite branch-selection problem has been reduced to three exact")
print("   multiplicative mismatch coordinates on the dependent triple.")
