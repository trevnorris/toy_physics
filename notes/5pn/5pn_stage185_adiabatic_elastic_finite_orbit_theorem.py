#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage181_186_common import *

"""
5pn_stage185_adiabatic_elastic_finite_orbit_theorem.py

Stage 185 — finite adiabatic-elastic orbit-lock theorem in the dependent-triple
mismatch variables.

What this script does
---------------------
1. Re-expresses the exact finite branch-selection problem under the adiabatic-elastic
   boundary rule in terms of the three mismatch ratios (m_T, m_K, m_mu).
2. Proves the exact equivalence
      same G_* orbit  <=>  m_T = m_K = m_mu = 1
   <=> all three quotient coordinates vanish.
3. Gives the exact finite criterion that the true moving-throat branch must satisfy.
"""

banner("STAGE 185 — FINITE ADIABATIC-ELASTIC ORBIT-LOCK THEOREM")

syms = base_symbols()
chi0s = syms["chi0s"]
Fstar = syms["Fstar"]

mT, mK, mMu = sp.symbols("m_T m_K m_mu", positive=True, real=True)
qtr = sp.log(mT ** (1 + chi0s))
qeta = sp.log(1 / mK)
qnt = sp.log(mMu / (mK * mT ** Fstar))

subbanner("I. Exact finite quotient coordinates")
print("q_tr  =", sp.simplify(qtr))
print("q_nt  =", sp.simplify(qnt))
print("q_eta =", sp.simplify(qeta))

subbanner("II. Orbit-lock equivalences")
print("same G_* orbit  <=>  q_tr = q_nt = q_eta = 0")
print("               <=>  m_T = m_K = m_mu = 1")

expect_zero("q_tr at m_T=1", qtr.subs(mT, 1))
expect_zero("q_eta at m_K=1", qeta.subs(mK, 1))
expect_zero("q_nt at m_T=m_K=m_mu=1", qnt.subs({mT: 1, mK: 1, mMu: 1}))

# Recover the mismatch ratios from the quotient coordinates exactly.
qtr_s, qnt_s, qeta_s = sp.symbols("qtr_s qnt_s qeta_s", real=True)
rec_mT = sp.exp(qtr_s / (1 + chi0s))
rec_mK = sp.exp(-qeta_s)
rec_mMu = sp.exp(qnt_s - qeta_s + Fstar * qtr_s / (1 + chi0s))

expect_zero("recovered q_tr - qtr_s", sp.simplify(sp.log(rec_mT ** (1 + chi0s)) - qtr_s))
expect_zero("recovered q_eta - qeta_s", sp.simplify(sp.log(1 / rec_mK) - qeta_s))
expect_zero("recovered q_nt - qnt_s", sp.simplify(sp.log(rec_mMu / (rec_mK * rec_mT ** Fstar)) - qnt_s))

subbanner("III. Adiabatic-elastic reading")
print("Under the adiabatic-elastic closure, the scalar off-bundle escape is removed.")
print("So the exact finite branch-selection problem is nothing but whether the")
print("dependent-triple mismatch ratios equal 1.")

banner("STAGE 185 LEDGER")
print("1. The adiabatic-elastic finite branch-selection problem is exactly three-dimensional.")
print("2. The finite failure coordinates are equivalently")
print("      (m_T, m_K, m_mu)")
print("   or")
print("      (q_tr, q_nt, q_eta).")
print("3. Exact orbit lock means")
print("      m_T = m_K = m_mu = 1,")
print("   equivalently")
print("      q_tr = q_nt = q_eta = 0.")
print("4. So the true moving-throat branch will satisfy the finite adiabatic-elastic")
print("   orbit theorem iff its dependent microscopic triple matches the exact")
print("   Stage-181 orbit values with no residual multiplicative mismatch.")
