#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage181_186_common import *

"""
5pn_stage183_exact_log_mismatch_coordinates.py

Stage 183 — exact logarithmic quotient coordinates of the finite dependent-triple
mismatch ratios.

What this script does
---------------------
1. Takes the exact finite mismatch laws from Stage 182.
2. Passes to logarithmic mismatch coordinates
      tau = ln m_T, kappa = ln m_K, mu = ln m_mu.
3. Proves that the exact quotient drifts are linear in those logarithmic mismatch
   coordinates, i.e. the Stage-179 formulas are not merely infinitesimal
   approximations but the exact log-chart of the finite mismatch ratios.
"""

banner("STAGE 183 — EXACT LOGARITHMIC QUOTIENT COORDINATES")

syms = base_symbols()
chi0s = syms["chi0s"]
Fstar = syms["Fstar"]

mT, mK, mMu = sp.symbols("m_T m_K m_mu", positive=True, real=True)

tau, kappa, mu = sp.symbols("tau kappa mu", real=True)

subbanner("I. Exact finite invariant ratios")
Ctr_ratio = mT ** (1 + chi0s)
Cnt_ratio = mMu / (mK * mT ** Fstar)
eps_ratio = 1 / mK

print("C_tr / C_tr,* =", Ctr_ratio)
print("C_nt / C_nt,* =", Cnt_ratio)
print("epsilon_eta / epsilon_eta,* =", eps_ratio)

subbanner("II. Exact logarithmic mismatch coordinates")
subs_logs = {mT: sp.exp(tau), mK: sp.exp(kappa), mMu: sp.exp(mu)}

qtr = sp.simplify(sp.log(Ctr_ratio.subs(subs_logs)))
qnt = sp.simplify(sp.log(Cnt_ratio.subs(subs_logs)))
qeta = sp.simplify(sp.log(eps_ratio.subs(subs_logs)))

print("q_tr  := ln(C_tr/C_tr,*) =", qtr)
print("q_nt  := ln(C_nt/C_nt,*) =", qnt)
print("q_eta := ln(epsilon_eta/epsilon_eta,*) =", qeta)

expect_zero("q_tr - (1+chi0_*) tau", qtr - (1 + chi0s) * tau)
expect_zero("q_nt - (mu - kappa - F_* tau)", qnt - (mu - kappa - Fstar * tau))
expect_zero("q_eta + kappa", qeta + kappa)

subbanner("III. Exact identification with the Stage-179 mismatch formulas")
print("tau   = ln(T_U / T_U^(orbit))")
print("kappa = ln(K_eta / K_eta^(orbit))")
print("mu    = ln(mu_W / mu_W^(orbit))")
print("Then the quotient drifts are exactly")
print("  d ln C_tr,*  = (1+chi_0,*) tau,")
print("  d ln C_nt,*  = mu - kappa - F_* tau,")
print("  d ln epsilon_eta = -kappa.")

banner("STAGE 183 LEDGER")
print("1. The Stage-179 mismatch formulas are the exact logarithmic chart of the")
print("   finite dependent-triple mismatch ratios, not merely a first-order guess.")
print("2. So the three exact quotient coordinates may be measured either as")
print("   finite invariant ratios or as linear functions of the log-mismatch triple")
print("      (tau, kappa, mu).")
print("3. This is the cleanest finite-to-infinitesimal bridge in the branch-selection")
print("   problem reached so far.")
