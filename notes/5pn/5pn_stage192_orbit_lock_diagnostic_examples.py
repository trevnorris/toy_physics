#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage187_192_common import *

"""
Stage 192 — concrete symbolic diagnostic examples for orbit lock and failure channels.

What this script does
---------------------
1. Gives exact symbolic examples of an orbit-locked candidate and simple off-orbit
   candidates.
2. Shows which quotient coordinate turns on for pure T_U, pure K_eta, and pure mu_W
   dependent residual mismatches.
3. Consolidates the practical diagnostic signatures of each failure channel.
"""

banner("STAGE 192 — ORBIT-LOCK DIAGNOSTIC EXAMPLES")

syms = base_symbols()
rat = free_ratio_symbols()
laws = orbit_ratio_laws(syms, rat)

mT, mK, mMu = sp.symbols("m_T m_K m_mu", positive=True, real=True)

# Build generic residual chart.
qtr = sp.simplify((1 + syms["chi0s"]) * sp.log(mT))
qeta = sp.simplify(-sp.log(mK))
qnt = sp.simplify(sp.log(mMu) - sp.log(mK) - syms["Fstar"] * sp.log(mT))

subbanner("I. Orbit-locked candidate")
expect_zero("q_tr(lock)", qtr.subs({mT: 1}))
expect_zero("q_eta(lock)", qeta.subs({mK: 1}))
expect_zero("q_nt(lock)", qnt.subs({mT: 1, mK: 1, mMu: 1}))
print("If m_T = m_K = m_mu = 1, then q_tr = q_eta = q_nt = 0 exactly.")

subbanner("II. Pure T_U mismatch channel")
zT = sp.symbols("z_T", positive=True, real=True)
print("q_tr(T-only) =")
sp.pprint(sp.simplify(qtr.subs({mT: zT})))
print("q_eta(T-only) =")
sp.pprint(sp.simplify(qeta.subs({mK: 1})))
print("q_nt(T-only) =")
sp.pprint(sp.simplify(qnt.subs({mT: zT, mK: 1, mMu: 1})))

subbanner("III. Pure K_eta mismatch channel")
zK = sp.symbols("z_K", positive=True, real=True)
print("q_tr(K-only) =")
sp.pprint(sp.simplify(qtr.subs({mT: 1})))
print("q_eta(K-only) =")
sp.pprint(sp.simplify(qeta.subs({mK: zK})))
print("q_nt(K-only) =")
sp.pprint(sp.simplify(qnt.subs({mT: 1, mK: zK, mMu: 1})))

subbanner("IV. Pure mu_W mismatch channel")
zM = sp.symbols("z_M", positive=True, real=True)
print("q_tr(mu-only) =")
sp.pprint(sp.simplify(qtr.subs({mT: 1})))
print("q_eta(mu-only) =")
sp.pprint(sp.simplify(qeta.subs({mK: 1})))
print("q_nt(mu-only) =")
sp.pprint(sp.simplify(qnt.subs({mT: 1, mK: 1, mMu: zM})))

subbanner("V. Practical signatures")
print("- Pure T_U residual mismatch turns on q_tr and, through F_*, also q_nt.")
print("- Pure K_eta residual mismatch turns on q_eta and also q_nt.")
print("- Pure mu_W residual mismatch turns on q_nt only.")
print("So the three quotient coordinates diagnose the dependent triple cleanly.")

banner("STAGE 192 LEDGER")
print("1. The exact residual mismatch ratios provide a direct physical interpretation of")
print("   the three quotient coordinates.")
print("2. Each component of the dependent triple has a distinct diagnostic signature.")
print("3. This makes the future PDE verdict sharper: once the actual branch is computed,")
print("   the pattern of (q_tr, q_eta, q_nt) immediately identifies which dependent")
print("   coevolution law failed, if any.")
