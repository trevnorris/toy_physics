#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage187_192_common import *

"""
Stage 190 — exact dependent residual coordinates relative to reference-orbit transport.

What this script does
---------------------
1. Factors a general candidate branch into exact reference-orbit transport plus three
   dependent residual mismatch ratios.
2. Proves that the three invariant ratios depend only on those residual mismatch ratios.
3. Rewrites the exact quotient coordinates as the logarithmic chart of those residuals.
"""

banner("STAGE 190 — EXACT DEPENDENT RESIDUAL COORDINATES")

syms = base_symbols()
rat = free_ratio_symbols()
laws = orbit_ratio_laws(syms, rat)
mis = mismatch_ratio_symbols()

RT, RK, RM = sp.symbols("R_T_actual R_K_actual R_mu_actual", positive=True, real=True)

# Total actual ratios are orbit transport times residual mismatch.
RT_total = sp.simplify(laws["RTU"] * mis["mT"])
RK_total = sp.simplify(laws["RKeta"] * mis["mK"])
RM_total = sp.simplify(laws["Rmu"] * mis["mMu"])

subbanner("I. Exact factorization of actual branch ratios")
print("R_T(actual) =")
sp.pprint(RT_total)
print("R_Keta(actual) =")
sp.pprint(RK_total)
print("R_muW(actual) =")
sp.pprint(RM_total)

# Build actual/reference ratio state.
ref = syms
act = dict(syms)
act["lamW"] = sp.simplify(rat["Rlam"] * ref["lamW"])
act["c_etaU"] = sp.simplify(rat["Rc"] * ref["c_etaU"])
act["gamma"] = sp.simplify(rat["Rgam"] * ref["gamma"])
act["KU"] = sp.simplify(rat["RU"] * ref["KU"])
act["KW"] = sp.simplify(rat["RW"] * ref["KW"])
act["Keta"] = sp.simplify(RK_total * ref["Keta"])
act["TU"] = sp.simplify(RT_total * ref["TU"])
act["muW"] = sp.simplify(RM_total * ref["muW"])

inv_ref = invariants_expr(ref)
inv_act = invariants_expr(act)

subbanner("II. Exact invariant ratios depend only on the residual triple")
expect_zero(
    "ln[(C_tr(actual)/C_tr(ref)) / m_T^(1+chi0_*)]",
    sp.expand_log(sp.log((inv_act["Ctr"] / inv_ref["Ctr"]) / (mis["mT"] ** (1 + syms["chi0s"]))), force=True),
)
expect_zero(
    "ln[(epsilon_eta(actual)/epsilon_eta(ref)) * m_K]",
    sp.expand_log(sp.log((inv_act["eps_eta"] / inv_ref["eps_eta"]) * mis["mK"]), force=True),
)
expect_zero(
    "ln[(C_nt(actual)/C_nt(ref)) / (m_mu/(m_K m_T^F_*))]",
    sp.expand_log(sp.log((inv_act["Cnt"] / inv_ref["Cnt"]) / (mis["mMu"] / (mis["mK"] * mis["mT"] ** syms["Fstar"]))), force=True),
)

subbanner("III. Exact logarithmic chart of the residual mismatch triple")
qtr = sp.simplify(sp.log(inv_act["Ctr"] / inv_ref["Ctr"]))
qeta = sp.simplify(sp.log(inv_act["eps_eta"] / inv_ref["eps_eta"]))
qnt = sp.simplify(sp.log(inv_act["Cnt"] / inv_ref["Cnt"]))

expect_zero("q_tr - (1+chi0_*) ln m_T", qtr - (1 + syms["chi0s"]) * sp.log(mis["mT"]))
expect_zero("q_eta + ln m_K", qeta + sp.log(mis["mK"]))
expect_zero("q_nt - (ln m_mu - ln m_K - F_* ln m_T)", qnt - (sp.log(mis["mMu"]) - sp.log(mis["mK"]) - syms["Fstar"] * sp.log(mis["mT"])))

print("q_tr =")
sp.pprint(qtr)
print("q_eta =")
sp.pprint(qeta)
print("q_nt =")
sp.pprint(qnt)

banner("STAGE 190 LEDGER")
print("1. Any actual branch can be factored exactly into reference-orbit transport plus")
print("   three dependent residual mismatch ratios (m_T, m_K, m_mu).")
print("2. The three invariant ratios depend only on that residual triple.")
print("3. So the exact quotient coordinates are the logarithmic chart of the dependent")
print("   residual mismatches, not of the free-coordinate transport along the orbit.")
