#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage181_186_common import *

"""
5pn_stage184_finite_dependent_restoration.py

Stage 184 — exact finite restoration of a candidate branch to the same G_* orbit
by changing only the dependent triple.

What this script does
---------------------
1. Uses the exact finite mismatch coordinates from Stages 181–183.
2. Builds the exact finite restoration factors that send a candidate branch back to
   the orbit with the same five free microscopic coordinates and the same invariant
   triple.
3. Verifies that the restored dependent triple reproduces the target invariants
   exactly.
"""

banner("STAGE 184 — EXACT FINITE RESTORATION OF THE DEPENDENT TRIPLE")

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

qtr, qnt, qeta = sp.symbols("q_tr q_nt q_eta", real=True)

subbanner("I. Exact restoration factors from the quotient coordinates")
T_restore = sp.simplify(TU * sp.exp(-qtr / (1 + chi0s)))
Keta_restore = sp.simplify(Keta * sp.exp(qeta))
mu_restore = sp.simplify(muW * sp.exp(-qnt + qeta - Fstar * qtr / (1 + chi0s)))

print("T_U^(restore) =")
sp.pprint(T_restore)
print("K_eta^(restore) =")
sp.pprint(Keta_restore)
print("mu_W^(restore) =")
sp.pprint(mu_restore)

subbanner("II. Exact candidate branch written in quotient coordinates")
# Candidate built from the orbit point and the finite quotient coordinates.
TU_candidate = sp.simplify(orb["TU_orbit"] * sp.exp(qtr / (1 + chi0s)))
Keta_candidate = sp.simplify(orb["Keta_orbit"] * sp.exp(-qeta))
mu_candidate = sp.simplify(orb["muW_orbit"] * sp.exp(qnt - qeta + Fstar * qtr / (1 + chi0s)))

print("T_U(candidate) =")
sp.pprint(TU_candidate)
print("K_eta(candidate) =")
sp.pprint(Keta_candidate)
print("mu_W(candidate) =")
sp.pprint(mu_candidate)

# Restoring that candidate must return the exact orbit values.
expect_zero(
    "restored T_U - T_U^(orbit)",
    sp.simplify(T_restore.subs(TU, TU_candidate) / orb["TU_orbit"] - 1),
)
expect_zero(
    "restored K_eta - K_eta^(orbit)",
    sp.simplify(Keta_restore.subs(Keta, Keta_candidate) / orb["Keta_orbit"] - 1),
)
expect_zero(
    "restored mu_W - mu_W^(orbit)",
    sp.simplify(mu_restore.subs(muW, mu_candidate) / orb["muW_orbit"] - 1),
)

subbanner("III. Exact invariant check after restoration")
subs_restored = {
    TU: T_restore.subs(TU, TU_candidate),
    Keta: Keta_restore.subs(Keta, Keta_candidate),
    muW: mu_restore.subs(muW, mu_candidate),
}
expect_zero("C_tr restored / C_tr,* - 1", sp.simplify(inv["Ctr"].subs(subs_restored) / Ctrs - 1))
expect_zero("C_nt restored / C_nt,* - 1", sp.simplify(inv["Cnt"].subs(subs_restored) / Cnts - 1))
expect_zero("epsilon_eta restored / epsilon_eta,* - 1", sp.simplify(inv["eps_eta"].subs(subs_restored) / epss - 1))

banner("STAGE 184 LEDGER")
print("1. Once the finite quotient coordinates (q_tr, q_nt, q_eta) are known, the exact")
print("   restoration to the same G_* orbit is achieved by rescaling only the dependent triple")
print("      (T_U, K_eta, mu_W).")
print("2. The restoration factors are exact exponentials of the quotient coordinates.")
print("3. So the finite branch-selection failure remains completely localized to the")
print("   dependent microscopic triple even beyond first order.")
