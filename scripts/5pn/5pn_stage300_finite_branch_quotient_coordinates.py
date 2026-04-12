#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage300_302_common import *

"""
Stage 300 — exact finite branch quotient coordinates.

This stage turns the Stage-299 microscopic monomials into exact finite quotient
coordinates between an actual branch state and a reference branch state, and
verifies that their first linearization is precisely the Stage-299 monomial-drift
vector.
"""

banner("STAGE 300 — FINITE BRANCH QUOTIENT COORDINATES")

# Actual microscopic variables.
lamW, cetaU, gam, KU, Keta, KW, muW, TU = sp.symbols(
    "lambda_W c_etaU gamma K_U K_etaeff K_Weff mu_W T_U", positive=True, real=True
)
# Reference branch variables.
lamWs, cetaUs, gams, KUs, Ketas, KWs, muWs, TUs = sp.symbols(
    "lambda_Ws c_etaUs gamma_s K_Us K_etaeffs K_Weffs mu_Ws T_Us", positive=True, real=True
)

L, sigma = sp.symbols("L sigma", positive=True, real=True)
chi0s, deltaUs, epsWs, epss = sp.symbols("chi0_star deltaU_star epsilonW_star epsilon_star", positive=True, real=True)

subbanner("I. Direct microscopic monomials on the actual and reference branches")
actual = microscopic_monomials(lamW, cetaU, gam, KU, Keta, KW, muW, TU, L, sigma, chi0s, deltaUs, epsWs, epss)
ref = microscopic_monomials(lamWs, cetaUs, gams, KUs, Ketas, KWs, muWs, TUs, L, sigma, chi0s, deltaUs, epsWs, epss)

print("C_tr,* (actual) =")
sp.pprint(actual["Ctr"])
print("C_nt,* (actual) =")
sp.pprint(actual["Cnt"])
print("epsilon_eta (actual) =")
sp.pprint(actual["epsEta"])

print("C_tr,* (reference) =")
sp.pprint(ref["Ctr"])
print("C_nt,* (reference) =")
sp.pprint(ref["Cnt"])
print("epsilon_eta (reference) =")
sp.pprint(ref["epsEta"])

subbanner("II. Exact finite quotient coordinates on the microscopic quotient space")
Qtr = sp.simplify(sp.log(actual["Ctr"] / ref["Ctr"]))
Qnt = sp.simplify(sp.log(actual["Cnt"] / ref["Cnt"]))
Qeta = sp.simplify(sp.log(actual["epsEta"] / ref["epsEta"]))

print("Q_tr = ln(C_tr/C_tr,ref) =")
sp.pprint(Qtr)
print("Q_nt = ln(C_nt/C_nt,ref) =")
sp.pprint(Qnt)
print("Q_eta = ln(epsilon_eta/epsilon_eta,ref) =")
sp.pprint(Qeta)

# Identity check: actual = reference implies zero quotient coordinates.
subs_equal = {
    lamW: lamWs,
    cetaU: cetaUs,
    gam: gams,
    KU: KUs,
    Keta: Ketas,
    KW: KWs,
    muW: muWs,
    TU: TUs,
}
expect_zero("Q_tr on the reference branch", sp.simplify(Qtr.subs(subs_equal)))
expect_zero("Q_nt on the reference branch", sp.simplify(Qnt.subs(subs_equal)))
expect_zero("Q_eta on the reference branch", sp.simplify(Qeta.subs(subs_equal)))

subbanner("III. Exact linearization reproduces the Stage-299 monomial-drift vector")
eps = sp.symbols("eps", real=True)
dl_lam, dl_c, dl_gam, dl_KU, dl_Keta, dl_KW, dl_muW, dl_TU = sp.symbols(
    "dl_lambda dl_c dl_gamma dl_KU dl_Keta dl_KW dl_muW dl_TU", real=True
)

lin_subs = {
    lamW: lamWs * sp.exp(eps * dl_lam),
    cetaU: cetaUs * sp.exp(eps * dl_c),
    gam: gams * sp.exp(eps * dl_gam),
    KU: KUs * sp.exp(eps * dl_KU),
    Keta: Ketas * sp.exp(eps * dl_Keta),
    KW: KWs * sp.exp(eps * dl_KW),
    muW: muWs * sp.exp(eps * dl_muW),
    TU: TUs * sp.exp(eps * dl_TU),
}

Qvec_lin = sp.Matrix([
    sp.simplify(sp.diff(Qtr.subs(lin_subs), eps).subs({eps: 0})),
    sp.simplify(sp.diff(Qnt.subs(lin_subs), eps).subs({eps: 0})),
    sp.simplify(sp.diff(Qeta.subs(lin_subs), eps).subs({eps: 0})),
])

Mstar, E_star, F_star = monomial_drift_matrix(chi0s, deltaUs, epsWs, epss)
dx = sp.Matrix([dl_lam, dl_c, dl_gam, dl_KU, dl_Keta, dl_KW, dl_muW, dl_TU])

print("Linearized quotient vector dQ =")
sp.pprint(Qvec_lin)
print("Stage-299 monomial-drift vector M_* dx =")
sp.pprint(Mstar * dx)
expect_zero("dQ - M_* dx", sp.simplify(Qvec_lin - Mstar * dx))

subbanner("IV. Interpretation")
print("The actual moving-throat branch can now be tested in the smallest exact finite language:")
print("  Q_tr  = ln(C_tr/C_tr,ref),")
print("  Q_nt  = ln(C_nt/C_nt,ref),")
print("  Q_eta = ln(epsilon_eta/epsilon_eta,ref).")
print("Their first linearization is exactly the Stage-299 monomial-drift vector.")
print("So the finite quotient coordinates and the infinitesimal defect ledger now fit together exactly.")
