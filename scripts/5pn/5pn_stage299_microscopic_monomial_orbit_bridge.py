#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage296_299_common import *

"""
Stage 299 — microscopic monomial / similarity-orbit bridge.

This stage packages the coherent weak-axisymmetric defect in the exact direct
microscopic monomials C_tr, C_nt, epsilon_eta and proves that the zero-defect
ledger is precisely the tangent-space condition of one exact rank-3 monomial-drift
map, i.e. a five-parameter microscopic similarity orbit.
"""

banner("STAGE 299 — MICROSCOPIC MONOMIAL / SIMILARITY-ORBIT BRIDGE")

# Microscopic variables.
lamW, cetaU, gam, KU, Keta, KW, muW, TU, L, sigma = sp.symbols(
    "lambda_W c_etaU gamma K_U K_etaeff K_Weff mu_W T_U L sigma", positive=True, real=True
)
# Reference-branch values used as frozen exponents.
chi0s, deltaUs, epsWs, epss = sp.symbols("chi0_star deltaU_star epsilonW_star epsilon_star", positive=True, real=True)

# Log-drift coordinates.
dl_lam, dl_c, dl_gam, dl_KU, dl_Keta, dl_KW, dl_muW, dl_TU = sp.symbols(
    "dl_lambda dl_c dl_gamma dl_KU dl_Keta dl_KW dl_muW dl_TU", real=True
)

subbanner("I. Direct microscopic monomials")
chi0 = sp.simplify(gam * cetaU / KU)
deltaU = sp.simplify(sp.pi**2 * TU / (L**2 * KU))
epsEta = sp.simplify(cetaU**2 / (KU * Keta))
epsW = sp.simplify(gam**2 * lamW**2 * sigma / (KU * KW))
ZW_over_OmW2 = sp.simplify(lamW**2 * muW / (Keta * KW**2))

E_star = sp.simplify(2 * epsWs * (11 + 9 * deltaUs) / (11 * (1 - epss) * (1 + deltaUs)))
F_star = sp.simplify(2 * chi0s / (1 + deltaUs) + 4 * epsWs * deltaUs / (11 * (1 - epss) * (1 + deltaUs) ** 2))

Ctr = sp.simplify(chi0 ** (1 + deltaUs) * deltaU ** (1 + chi0s))
Cnt = sp.simplify(ZW_over_OmW2 * epsW ** E_star * deltaU ** (-F_star))

print("C_tr,* =")
sp.pprint(Ctr)
print("C_nt,* =")
sp.pprint(Cnt)
print("epsilon_eta =")
sp.pprint(epsEta)

subbanner("II. Exact direct logarithmic drifts")
dlnCtr = sp.simplify(
    (1 + deltaUs) * (dl_gam + dl_c - dl_KU)
    + (1 + chi0s) * (dl_TU - dl_KU)
)
dlnCnt = sp.simplify(
    (2 * dl_lam + dl_muW - dl_Keta - 2 * dl_KW)
    + E_star * (2 * dl_gam + 2 * dl_lam - dl_KU - dl_KW)
    - F_star * (dl_TU - dl_KU)
)
dlnEta = sp.simplify(2 * dl_c - dl_KU - dl_Keta)

print("delta ln C_tr,* =")
sp.pprint(dlnCtr)
print("delta ln C_nt,* =")
sp.pprint(dlnCnt)
print("delta ln epsilon_eta =")
sp.pprint(dlnEta)

subbanner("III. Exact rank-3 monomial-drift map")
Mstar = sp.Matrix([
    [0, 1 + deltaUs, 1 + deltaUs, -(2 + chi0s + deltaUs), 0, 0, 0, 1 + chi0s],
    [2 * (1 + E_star), 0, 2 * E_star, F_star - E_star, -1, -(2 + E_star), 1, -F_star],
    [0, 2, 0, -1, -1, 0, 0, 0],
])
dx = sp.Matrix([dl_lam, dl_c, dl_gam, dl_KU, dl_Keta, dl_KW, dl_muW, dl_TU])
drift_vec = sp.Matrix([dlnCtr, dlnCnt, dlnEta])

print("M_* =")
sp.pprint(Mstar)
expect_zero("M_* dx - drift_vec", sp.simplify(Mstar * dx - drift_vec))

minor = Mstar[:, [7, 4, 6]].det()
print("det M_*^(dl_TU, dl_Keta, dl_muW) =")
sp.pprint(sp.simplify(minor))
expect_zero("det - (1+chi0_star)", sp.simplify(minor - (1 + chi0s)))
print("rank(M_*) =", Mstar.rank())

subbanner("IV. Exact finite similarity orbit")
Delta_lam, Delta_c, Delta_gam, Delta_KU, Delta_KW = sp.symbols("Delta_lambda Delta_c Delta_gamma Delta_KU Delta_KW", real=True)
Delta_Keta = sp.simplify(2 * Delta_c - Delta_KU)
Delta_TU = sp.simplify(Delta_KU - (1 + deltaUs) * (Delta_gam + Delta_c - Delta_KU) / (1 + chi0s))
Delta_muW = sp.simplify(
    2 * Delta_c - Delta_KU + 2 * Delta_KW - 2 * Delta_lam
    - E_star * (2 * Delta_gam + 2 * Delta_lam - Delta_KU - Delta_KW)
    - F_star * (1 + deltaUs) * (Delta_gam + Delta_c - Delta_KU) / (1 + chi0s)
)

Dvec = sp.Matrix([Delta_lam, Delta_c, Delta_gam, Delta_KU, Delta_Keta, Delta_KW, Delta_muW, Delta_TU])
expect_zero("M_* Delta_x", sp.simplify(Mstar * Dvec))

print("Finite similarity-orbit dependent exponents:")
print("Delta_Keta ="); sp.pprint(Delta_Keta)
print("Delta_TU ="); sp.pprint(Delta_TU)
print("Delta_muW ="); sp.pprint(Delta_muW)

subbanner("V. Final zero-defect theorem")
print("Theta_1 = Xi_1 = R_1 = 0  iff  delta ln C_tr,* = delta ln C_nt,* = delta ln epsilon_eta = 0")
print("Equivalently, the grouped weak-axisymmetric branch tangent lies in ker(M_*),")
print("the exact five-dimensional monomial-preserving similarity orbit.")
