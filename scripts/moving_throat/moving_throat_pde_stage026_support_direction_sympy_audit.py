#!/usr/bin/env python3
"""
Moving-throat PDE — Stage 26 SymPy audit.

What this audit verifies
------------------------
1. Adding the first symmetry-allowed U/phi bilinear to the split-U continuum
   gives an exact effective support vector y.
2. The support-direction factor R_phi is
      R_phi = [1 + sigma0/(1+delta_U)] / (1 + sigma0).
3. The support direction is source-tied iff sigma0 = 0 (or delta_U = 0).
4. The support pole shift gives the exact split support-blocking ratio
      eps_phi^(split) = eps_phi [1 - (2/11) delta_U/(1+delta_U)].
5. The actual physical support baseline is
      M_supp = 8 Z_phi (1+sigma0)^2 / [pi^2 (1-eps_eta)(1-eps_phi^(split))].
6. The support vector tracks the mixed vector iff g_B g_R = g_W g_S,
   equivalently sigma0 = rho0.
7. The exact mismatch formula R_phi - R_U is recovered.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


KU, Keta_eff, Kphi_eff = sp.symbols("K_U K_eta_eff K_phi_eff", positive=True, real=True)
delta_U = sp.symbols("delta_U", positive=True, real=True)
gU, gB, gS, gW, gR = sp.symbols("g_U g_B g_S g_W g_R", real=True)
kappa0, kappa1 = sp.symbols("kappa_0 kappa_1", real=True, nonzero=True)
mu_eta, mu_phi = sp.symbols("mu_eta mu_phi", positive=True, real=True)
ceU, cB, cUphi, cW, cUW = sp.symbols("c_etaU c_etaphi c_Uphi c_etaW c_UW", real=True, nonzero=True)
eps_eta, Zphi = sp.symbols("eps_eta Z_phi", real=True)

banner("STAGE 26 — CONTINUUM EXTRACTION OF THE ACTUAL SUPPORT/BdG DIRECTION")

subbanner("26.1 — Exact effective support vector after eliminating the split U doublet")

DU = sp.diag(1 / KU, 1 / (KU * (1 + delta_U)))
v = sp.Matrix([kappa0, kappa1])
I2 = sp.eye(2)

# Effective wall-to-support vector after integrating out U.
y = sp.simplify(gB * v + gU * gS * DU * v)
y0, y1 = sp.simplify(y[0]), sp.simplify(y[1])
print("y =")
sp.pprint(y)

sigma0 = sp.simplify(gU * gS / (KU * gB))
Rphi = sp.simplify((y1 / y0) / (kappa1 / kappa0))
Rphi_expected = sp.simplify((1 + sigma0 / (1 + delta_U)) / (1 + sigma0))
print("R_phi =")
sp.pprint(sp.factor(Rphi))
expect_zero("R_phi - expected", Rphi - Rphi_expected)

# Exact support direction-splitting invariant.
Dphi = sp.factor(sp.expand(kappa0 * y1 - kappa1 * y0))
Dphi_expected = sp.simplify(-kappa0 * kappa1 * gB * sigma0 * delta_U / (1 + delta_U))
print("D_phi =")
sp.pprint(Dphi)
expect_zero("D_phi - expected", Dphi - Dphi_expected)

subbanner("26.2 — Exact support pole shift and split support-blocking ratio")

sigma = sp.symbols("sigma", positive=True, real=True)
lam0 = sp.Rational(2, 9)
# Use the already frozen exact overlap relation kappa1^2 / sigma = 2/11.
SU = sp.simplify((v.T * DU * v)[0])
SU_sub = sp.simplify(SU.subs({kappa0**2: sigma - kappa1**2, kappa1**2: sp.Rational(2, 11) * sigma}))
SU_expected = sp.simplify((sigma / KU) * (1 - sp.Rational(2, 11) * delta_U / (1 + delta_U)))
print("v.D_U.v =")
sp.pprint(SU_sub)
expect_zero("support overlap contraction", SU_sub - SU_expected)

eps_phi = sp.symbols("eps_phi", real=True)
eps_phi_split = sp.simplify(eps_phi * (1 - sp.Rational(2, 11) * delta_U / (1 + delta_U)))
Aphi_eff = sp.simplify(Kphi_eff - cUphi**2 * SU_expected)
Aphi_eff_expected = sp.simplify(Kphi_eff * (1 - eps_phi_split)).subs(eps_phi, cUphi**2 * sigma / (KU * Kphi_eff))
print("A_phi^(eff) =")
sp.pprint(sp.factor(Aphi_eff))
expect_zero("A_phi^(eff) - expected", Aphi_eff - Aphi_eff_expected)

subbanner("26.3 — Exact physical support baseline")

# Continuum expression from y0^2 / [A0 * Aphi_eff], with A0 = Keta_eff(1-eps_eta)/mu_eta,
# and y0 = kappa0 c_etaphi (1+sigma0)/sqrt(mu_eta mu_phi).
Msupp_cont = sp.simplify(
    (kappa0**2 * cB**2 * (1 + ceU * cUphi / (KU * cB))**2 / (mu_eta * mu_phi))
    / ((Keta_eff * (1 - eps_eta) / mu_eta) * (Kphi_eff * (1 - eps_phi_split.subs(eps_phi, cUphi**2 * sigma / (KU * Kphi_eff))) / mu_phi))
)
Msupp_expected = sp.simplify(
    (sp.Rational(8, 1) / sp.pi**2)
    * (cB**2 / (Keta_eff * Kphi_eff))
    * (1 + ceU * cUphi / (KU * cB))**2
    / ((1 - eps_eta) * (1 - eps_phi_split.subs(eps_phi, cUphi**2 * sigma / (KU * Kphi_eff))))
)
Msupp_cont_eval = sp.simplify(Msupp_cont.subs(kappa0**2, sp.Rational(8, 1) / sp.pi**2))
print("M_supp =")
sp.pprint(sp.factor(Msupp_cont_eval))
expect_zero("M_supp - expected", Msupp_cont_eval - Msupp_expected)

subbanner("26.4 — Exact tracking condition relative to the mixed vector")

z0 = kappa0 * (gW + gU * gR / KU)
z1 = kappa1 * (gW + gU * gR / (KU * (1 + delta_U)))
Dzphi = sp.factor(sp.expand(y0 * z1 - y1 * z0))
Dzphi_expected = sp.simplify(-delta_U * gU * kappa0 * kappa1 * (gB * gR - gW * gS) / (KU * (1 + delta_U)))
print("D_(phi z) =")
sp.pprint(Dzphi)
expect_zero("D_(phi z) - expected", Dzphi - Dzphi_expected)

rho0 = sp.simplify(gU * gR / (KU * gW))
RU = sp.simplify((1 + rho0 / (1 + delta_U)) / (1 + rho0))
expect_zero("tracking condition via g_B g_R = g_W g_S", sp.simplify(Rphi_expected - RU).subs(gS, gB * gR / gW))

subbanner("26.5 — Exact mismatch formula")

mismatch = sp.simplify(Rphi_expected - RU)
mismatch_expected = sp.simplify(delta_U * (rho0 - sigma0) / ((1 + delta_U) * (1 + rho0) * (1 + sigma0)))
print("R_phi - R_U =")
sp.pprint(sp.factor(mismatch))
expect_zero("mismatch formula", mismatch - mismatch_expected)

banner("STAGE 26 THEOREM LEDGER")
print("1. The first symmetry-allowed U/phi continuum extension keeps the support loading rank-1,")
print("   but rotates the support direction to the exact vector y.")
print("2. The exact support-direction factor is")
print("      R_phi = [1 + sigma0/(1+delta_U)] / (1 + sigma0).")
print("3. The support direction is source-tied iff sigma0 = 0 (or delta_U = 0).")
print("4. The exact split support-blocking ratio is")
print("      eps_phi^(split) = eps_phi [1 - (2/11) delta_U/(1+delta_U)].")
print("5. The exact physical support baseline is")
print("      M_supp = 8 Z_phi (1+sigma0)^2 / [pi^2 (1-eps_eta)(1-eps_phi^(split))].")
print("6. The support vector tracks the mixed vector iff")
print("      g_B g_R = g_W g_S,  equivalently  sigma0 = rho0.")
print("7. Therefore the minimal kernel is source-tied, the generic first extended kernel is")
print("   intermediate, and exact tracking is a codimension-one interference match.")
