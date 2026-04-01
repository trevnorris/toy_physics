#!/usr/bin/env python3
"""
moving_throat_pde_stage168_microscopic_monomials_sympy_audit.py

SymPy-backed audit for Stage 168.

Checks:
1. Direct microscopic ratio drifts for chi_0, delta_U, epsilon_W, Z_W/Omega_W^2,
   and epsilon_eta.
2. Direct microscopic tracking monomial drift = Sigma_tr.
3. Direct microscopic nontracking monomial drift = Sigma_nt.
4. Dressing monomial drift = Sigma_eta and the observable complement law.
5. Exact zero-defect compatibility solve for tau1, kappa_eta, mu1.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


def first_ratio_drift(ratio: sp.Expr, eps: sp.Symbol, lam: sp.Symbol) -> sp.Expr:
    """Return the first grouped drift of a multiplicative ratio with ratio(0)=1."""
    return sp.simplify(sp.diff(ratio, eps).subs(eps, 0) / lam)


banner("STAGE 168 — DIRECT MICROSCOPIC MONOMIALS")

chi0s, deltaUs, epsWs, epss = sp.symbols("chi0s deltaUs epsWs epss", positive=True, real=True)
lam1, c1, gam1, kU, keta, kW, mu1, tau1 = sp.symbols(
    "lam1 c1 gam1 kU keta kW mu1 tau1", real=True
)
eps, lam = sp.symbols("eps lam", real=True)

gamma_ref, cetaU_ref, TU_ref = sp.symbols("gamma c_etaU T_U", positive=True, real=True)
KU_ref, Keta_ref, KWeff_ref = sp.symbols("K_U K_eta_eff K_W_eff", positive=True, real=True)
lamW_ref, muW_ref = sp.symbols("lambda_W mu_W", positive=True, real=True)

# Microscopic logarithmic slippages from Stage 165.
Sigma_chi = gam1 + c1 - kU
Sigma_delta = tau1 - kU
Sigma_Z = 2 * lam1 + mu1 - keta - 2 * kW
Sigma_eps = 2 * gam1 + 2 * lam1 - kU - kW
Sigma_eta = 2 * c1 - kU - keta

# Stage-166 branch-adapted coordinates at the frozen reference branch.
Sigma_tr = (1 + chi0s) * Sigma_delta + (1 + deltaUs) * Sigma_chi
E_star = 2 * epsWs / (1 - epss) * (11 + 9 * deltaUs) / (11 * (1 + deltaUs))
F_star = 2 * chi0s / (1 + deltaUs) + 4 * epsWs * deltaUs / (11 * (1 - epss) * (1 + deltaUs) ** 2)
Sigma_nt = Sigma_Z + E_star * Sigma_eps - F_star * Sigma_delta

print("Sigma_tr =", sp.simplify(Sigma_tr))
print("Sigma_nt =", sp.simplify(Sigma_nt))
print("Sigma_eta =", sp.simplify(Sigma_eta))

# Independent multiplicative perturbations of the microscopic variables.
gamma_p = gamma_ref * sp.exp(eps * lam * gam1)
cetaU_p = cetaU_ref * sp.exp(eps * lam * c1)
TU_p = TU_ref * sp.exp(eps * lam * tau1)
KU_p = KU_ref * sp.exp(eps * lam * kU)
Keta_p = Keta_ref * sp.exp(eps * lam * keta)
KWeff_p = KWeff_ref * sp.exp(eps * lam * kW)
lamW_p = lamW_ref * sp.exp(eps * lam * lam1)
muW_p = muW_ref * sp.exp(eps * lam * mu1)

chi0_0 = gamma_ref * cetaU_ref / KU_ref
chi0_p = gamma_p * cetaU_p / KU_p
chi_ratio = sp.cancel(chi0_p / chi0_0)

deltaU_0 = TU_ref / KU_ref
deltaU_p = TU_p / KU_p
deltaU_ratio = sp.cancel(deltaU_p / deltaU_0)

epsW_0 = gamma_ref**2 * lamW_ref**2 / (KU_ref * KWeff_ref)
epsW_p = gamma_p**2 * lamW_p**2 / (KU_p * KWeff_p)
epsW_ratio = sp.cancel(epsW_p / epsW_0)

Zratio_0 = lamW_ref**2 * muW_ref / (Keta_ref * KWeff_ref**2)
Zratio_p = lamW_p**2 * muW_p / (Keta_p * KWeff_p**2)
Zratio = sp.cancel(Zratio_p / Zratio_0)

epseta0 = sp.simplify(cetaU_ref**2 / (KU_ref * Keta_ref))
epseta_p = cetaU_p**2 / (KU_p * Keta_p)
epseta_ratio = sp.cancel(epseta_p / epseta0)

Sigma_chi_direct = first_ratio_drift(chi_ratio, eps, lam)
Sigma_delta_direct = first_ratio_drift(deltaU_ratio, eps, lam)
Sigma_eps_direct = first_ratio_drift(epsW_ratio, eps, lam)
Sigma_Z_direct = first_ratio_drift(Zratio, eps, lam)
Sigma_eta_direct = first_ratio_drift(epseta_ratio, eps, lam)

banner("Microscopic ratio drifts")
expect_zero("d ln chi0 - Sigma_chi", Sigma_chi_direct - Sigma_chi)
expect_zero("d ln delta_U - Sigma_delta", Sigma_delta_direct - Sigma_delta)
expect_zero("d ln epsilon_W - Sigma_eps", Sigma_eps_direct - Sigma_eps)
expect_zero("d ln (Z_W/Omega_W^2) - Sigma_Z", Sigma_Z_direct - Sigma_Z)
expect_zero("d ln epsilon_eta - Sigma_eta", Sigma_eta_direct - Sigma_eta)

banner("Tracking monomial")
# C_tr,* = chi0^(1+deltaU*) deltaU^(1+chi0*)
Ctr_ratio = sp.simplify(chi_ratio ** (1 + deltaUs) * deltaU_ratio ** (1 + chi0s))
Sigma_tr_direct = first_ratio_drift(Ctr_ratio, eps, lam)
expect_zero("d ln C_tr,* - Sigma_tr", Sigma_tr_direct - Sigma_tr)

banner("Nontracking monomial")
# C_nt,* = (Z_W/Omega_W^2) * eps_W^(E_*) * delta_U^(-F_*)
Cnt_ratio = sp.simplify(Zratio * epsW_ratio ** E_star * deltaU_ratio ** (-F_star))
Sigma_nt_direct = first_ratio_drift(Cnt_ratio, eps, lam)
expect_zero("d ln C_nt,* - Sigma_nt", Sigma_nt_direct - Sigma_nt)

banner("Dressing monomial")
expect_zero("d ln epsilon_eta - Sigma_eta", Sigma_eta_direct - Sigma_eta)

banner("Observable triangular law in microscopic monomials")
C_tr_star = chi0s * deltaUs / ((1 + chi0s) * (1 + deltaUs) * (1 + chi0s + deltaUs))
A_tr_star = 2 * chi0s / ((1 + chi0s) * (1 + deltaUs))
Theta1 = sp.simplify(-C_tr_star * Sigma_tr)
Xi1 = sp.simplify(A_tr_star * Sigma_tr + Sigma_nt)
Rcombo = sp.simplify(-epseta0 / (1 - epseta0) * Sigma_eta)
expect_zero("Theta_1 monomial law", Theta1 - (-C_tr_star * Sigma_tr_direct))
expect_zero("Xi_1 monomial law", Xi1 - (A_tr_star * Sigma_tr_direct + Sigma_nt_direct))
Rcombo_ratio = sp.simplify((1 - epseta0 * epseta_ratio) / (1 - epseta0))
Rcombo_direct = first_ratio_drift(Rcombo_ratio, eps, lam)
expect_zero("R_1 + Xi_1 complement law", Rcombo_direct - Rcombo)
print("Theta1 =", sp.simplify(Theta1))
print("Xi1 =", sp.simplify(Xi1))
print("R1 + Xi1 =", sp.simplify(Rcombo))

banner("Exact zero-defect compatibility solve")
tau_sol = sp.solve(sp.Eq(Sigma_tr, 0), tau1)[0]
keta_sol = sp.solve(sp.Eq(Sigma_eta, 0), keta)[0]
mu_sol = sp.simplify(keta + 2 * kW - 2 * lam1 - E_star * Sigma_eps + F_star * Sigma_delta)

print("tau1 =", sp.simplify(tau_sol))
print("kappa_eta =", sp.simplify(keta_sol))
print("mu1 =", sp.simplify(mu_sol))
mu_sol_full = sp.simplify(mu_sol.subs({tau1: tau_sol, keta: keta_sol}))
print("mu1 on full zero-defect branch =", mu_sol_full)

expect_zero("tracking substitution", Sigma_tr.subs(tau1, tau_sol))
expect_zero("dressing substitution", Sigma_eta.subs(keta, keta_sol))
expect_zero(
    "nontracking substitution",
    Sigma_nt.subs(
        {
            tau1: tau_sol,
            keta: keta_sol,
            mu1: sp.simplify(mu_sol.subs({tau1: tau_sol, keta: keta_sol})),
        }
    ),
)

print("\nCarry-forward formulas:")
print("  C_tr,*  = chi0^(1+deltaU*) deltaU^(1+chi0*)")
print("  C_nt,*  = (Z_W/Omega_W^2) eps_W^(E_*) delta_U^(-F_*)")
print("  epsilon_eta = c_{etaU}^2 / (K_U K_eta^(eff))")
print("  Zero defect iff these three microscopic monomials are invariant at first grouped order.")
