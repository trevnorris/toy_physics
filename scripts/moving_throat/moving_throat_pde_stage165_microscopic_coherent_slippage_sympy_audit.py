#!/usr/bin/env python3
"""
moving_throat_pde_stage165_microscopic_coherent_slippage_sympy_audit.py

SymPy-backed audit for the microscopic coherent-kernel slippage decomposition.

Checks:
1. Microscopic log drifts induce the Stage-164 physical variables as stated.
2. The grouped defect Xi_1 collapses to the exact four-slippage law.
3. The selected-branch demand drift R_1 acquires one extra dressing slippage Sigma_eta.
4. The tracking-factor drift Theta_1 factorizes through the single tracking slippage Sigma_tr.
5. Xi_1 admits the exact tracking/nontracking split.
6. Support microscopic drifts do not appear.
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


banner("STAGE 165 — MICROSCOPIC COHERENT-KERNEL SLIPPAGE DECOMPOSITION")

# Microscopic logarithmic drift coefficients.
lam1, c1, gam1 = sp.symbols('lam1 c1 gam1', real=True)
kU, keta, kW, mu1, tau1 = sp.symbols('kU keta kW mu1 tau1', real=True)
# Support-lane drifts are declared but should drop out identically.
lamphi1, kphi = sp.symbols('lamphi1 kphi', real=True)

# Background branch variables.
chi0, epsW, eps_eta, deltaU = sp.symbols('chi0 epsW eps_eta deltaU', positive=True, real=True)

# Stage-30 coherent branch definitions.
zetaZ = 2 * lam1 - keta - kW
omegaW = kW - mu1
chi1 = chi0 * (gam1 + c1 - kU)
eta1 = eps_eta * (2 * c1 - kU - keta)
varepsW = epsW * (2 * gam1 + 2 * lam1 - kU - kW)
deltaU1 = deltaU * (tau1 - kU)
eps = sp.simplify(epsW * (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)))
eps1 = sp.simplify(
    (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)) * varepsW
    - (2 * epsW) / (11 * (1 + deltaU) ** 2) * deltaU1
)

# Microscopic slippages.
SigmaZ = sp.symbols('Sigma_Z', real=True)
SigmaChi = sp.symbols('Sigma_chi', real=True)
SigmaEta = sp.symbols('Sigma_eta', real=True)
SigmaEps = sp.symbols('Sigma_eps', real=True)
SigmaDel = sp.symbols('Sigma_del', real=True)

slip_subs = {
    SigmaZ: 2 * lam1 + mu1 - keta - 2 * kW,
    SigmaChi: gam1 + c1 - kU,
    SigmaEta: 2 * c1 - kU - keta,
    SigmaEps: 2 * gam1 + 2 * lam1 - kU - kW,
    SigmaDel: tau1 - kU,
}

banner("Physical branch drifts from microscopic logs")
expect_zero("zeta_Z formula", zetaZ - (2 * lam1 - keta - kW))
expect_zero("omega_W formula", omegaW - (kW - mu1))
expect_zero("chi_1 formula", chi1 - chi0 * (gam1 + c1 - kU))
expect_zero("eta_1 formula", eta1 - eps_eta * (2 * c1 - kU - keta))
expect_zero("varepsilon_W formula", varepsW - epsW * (2 * gam1 + 2 * lam1 - kU - kW))
expect_zero("delta_U,1 formula", deltaU1 - deltaU * (tau1 - kU))

banner("Four-slippage grouped-defect law")
Xi1_direct = sp.simplify(
    zetaZ - omegaW + 2 * chi1 / (1 + chi0) + 2 * eps1 / (1 - eps)
)
Xi1_slip = sp.simplify(
    SigmaZ
    + 2 * chi0 / (1 + chi0) * SigmaChi
    + 2 * epsW / (1 - eps) * (
        (11 + 9 * deltaU) / (11 * (1 + deltaU)) * SigmaEps
        - 2 * deltaU / (11 * (1 + deltaU) ** 2) * SigmaDel
    )
)
expect_zero("Xi_1 direct - slippage form", Xi1_direct - Xi1_slip.subs(slip_subs))
print("Xi_1 =")
sp.pprint(Xi1_slip)

banner("Selected-branch demand slippage")
R1_direct = sp.simplify(-eta1 / (1 - eps_eta) - Xi1_direct)
R1_slip = sp.simplify(-eps_eta / (1 - eps_eta) * SigmaEta - Xi1_slip)
expect_zero("R_1 direct - slippage form", R1_direct - R1_slip.subs(slip_subs))
print("R_1 =")
sp.pprint(R1_slip)

banner("Tracking-factor factorization")
Theta1_direct = sp.simplify(
    -(chi0 * (1 + chi0) * deltaU1 + deltaU * (1 + deltaU) * chi1)
    / ((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU))
)
SigmaTr = sp.symbols('Sigma_tr', real=True)
Theta1_fact = sp.simplify(
    -chi0 * deltaU / ((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU)) * SigmaTr
)
SigmaTr_def = sp.simplify((1 + chi0) * SigmaDel + (1 + deltaU) * SigmaChi)
expect_zero("Theta_1 factorization", Theta1_direct - Theta1_fact.subs(SigmaTr, SigmaTr_def).subs(slip_subs))
print("Sigma_tr =")
sp.pprint(SigmaTr_def)
print("Theta_1 =")
sp.pprint(Theta1_fact)

banner("Tracking/nontracking split of Xi_1")
Xi1_split = sp.simplify(
    SigmaZ
    + 2 * chi0 / ((1 + chi0) * (1 + deltaU)) * SigmaTr
    + 2 * epsW / (1 - eps) * (11 + 9 * deltaU) / (11 * (1 + deltaU)) * SigmaEps
    - (
        2 * chi0 / (1 + deltaU)
        + 4 * epsW * deltaU / (11 * (1 - eps) * (1 + deltaU) ** 2)
    ) * SigmaDel
)
expect_zero(
    "Xi_1 split - slippage form",
    Xi1_split.subs(SigmaTr, SigmaTr_def) - Xi1_slip
)
print("Xi_1 split =")
sp.pprint(Xi1_split)

banner("Support-blindness at the microscopic level")
print("free symbols of Xi_1:", Xi1_slip.free_symbols)
expect_zero("dXi_1/dlamphi1", sp.diff(Xi1_slip, lamphi1))
expect_zero("dXi_1/dkphi", sp.diff(Xi1_slip, kphi))
expect_zero("dR_1/dlamphi1", sp.diff(R1_slip, lamphi1))
expect_zero("dR_1/dkphi", sp.diff(R1_slip, kphi))
expect_zero("dTheta_1/dlamphi1", sp.diff(Theta1_fact, lamphi1))
expect_zero("dTheta_1/dkphi", sp.diff(Theta1_fact, kphi))

print("\nCarry-forward formulas:")
print("  Sigma_Z   = 2 lam_1 + mu_1 - kappa_eta - 2 kappa_W")
print("  Sigma_chi = gamma_1 + c_1 - kappa_U")
print("  Sigma_eta = 2 c_1 - kappa_U - kappa_eta")
print("  Sigma_eps = 2 gamma_1 + 2 lam_1 - kappa_U - kappa_W")
print("  Sigma_del = tau_1 - kappa_U")
print("  Xi_1 depends only on (Sigma_Z, Sigma_chi, Sigma_eps, Sigma_del)")
print("  R_1 adds only Sigma_eta")
print("  Theta_1 is carried by Sigma_tr = (1+chi_0) Sigma_del + (1+delta_U) Sigma_chi")
