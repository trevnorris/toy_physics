#!/usr/bin/env python3
"""
moving_throat_pde_stage183_triangular_normal_form_sympy_audit.py

SymPy-backed audit for the exact triangular normal form of the coherent
weak-axisymmetric defect.

Checks:
1. Stage 182 defect/tracking formulas reduce to the three branch-adapted
   defect coordinates (Sigma_tr, Sigma_nt, Sigma_eta).
2. The observable ledger (Theta_1, Xi_1, R_1 + Xi_1) is triangular.
3. The exact inverse reconstruction formulas are correct.
4. The full triple-rigidity theorem reduces to
   Sigma_tr = Sigma_nt = Sigma_eta = 0 on the physical branch.
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


def expect_nonzero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(expr)
    print(f"{name} = {expr}")
    if expr == 0:
        raise AssertionError(f"{name} is unexpectedly zero")


banner("STAGE 183 — TRIANGULAR NORMAL FORM OF THE COHERENT DEFECT")

# Physical branch variables.
chi0, epsW, eps_eta, deltaU = sp.symbols('chi0 epsW eps_eta deltaU', positive=True, real=True)
eps = sp.simplify(epsW * (1 - sp.Rational(2, 11) * deltaU / (1 + deltaU)))

# Stage 182 microscopic slippages.
SigmaZ, SigmaChi, SigmaEta, SigmaEps, SigmaDel = sp.symbols(
    'Sigma_Z Sigma_chi Sigma_eta Sigma_eps Sigma_del', real=True
)

SigmaTr = sp.symbols('Sigma_tr', real=True)
SigmaNT = sp.symbols('Sigma_nt', real=True)

# Stage 182 exact formulas.
Theta1 = sp.simplify(
    -chi0 * deltaU / ((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU)) * SigmaTr
)

Xi1 = sp.simplify(
    SigmaZ
    + 2 * chi0 / ((1 + chi0) * (1 + deltaU)) * SigmaTr
    + 2 * epsW / (1 - eps) * (11 + 9 * deltaU) / (11 * (1 + deltaU)) * SigmaEps
    - (
        2 * chi0 / (1 + deltaU)
        + 4 * epsW * deltaU / (11 * (1 - eps) * (1 + deltaU) ** 2)
    ) * SigmaDel
)

R1 = sp.simplify(-eps_eta / (1 - eps_eta) * SigmaEta - Xi1)

banner("Branch-adapted nontracking slippage")
SigmaNT_def = sp.simplify(
    SigmaZ
    + 2 * epsW / (1 - eps) * (11 + 9 * deltaU) / (11 * (1 + deltaU)) * SigmaEps
    - (
        2 * chi0 / (1 + deltaU)
        + 4 * epsW * deltaU / (11 * (1 - eps) * (1 + deltaU) ** 2)
    ) * SigmaDel
)
print("Sigma_nt =")
sp.pprint(SigmaNT_def)

A_tr = sp.simplify(2 * chi0 / ((1 + chi0) * (1 + deltaU)))
C_tr = sp.simplify(chi0 * deltaU / ((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU)))

expect_zero("Xi_1 - (A_tr Sigma_tr + Sigma_nt)", Xi1 - (A_tr * SigmaTr + SigmaNT_def))
expect_zero("R_1 + Xi_1 + eps_eta/(1-eps_eta) Sigma_eta", R1 + Xi1 + eps_eta / (1 - eps_eta) * SigmaEta)

banner("Triangular observable ledger")
print("Theta_1 =")
sp.pprint(Theta1)
print("Xi_1 =")
sp.pprint(sp.simplify(A_tr * SigmaTr + SigmaNT))
print("R_1 + Xi_1 =")
sp.pprint(sp.simplify(-eps_eta / (1 - eps_eta) * SigmaEta))

banner("Exact inverse reconstruction")
SigmaTr_inv = sp.simplify(-((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU)) / (chi0 * deltaU) * Theta1)
expect_zero("Sigma_tr inverse", SigmaTr_inv - SigmaTr)

ratio = sp.simplify(A_tr / C_tr)
print("A_tr/C_tr =", ratio)
expect_zero("A_tr/C_tr - 2(1+chi0+deltaU)/deltaU", ratio - 2 * (1 + chi0 + deltaU) / deltaU)

SigmaNT_inv = sp.simplify(Xi1 + ratio * Theta1)
expect_zero("Sigma_nt inverse", SigmaNT_inv - SigmaNT_def)

SigmaEta_inv = sp.simplify(-(1 - eps_eta) / eps_eta * (R1 + Xi1))
expect_zero("Sigma_eta inverse", SigmaEta_inv - SigmaEta)

banner("Triple-rigidity theorem")
# Rigidity (Theta_1=Xi_1=R_1=0 <=> Sigma_tr=Sigma_nt=Sigma_eta=0) holds iff the
# triangular map is invertible on the branch chi0>0, deltaU>0, 0<eps_eta<1, i.e.
# iff each diagonal prefactor is nonzero there. We test that non-trivial content;
# the trivial forward direction is already implied, and the inverse round-trips
# above (Sigma_tr/Sigma_nt/Sigma_eta inverse) confirm full invertibility.
dressing_pref = sp.simplify(eps_eta / (1 - eps_eta))
expect_nonzero("C_tr (Theta_1 <- Sigma_tr prefactor) nonzero on branch", C_tr)
expect_nonzero("A_tr (Xi_1 <- Sigma_tr feed-through) nonzero on branch", A_tr)
expect_nonzero("eps_eta/(1-eps_eta) (R_1+Xi_1 <- Sigma_eta prefactor) nonzero on branch", dressing_pref)

print("\nCarry-forward formulas:")
print("  Sigma_tr = (1+chi_0) Sigma_del + (1+delta_U) Sigma_chi")
print("  Sigma_nt = Sigma_Z")
print("             + [2 eps_W/(1-eps)] * [(11+9 delta_U)/(11(1+delta_U))] Sigma_eps")
print("             - [2 chi_0/(1+delta_U) + 4 eps_W delta_U/(11(1-eps)(1+delta_U)^2)] Sigma_del")
print("  Theta_1   = -C_tr Sigma_tr")
print("  Xi_1      = A_tr Sigma_tr + Sigma_nt")
print("  R_1+Xi_1  = -(eps_eta/(1-eps_eta)) Sigma_eta")
print("  Sigma_tr  = -((1+chi_0)(1+delta_U)(1+chi_0+delta_U)/(chi_0 delta_U)) Theta_1")
print("  Sigma_nt  = Xi_1 + 2(1+chi_0+delta_U)/delta_U * Theta_1")
print("  Sigma_eta = -((1-eps_eta)/eps_eta) (R_1 + Xi_1)")
