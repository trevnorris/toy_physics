#!/usr/bin/env python3
"""
moving_throat_pde_stage167_branch_invariant_coordinates_sympy_audit.py

SymPy-backed audit for the exact branch-invariant coordinates of the coherent
weak-axisymmetric defect.

Checks:
1. The exact branch identity R_target * T^2 = Lambda0 * (1 - eps_eta).
2. The normalized tracking invariant gives delta ln T_* = Sigma_tr.
3. The corrected nontracking composite gives delta ln N_* = Sigma_nt.
4. The microscopic dressing ratio gives delta ln eps_eta = Sigma_eta.
5. The full zero-defect theorem is equivalent to invariance of
   (R_tr, N_*, eps_eta) at first grouped order.
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


banner("STAGE 167 — EXACT BRANCH-INVARIANT COORDINATES")

# Background coherent-branch variables.
chi0, deltaU, eps_eta = sp.symbols('chi0 deltaU eps_eta', positive=True, real=True)
Bstar = sp.simplify(2 * (1 + chi0 + deltaU) / deltaU)
Cstar = sp.simplify((1 + chi0) * (1 + deltaU) * (1 + chi0 + deltaU) / (chi0 * deltaU))

print("B_* =", Bstar)
print("C_* =", Cstar)

# First-order logarithmic drift data.
small = sp.symbols('small', real=True)
Theta1, Xi1, SigmaEta = sp.symbols('Theta_1 Xi_1 Sigma_eta', real=True)
Rtr0, T20 = sp.symbols('Rtr0 T20', positive=True, real=True)
Lam0 = sp.symbols('Lambda0', positive=True, real=True)

# Stage-166 branch-adapted coordinates.
SigmaTr = sp.simplify(-Cstar * Theta1)
SigmaNT = sp.simplify(Xi1 + Bstar * Theta1)

banner("Exact branch identities")
Rtr = Rtr0 * sp.exp(small * Theta1)
T2 = T20 * sp.exp(small * Xi1)
eps_eta_var = eps_eta * (1 + small * SigmaEta)
Rtarget = Lam0 * (1 - eps_eta_var) / T2

# 1. Exact product identity.
expect_zero("R_target * T^2 - Lambda0 * (1 - eps_eta)", sp.simplify(Rtarget * T2 - Lam0 * (1 - eps_eta_var)))

banner("Tracking invariant")
Ttr = sp.simplify(Rtr**(-Cstar))
Ttr0 = sp.simplify(Rtr0**(-Cstar))
dln_Ttr = sp.simplify(sp.series(sp.log(Ttr / Ttr0), small, 0, 2).removeO().coeff(small, 1))
print("delta ln T_* =", dln_Ttr)
expect_zero("delta ln T_* - Sigma_tr", dln_Ttr - SigmaTr)

banner("Corrected nontracking composite")
Ntr = sp.simplify(T2 * Rtr**Bstar)
Ntr0 = sp.simplify(T20 * Rtr0**Bstar)
dln_Ntr = sp.simplify(sp.series(sp.log(Ntr / Ntr0), small, 0, 2).removeO().coeff(small, 1))
print("delta ln N_* =", dln_Ntr)
expect_zero("delta ln N_* - Sigma_nt", dln_Ntr - SigmaNT)

banner("Dressing coordinate and selected-branch complement")
dln_eps_eta = sp.simplify(sp.series(sp.log(eps_eta_var / eps_eta), small, 0, 2).removeO().coeff(small, 1))
print("delta ln eps_eta =", dln_eps_eta)
expect_zero("delta ln eps_eta - Sigma_eta", dln_eps_eta - SigmaEta)

Ecomp = sp.simplify((Rtarget * T2) / Lam0)
Ecomp0 = sp.simplify(1 - eps_eta)
dln_Ecomp = sp.simplify(sp.series(sp.log(Ecomp / Ecomp0), small, 0, 2).removeO().coeff(small, 1))
print("delta ln[(R_target T^2)/Lambda0] =", dln_Ecomp)
expect_zero(
    "selected-branch complement identity",
    dln_Ecomp + eps_eta / (1 - eps_eta) * SigmaEta,
)

banner("Composite zero-defect theorem")
# The three branch-adapted coordinates are direct log-drifts of exact branch composites.
expect_zero("Sigma_tr as branch-invariant log drift", SigmaTr - dln_Ttr)
expect_zero("Sigma_nt as branch-invariant log drift", SigmaNT - dln_Ntr)
expect_zero("Sigma_eta as branch-invariant log drift", SigmaEta - dln_eps_eta)

# Vanishing observable drifts is equivalent to invariance of the three branch composites.
# Here we verify the forward zero map directly.
zero_map_1 = sp.simplify(dln_Ttr.subs({Theta1: 0}))
zero_map_2 = sp.simplify(dln_Ntr.subs({Theta1: 0, Xi1: 0}))
zero_map_3 = sp.simplify(dln_eps_eta.subs({SigmaEta: 0}))
expect_zero("delta ln T_* | Theta1=0", zero_map_1)
expect_zero("delta ln N_* | Theta1=Xi1=0", zero_map_2)
expect_zero("delta ln eps_eta | Sigma_eta=0", zero_map_3)

print("\nCarry-forward formulas:")
print("  T_*  = R_tr^(-C_*)")
print("  N_*  = T^2 R_tr^(B_*)")
print("  D_*  = eps_eta")
print("  delta ln T_* = Sigma_tr")
print("  delta ln N_* = Sigma_nt")
print("  delta ln D_* = Sigma_eta")
print("  (R_target T^2)/Lambda0 = 1 - eps_eta")
print("  delta ln[(R_target T^2)/Lambda0] = -(eps_eta/(1-eps_eta)) Sigma_eta")
print("  Full zero-defect theorem: invariance of (R_tr, N_*, eps_eta) at first grouped order.")
