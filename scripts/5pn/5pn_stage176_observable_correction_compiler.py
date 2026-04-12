#!/usr/bin/env python3
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


def expect_zero(name: str, expr):
    if isinstance(expr, sp.MatrixBase):
        expr = expr.applyfunc(lambda z: sp.simplify(sp.expand(z)))
        print(f"{name} =")
        sp.pprint(expr)
        if any(entry != 0 for entry in expr):
            raise AssertionError(f"{name} is not zero")
    else:
        expr = sp.simplify(sp.expand(expr))
        print(f"{name} = {expr}")
        if expr != 0:
            raise AssertionError(f"{name} is not zero")


"""
5pn_stage176_observable_correction_compiler.py

Stage 176 — exact compiler from the weak-axisymmetric observable triple
(Theta_1, Xi_1, R_1 + Xi_1) to the microscopic dependent-coordinate correction
in (Delta_T, Delta_Keta, Delta_mu).

What this script does
---------------------
1. Rebuilds the Stage-166/167 inverse observable map to the quotient drifts
      (d ln C_tr,* , d ln C_nt,* , d ln epsilon_eta).
2. Composes that inverse with the Stage-175 complementary projector.
3. Produces the exact dependent-coordinate correction needed to represent a given
   observable defect triple as a pure quotient motion.
4. Verifies that applying M_* to this correction reproduces the exact quotient
   drifts implied by the observables.
"""

banner("STAGE 176 — OBSERVABLE-TO-MICROSCOPIC CORRECTION COMPILER")

# ---------------------------------------------------------------------------
# I. Inverse observable map
# ---------------------------------------------------------------------------

subbanner("I. Exact inverse map from observables to quotient drifts")

chi0s, deltaUs = sp.symbols("chi0_star deltaU_star", positive=True, real=True)
eps_eta_star = sp.symbols("epsilon_eta_star", positive=True, real=True)

Theta1, Xi1, Rsum = sp.symbols("Theta_1 Xi_1 Rsum", real=True)

dlnCtr = sp.simplify(- (1 + chi0s) * (1 + deltaUs) * (1 + chi0s + deltaUs) / (chi0s * deltaUs) * Theta1)
dlnCnt = sp.simplify(Xi1 + 2 * (1 + chi0s + deltaUs) / deltaUs * Theta1)
dlnEta = sp.simplify(- (1 - eps_eta_star) / eps_eta_star * Rsum)

print("d ln C_tr,* =", dlnCtr)
print("d ln C_nt,* =", dlnCnt)
print("d ln epsilon_eta =", dlnEta)

# ---------------------------------------------------------------------------
# II. Exact quotient correction in the dependent triple
# ---------------------------------------------------------------------------

subbanner("II. Dependent-coordinate quotient correction")

Estar, Fstar = sp.symbols("E_star F_star", real=True)

# Same pivot inverse as Stage 175.
Pinv = sp.Matrix([
    [1 / (1 + chi0s), 0, 0],
    [0, 0, -1],
    [Fstar / (1 + chi0s), 1, -1],
])
q = sp.Matrix([dlnCtr, dlnCnt, dlnEta])

DeltaT_q, DeltaKeta_q, DeltaMu_q = sp.simplify(Pinv * q)

print("Delta_T^(q) =", DeltaT_q)
print("Delta_Keta^(q) =", DeltaKeta_q)
print("Delta_mu^(q) =", DeltaMu_q)

# More explicit observable forms.
DeltaT_obs = sp.simplify(DeltaT_q)
DeltaKeta_obs = sp.simplify(DeltaKeta_q)
DeltaMu_obs = sp.simplify(DeltaMu_q)

print("\nObservable forms:")
print("Delta_T^(q) =", DeltaT_obs)
print("Delta_Keta^(q) =", DeltaKeta_obs)
print("Delta_mu^(q) =", DeltaMu_obs)

# ---------------------------------------------------------------------------
# III. Compiler check against the exact monomial-drift map
# ---------------------------------------------------------------------------

subbanner("III. Exact compiler check")

M = sp.Matrix([
    [0, 1 + deltaUs, 1 + deltaUs, -(2 + chi0s + deltaUs), 0, 0, 0, 1 + chi0s],
    [2 * (1 + Estar), 0, 2 * Estar, Fstar - Estar, -1, -(2 + Estar), 1, -Fstar],
    [0, 2, 0, -1, -1, 0, 0, 0],
])

Delta_x_obs = sp.Matrix([0, 0, 0, 0, DeltaKeta_obs, 0, DeltaMu_obs, DeltaT_obs])
print("Delta_x^(obs) =")
sp.pprint(Delta_x_obs)

expect_zero("M_* Delta_x^(obs) - (dlnCtr,dlnCnt,dlnEta)", sp.simplify(M * Delta_x_obs - q))

banner("STAGE 176 LEDGER")
print("1. The observable defect triple (Theta_1, Xi_1, R_1 + Xi_1) fixes the quotient drifts")
print("      (d ln C_tr,* , d ln C_nt,* , d ln epsilon_eta) exactly.")
print("2. Those quotient drifts can be represented by a pure dependent-coordinate correction")
print("   supported only in (Delta_T, Delta_Keta, Delta_mu).")
print("3. The exact compiler is")
print("      Delta_T^(q)    = (d ln C_tr,*)/(1+chi_0,*),")
print("      Delta_Keta^(q) = - d ln epsilon_eta,")
print("      Delta_mu^(q)   = d ln C_nt,* - d ln epsilon_eta + F_* d ln C_tr,* /(1+chi_0,*).")
print("4. So once the observable defect is known, the exact microscopic dependent-coordinate")
print("   correction needed to account for it is already fixed.")
