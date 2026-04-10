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
5pn_stage179_dependent_coordinate_mismatch.py

Stage 179 — exact mismatch formulas for the dependent microscopic coordinates.

What this script does
---------------------
1. Starts from the exact Stage-170 single-orbit dependent-coordinate laws.
2. Compares an arbitrary candidate branch to those laws.
3. Shows that the mismatch of (Delta_T, Delta_Keta, Delta_mu) is in one-to-one
   correspondence with the quotient drift triple and therefore with
   (Theta_1, Xi_1, R_1+Xi_1).
"""

banner("STAGE 179 — DEPENDENT-COORDINATE MISMATCH FORMULAS")

chi0s, deltaUs = sp.symbols("chi0_star deltaU_star", positive=True, real=True)
Estar, Fstar = sp.symbols("E_star F_star", real=True)
eps_eta_star = sp.symbols("epsilon_eta_star", positive=True, real=True)

Dl, Dc, Dg, DUf, DKe, DWf, Dmu, DTf = sp.symbols(
    "Delta_lambda Delta_c Delta_gamma Delta_U Delta_Keta Delta_KW Delta_mu Delta_T",
    real=True,
)

# Exact single-orbit laws from Stage 170.
DKe_orbit = sp.simplify(2 * Dc - DUf)
DT_orbit = sp.simplify(DUf - (1 + deltaUs) / (1 + chi0s) * (Dg + Dc - DUf))
Dmu_orbit = sp.simplify(
    2 * Dc - DUf + 2 * DWf - 2 * Dl
    - Estar * (2 * Dg + 2 * Dl - DUf - DWf)
    - Fstar * (1 + deltaUs) / (1 + chi0s) * (Dg + Dc - DUf)
)

mT = sp.simplify(DTf - DT_orbit)
mKeta = sp.simplify(DKe - DKe_orbit)
mMu = sp.simplify(Dmu - Dmu_orbit)

print("mismatch_T =", mT)
print("mismatch_Keta =", mKeta)
print("mismatch_mu =", mMu)

# Quotient drifts from the exact matrix.
q1 = sp.simplify((1 + deltaUs) * (Dc + Dg - DUf) + (1 + chi0s) * (DTf - DUf + Dg + Dc - DUf))
# easier to use exact formulas directly
q1 = sp.simplify((1 + deltaUs) * (Dc + Dg - DUf) + (1 + chi0s) * mT)
q1 = sp.simplify((1 + chi0s) * mT)
q3 = sp.simplify(-mKeta)
q2 = sp.simplify(mMu - Fstar * mT - mKeta)

print("q1 =", q1)
print("q2 =", q2)
print("q3 =", q3)

# Rebuild the exact monomial-drift map to verify the identification.
M = sp.Matrix([
    [0, 1 + deltaUs, 1 + deltaUs, -(2 + chi0s + deltaUs), 0, 0, 0, 1 + chi0s],
    [2 * (1 + Estar), 0, 2 * Estar, Fstar - Estar, -1, -(2 + Estar), 1, -Fstar],
    [0, 2, 0, -1, -1, 0, 0, 0],
])
q_direct = sp.simplify(M * sp.Matrix([Dl, Dc, Dg, DUf, DKe, DWf, Dmu, DTf]))
expect_zero("q1 direct - (1+chi0_*) mismatch_T", q_direct[0] - q1)
expect_zero("q2 direct - (mismatch_mu - F_* mismatch_T - mismatch_Keta)", q_direct[1] - q2)
expect_zero("q3 direct + mismatch_Keta", q_direct[2] - q3)

Theta1 = sp.simplify(-chi0s * deltaUs / ((1 + chi0s) * (1 + deltaUs) * (1 + chi0s + deltaUs)) * q1)
Xi1 = sp.simplify(2 * chi0s / ((1 + chi0s) * (1 + deltaUs)) * q1 + q2)
Rsum = sp.simplify(-eps_eta_star / (1 - eps_eta_star) * q3)

print("\nObservable form:")
print("Theta_1 =", Theta1)
print("Xi_1 =", Xi1)
print("R_1 + Xi_1 =", Rsum)

banner("STAGE 179 LEDGER")
print("1. The exact single-orbit law fixes the dependent triple")
print("      (Delta_T, Delta_Keta, Delta_mu)")
print("   as functions of the five free microscopic drifts.")
print("2. Any failure of the actual branch to stay on one G_* orbit is therefore exactly")
print("   the mismatch of that dependent triple from the Stage-170 orbit law.")
print("3. The quotient drifts are simply")
print("      d ln C_tr,* = (1+chi_0,*) mismatch_T,")
print("      d ln epsilon_eta = - mismatch_Keta,")
print("      d ln C_nt,* = mismatch_mu - F_* mismatch_T - mismatch_Keta.")
print("4. So the remaining dynamical theorem gap has been localized completely:")
print("   the PDE must show that the dependent microscopic coordinates follow the exact orbit law.")
