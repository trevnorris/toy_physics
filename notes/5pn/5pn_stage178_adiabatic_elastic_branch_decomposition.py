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
5pn_stage178_adiabatic_elastic_branch_decomposition.py

Stage 178 — exact decomposition of an adiabatic-elastic candidate branch into
single-orbit motion plus the three exact quotient-failure coordinates.
"""

banner("STAGE 178 — ADIABATIC-ELASTIC BRANCH DECOMPOSITION")

# ---------------------------------------------------------------------------
# I. Adiabatic-elastic closure removes only the scalar off-bundle lane
# ---------------------------------------------------------------------------

subbanner("I. Adiabatic-elastic closure")

print("Assumed closure:")
print("  delta ln Theta_w = 0   (adiabatic wall)")
print("  eps_L = eps_v = eps_T = 0   (elastic / no-entropy / no-fraying)")
print("These conditions remove the scalar off-bundle escape, so the remaining")
print("first-order branch-selection problem is entirely microscopic and 8-dimensional.")

# ---------------------------------------------------------------------------
# II. Candidate microscopic branch and exact quotient split
# ---------------------------------------------------------------------------

subbanner("II. Exact microscopic split of an adiabatic-elastic candidate branch")

chi0s, deltaUs = sp.symbols("chi0_star deltaU_star", positive=True, real=True)
Estar, Fstar = sp.symbols("E_star F_star", real=True)

Dl, Dc, Dg, DUf, DKe, DWf, Dmu, DTf = sp.symbols(
    "Delta_lambda Delta_c Delta_gamma Delta_U Delta_Keta Delta_KW Delta_mu Delta_T",
    real=True,
)
Delta_x_AE = sp.Matrix([Dl, Dc, Dg, DUf, DKe, DWf, Dmu, DTf])

M = sp.Matrix([
    [0, 1 + deltaUs, 1 + deltaUs, -(2 + chi0s + deltaUs), 0, 0, 0, 1 + chi0s],
    [2 * (1 + Estar), 0, 2 * Estar, Fstar - Estar, -1, -(2 + Estar), 1, -Fstar],
    [0, 2, 0, -1, -1, 0, 0, 0],
])
Pinv = sp.Matrix([
    [1 / (1 + chi0s), 0, 0],
    [0, 0, -1],
    [Fstar / (1 + chi0s), 1, -1],
])
E = sp.zeros(8, 3)
E[7, 0] = 1
E[4, 1] = 1
E[6, 2] = 1
Qproj = sp.simplify(E * Pinv * M)
Oproj = sp.simplify(sp.eye(8) - Qproj)

Delta_x_fail = sp.simplify(Qproj * Delta_x_AE)
Delta_x_orbit = sp.simplify(Oproj * Delta_x_AE)
q = sp.simplify(M * Delta_x_AE)

print("Delta_x_fail =")
sp.pprint(Delta_x_fail)
print("\nDelta_x_orbit =")
sp.pprint(Delta_x_orbit)
print("\nq = M_* Delta_x_AE =")
sp.pprint(q)

expect_zero("M_* Delta_x_orbit", sp.simplify(M * Delta_x_orbit))
expect_zero("M_* Delta_x_fail - q", sp.simplify(M * Delta_x_fail - q))

# ---------------------------------------------------------------------------
# III. Observable defect depends only on the quotient piece
# ---------------------------------------------------------------------------

subbanner("III. Observable defect depends only on the quotient piece")

eps_eta_star = sp.symbols("epsilon_eta_star", positive=True, real=True)
q1, q2, q3 = q
Theta1 = sp.simplify(-chi0s * deltaUs / ((1 + chi0s) * (1 + deltaUs) * (1 + chi0s + deltaUs)) * q1)
Xi1 = sp.simplify(2 * chi0s / ((1 + chi0s) * (1 + deltaUs)) * q1 + q2)
Rsum = sp.simplify(-eps_eta_star / (1 - eps_eta_star) * q3)

print("Theta_1 =", Theta1)
print("Xi_1 =", Xi1)
print("R_1 + Xi_1 =", Rsum)

# Replace Delta_x_AE by its orbit part only: the observables must vanish.
q_orbit = sp.simplify(M * Delta_x_orbit)
expect_zero("quotient drift of orbit part", q_orbit)
expect_zero("Theta_1 on orbit part", Theta1.subs({q1: 0, q2: 0, q3: 0}))
expect_zero("Xi_1 on orbit part", Xi1.subs({q1: 0, q2: 0, q3: 0}))
expect_zero("R_1 + Xi_1 on orbit part", Rsum.subs({q1: 0, q2: 0, q3: 0}))

banner("STAGE 178 LEDGER")
print("1. Under the adiabatic-elastic closure, the remaining branch-selection problem is exactly")
print("   the microscopic 8-vector Delta_x_AE modulo the Stage-169 similarity orbit.")
print("2. Any candidate branch splits uniquely as")
print("      Delta_x_AE = Delta_x_orbit + Delta_x_fail,")
print("   with Delta_x_orbit in ker M_* and Delta_x_fail carrying the entire quotient motion.")
print("3. The weak-axisymmetric observables (Theta_1, Xi_1, R_1+Xi_1) depend only on")
print("   Delta_x_fail, not on the orbit part.")
print("4. Therefore the adiabatic-elastic branch is orbit-locked iff Delta_x_fail = 0,")
print("   equivalently iff M_* Delta_x_AE = 0.")
