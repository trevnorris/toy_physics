#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp


def banner(title: str) -> None:
    line = '=' * 88
    print('\n' + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = '-' * 88
    print('\n' + line)
    print(title)
    print(line)


def expect_zero(name: str, expr) -> None:
    if isinstance(expr, sp.MatrixBase):
        simp = expr.applyfunc(lambda e: sp.simplify(sp.expand(e)))
        print(f"{name} =")
        sp.pprint(simp)
        if any(e != 0 for e in simp):
            raise AssertionError(f"{name} is not zero")
    else:
        simp = sp.simplify(sp.expand(expr))
        print(f"{name} = {simp}")
        if simp != 0:
            raise AssertionError(f"{name} is not zero")


"""
5pn_stage173_adiabatic_wall_not_sufficient_for_orbit_lock.py

Stage 173 — test whether the adiabatic wall condition by itself forces the
Stage-170 quotient coordinates to be preserved.

What this script does
---------------------
1. Rebuilds the Stage-169/170 monomial-drift matrix M_*.
2. Shows that the quotient drifts
      (d ln C_tr,* , d ln C_nt,* , d ln epsilon_eta)
   are exactly M_* delta x.
3. Exhibits explicit microscopic drift directions for which the quotient drifts
   are nonzero.
4. Concludes that freezing the wall-depth datum alone is not enough to force the
   branch onto a single G_* orbit.

Interpretation
--------------
The adiabatic wall removes one isotropic boundary channel, but the actual orbit-lock
condition still requires preserving the three quotient coordinates themselves.
"""

banner("STAGE 173 — ADIABATIC WALL ALONE IS NOT SUFFICIENT FOR ORBIT LOCK")

chi0s, deltaUs = sp.symbols('chi0_star deltaU_star', positive=True, real=True)
Estar, Fstar = sp.symbols('E_star F_star', real=True)

Dl, Dc, Dg, DUf, DKe, DWf, Dmu, DTf = sp.symbols(
    'Delta_lambda Delta_c Delta_gamma Delta_U Delta_Keta Delta_KW Delta_mu Delta_T',
    real=True,
)
Delta_x = sp.Matrix([Dl, Dc, Dg, DUf, DKe, DWf, Dmu, DTf])

M = sp.Matrix([
    [0, 1 + deltaUs, 1 + deltaUs, -(2 + chi0s + deltaUs), 0, 0, 0, 1 + chi0s],
    [2 * (1 + Estar), 0, 2 * Estar, Fstar - Estar, -1, -(2 + Estar), 1, -Fstar],
    [0, 2, 0, -1, -1, 0, 0, 0],
])

subbanner("1. Exact quotient-drift map")
print("M_* =")
sp.pprint(M)
q = sp.simplify(M * Delta_x)
print("(dln C_tr,* , dln C_nt,* , dln epsilon_eta)^T = M_* Delta_x =")
sp.pprint(q)

subbanner("2. There are explicit microscopic drifts with nonzero quotient motion")

# Example 1: pure Delta_lambda drift
ex1 = sp.Matrix([1, 0, 0, 0, 0, 0, 0, 0])
q1 = sp.simplify(M * ex1)
print("Example 1: Delta_x = e_lambda gives quotient drift")
sp.pprint(q1)
expect_zero("Example 1 first quotient coordinate", q1[0])
# second should be 2(1+E_*), third zero; show not all zero generically
if sp.simplify(q1[1]) == 0:
    raise AssertionError("Expected nonzero nontracking quotient drift for the e_lambda direction.")
expect_zero("Example 1 third quotient coordinate", q1[2])

# Example 2: pure Delta_c drift
ex2 = sp.Matrix([0, 1, 0, 0, 0, 0, 0, 0])
q2 = sp.simplify(M * ex2)
print("Example 2: Delta_x = e_c gives quotient drift")
sp.pprint(q2)
if sp.simplify(q2[0]) == 0 or sp.simplify(q2[2]) == 0:
    raise AssertionError("Expected nonzero tracking and dressing quotient drift for the e_c direction.")

subbanner("3. The exact orbit-lock condition is still M_* Delta_x = 0")

# Solve exact finite fibre equations for the dependent variables.
sol = sp.solve(list(q), [DTf, DKe, Dmu], dict=True)
if len(sol) != 1:
    raise AssertionError("Expected a unique solve for (Delta_T, Delta_Keta, Delta_mu).")
sol = sol[0]
print("Finite orbit-lock solution for (Delta_T, Delta_Keta, Delta_mu):")
for k, v in sol.items():
    print(f"{k} = {sp.simplify(v)}")

banner("FINAL STAGE-173 LEDGER")
print("1. The exact quotient drifts are always M_* Delta_x.")
print("2. The adiabatic wall condition freezes the isotropic wall-depth datum, but M_* has no Theta_w coordinate.")
print("3. Therefore freezing Theta_w alone does not force")
print("      d ln C_tr,* = d ln C_nt,* = d ln epsilon_eta = 0.")
print("4. Explicit counterexamples exist: for example, pure Delta_lambda or pure Delta_c motion produces nonzero quotient drift.")
print("5. So the wall-freeze condition is necessary for removing one isotropic boundary channel, but it is not sufficient for single-orbit locking.")
print("6. The actual orbit-lock test remains exactly")
print("      M_* Delta_x = 0,  equivalently  preservation of (C_tr,* , C_nt,* , epsilon_eta).")
