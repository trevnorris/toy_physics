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
5pn_stage177_finite_orbit_restoration_operator.py

Stage 177 — exact finite restoration operator to return any candidate microscopic
branch to a single G_* orbit by adjusting only the dependent triple
(Delta_T, Delta_Keta, Delta_mu).
"""

banner("STAGE 177 — FINITE ORBIT RESTORATION OPERATOR")

# ---------------------------------------------------------------------------
# I. General candidate drift and its quotient mismatch
# ---------------------------------------------------------------------------

subbanner("I. General candidate finite log-ratio vector")

chi0s, deltaUs = sp.symbols("chi0_star deltaU_star", positive=True, real=True)
Estar, Fstar = sp.symbols("E_star F_star", real=True)

Dl, Dc, Dg, DUf, DKe, DWf, Dmu, DTf = sp.symbols(
    "Delta_lambda Delta_c Delta_gamma Delta_U Delta_Keta Delta_KW Delta_mu Delta_T",
    real=True,
)
Delta_x = sp.Matrix([Dl, Dc, Dg, DUf, DKe, DWf, Dmu, DTf])

M = sp.Matrix([
    [0, 1 + deltaUs, 1 + deltaUs, -(2 + chi0s + deltaUs), 0, 0, 0, 1 + chi0s],
    [2 * (1 + Estar), 0, 2 * Estar, Fstar - Estar, -1, -(2 + Estar), 1, -Fstar],
    [0, 2, 0, -1, -1, 0, 0, 0],
])
q = sp.simplify(M * Delta_x)
print("q = M_* Delta_x =")
sp.pprint(q)

# ---------------------------------------------------------------------------
# II. Exact restoration operator
# ---------------------------------------------------------------------------

subbanner("II. Exact restoration by dependent-coordinate shift only")

Pinv = sp.Matrix([
    [1 / (1 + chi0s), 0, 0],
    [0, 0, -1],
    [Fstar / (1 + chi0s), 1, -1],
])

E = sp.zeros(8, 3)
E[7, 0] = 1
E[4, 1] = 1
E[6, 2] = 1

Delta_x_fail = sp.simplify(E * Pinv * q)
Delta_x_restore = sp.simplify(-Delta_x_fail)
Delta_x_orbit = sp.simplify(Delta_x + Delta_x_restore)

print("Delta_x_fail =")
sp.pprint(Delta_x_fail)
print("\nDelta_x_restore =")
sp.pprint(Delta_x_restore)
print("\nDelta_x_orbit = Delta_x + Delta_x_restore =")
sp.pprint(Delta_x_orbit)

expect_zero("M_* Delta_x_restore + q", sp.simplify(M * Delta_x_restore + q))
expect_zero("M_* Delta_x_orbit", sp.simplify(M * Delta_x_orbit))

# Explicit dependent-coordinate restoration laws.
DeltaT_rest = sp.simplify(Delta_x_restore[7])
DeltaKeta_rest = sp.simplify(Delta_x_restore[4])
DeltaMu_rest = sp.simplify(Delta_x_restore[6])

print("\nDependent-coordinate restoration laws:")
print("Delta_T^(rest) =", DeltaT_rest)
print("Delta_Keta^(rest) =", DeltaKeta_rest)
print("Delta_mu^(rest) =", DeltaMu_rest)

# ---------------------------------------------------------------------------
# III. Agreement with the exact orbit law when q = 0
# ---------------------------------------------------------------------------

subbanner("III. Zero-mismatch specialization recovers the Stage-170 orbit law")

# Exact Stage-170 dependent-coordinate orbit law.
DKe_orbit = sp.simplify(2 * Dc - DUf)
DT_orbit = sp.simplify(DUf - (1 + deltaUs) / (1 + chi0s) * (Dg + Dc - DUf))
Dmu_orbit = sp.simplify(
    2 * Dc - DUf + 2 * DWf - 2 * Dl
    - Estar * (2 * Dg + 2 * Dl - DUf - DWf)
    - Fstar * (1 + deltaUs) / (1 + chi0s) * (Dg + Dc - DUf)
)

x_on_orbit = sp.Matrix([Dl, Dc, Dg, DUf, DKe_orbit, DWf, Dmu_orbit, DT_orbit])
expect_zero("M_* x_on_orbit", sp.simplify(M * x_on_orbit))

banner("STAGE 177 LEDGER")
print("1. For any candidate finite log-ratio vector Delta_x, the exact quotient mismatch is q = M_* Delta_x.")
print("2. There is a unique finite restoration supported only in the dependent triple")
print("      (Delta_T, Delta_Keta, Delta_mu)")
print("   that sends the candidate back onto a single G_* orbit.")
print("3. The restored vector Delta_x_orbit satisfies M_* Delta_x_orbit = 0 exactly.")
print("4. When q = 0, the restoration vanishes and one recovers the Stage-170 finite orbit law directly.")
