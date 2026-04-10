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
5pn_stage175_orbit_quotient_projectors.py

Stage 175 — exact orbit/quotient projectors for the Stage-170 finite fibre theorem.

What this script does
---------------------
1. Rebuilds the exact monomial-drift matrix M_* from Stages 169–170.
2. Uses the invertible pivot block in the dependent coordinates
      (Delta_T, Delta_Keta, Delta_mu)
   to construct an exact complementary quotient projector Q and orbit projector O.
3. Proves
      Q^2 = Q,
      O^2 = O,
      QO = OQ = 0,
      M_* O = 0,
      M_* Q = M_*.
4. Shows that any microscopic drift splits uniquely as
      Delta x = Delta x_orbit + Delta x_fail,
   where Delta x_orbit lies in ker M_* and Delta x_fail has support only in the
   dependent triple (T_U, K_eta, mu_W).
"""

banner("STAGE 175 — EXACT ORBIT/QUOTIENT PROJECTORS")

# ---------------------------------------------------------------------------
# I. Exact drift matrix and pivot complement
# ---------------------------------------------------------------------------

subbanner("I. Exact monomial-drift matrix and pivot block")

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

print("M_* =")
sp.pprint(M)

# Use the exact Stage-170 pivot block in the dependent coordinates
# (Delta_T, Delta_Keta, Delta_mu), i.e. columns (7,4,6).
P = M[:, [7, 4, 6]]
print("\nPivot block P_(T,Keta,mu) =")
sp.pprint(P)

expect_zero("det(P) - (1+chi0_*)", sp.factor(P.det()) - (1 + chi0s))
Pinv = sp.simplify(P.inv())
print("P^{-1} =")
sp.pprint(Pinv)

# Embedding of the dependent triple back into the 8-dimensional microscopic space.
E = sp.zeros(8, 3)
E[7, 0] = 1  # Delta_T
E[4, 1] = 1  # Delta_Keta
E[6, 2] = 1  # Delta_mu

# ---------------------------------------------------------------------------
# II. Exact complementary projectors
# ---------------------------------------------------------------------------

subbanner("II. Exact quotient and orbit projectors")

Qproj = sp.simplify(E * Pinv * M)
Oproj = sp.simplify(sp.eye(8) - Qproj)

print("Q_proj =")
sp.pprint(Qproj)
print("\nO_proj =")
sp.pprint(Oproj)

expect_zero("Q^2 - Q", sp.simplify(Qproj * Qproj - Qproj))
expect_zero("O^2 - O", sp.simplify(Oproj * Oproj - Oproj))
expect_zero("Q O", sp.simplify(Qproj * Oproj))
expect_zero("O Q", sp.simplify(Oproj * Qproj))
expect_zero("M_* O", sp.simplify(M * Oproj))
expect_zero("M_* Q - M_*", sp.simplify(M * Qproj - M))

# ---------------------------------------------------------------------------
# III. Unique split of an arbitrary microscopic drift
# ---------------------------------------------------------------------------

subbanner("III. Unique decomposition Delta x = Delta x_orbit + Delta x_fail")

Delta_x_orbit = sp.simplify(Oproj * Delta_x)
Delta_x_fail = sp.simplify(Qproj * Delta_x)

print("Delta_x_orbit =")
sp.pprint(Delta_x_orbit)
print("\nDelta_x_fail =")
sp.pprint(Delta_x_fail)

expect_zero("Delta_x - Delta_x_orbit - Delta_x_fail", Delta_x - Delta_x_orbit - Delta_x_fail)
expect_zero("M_* Delta_x_orbit", sp.simplify(M * Delta_x_orbit))
expect_zero("M_* Delta_x_fail - M_* Delta_x", sp.simplify(M * Delta_x_fail - M * Delta_x))

# The complementary part lives only in the dependent triple.
expect_zero("quotient support in free coordinates", sp.Matrix([
    Delta_x_fail[0], Delta_x_fail[1], Delta_x_fail[2], Delta_x_fail[3], Delta_x_fail[5]
]))

q1, q2, q3 = sp.symbols("q1 q2 q3", real=True)
q = sp.Matrix([q1, q2, q3])
Delta_x_fail_from_q = sp.simplify(E * Pinv * q)
print("\nDelta_x_fail(q) =")
sp.pprint(Delta_x_fail_from_q)
expect_zero("M_* Delta_x_fail(q) - q", sp.simplify(M * Delta_x_fail_from_q - q))

banner("STAGE 175 LEDGER")
print("1. The exact quotient projector is Q = E P^{-1} M, where P is the invertible")
print("   Stage-170 pivot block in (Delta_T, Delta_Keta, Delta_mu).")
print("2. The complementary orbit projector is O = I - Q.")
print("3. Any microscopic drift splits uniquely as")
print("      Delta x = Delta x_orbit + Delta x_fail,")
print("   with Delta x_orbit in ker M_* and Delta x_fail supported only on")
print("      (Delta_T, Delta_Keta, Delta_mu).")
print("4. So the exact finite branch-selection failure coordinates live entirely in the")
print("   dependent triple once the five free similarity directions are quotiented out.")
