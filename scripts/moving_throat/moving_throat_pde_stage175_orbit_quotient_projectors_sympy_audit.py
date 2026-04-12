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


def expect_zero(name: str, expr) -> None:
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


banner("STAGE 175 — EXACT ORBIT/QUOTIENT PROJECTORS AND THE MICROSCOPIC ORBIT-LOCK SPLIT")

chi0s, deltaUs, Estar, Fstar = sp.symbols(
    "chi0_star deltaU_star E_star F_star", positive=True, real=True
)

# Ordered drift basis:
# (Delta_lambda, Delta_c, Delta_gamma, Delta_U, Delta_Keta, Delta_W, Delta_mu, Delta_T)
Mstar = sp.Matrix(
    [
        [0, 1 + deltaUs, 1 + deltaUs, -(2 + chi0s + deltaUs), 0, 0, 0, 1 + chi0s],
        [2 * (1 + Estar), 0, 2 * Estar, Fstar - Estar, -1, -(2 + Estar), 1, -Fstar],
        [0, 2, 0, -1, -1, 0, 0, 0],
    ]
)

subbanner("I. Exact monomial-drift map and dependent pivot block")
print("M_* =")
sp.pprint(Mstar)

# Dependent coordinates in the order (T, K_eta, mu)
T_idx = 7
Keta_idx = 4
mu_idx = 6
Pdep = sp.Matrix.hstack(Mstar[:, T_idx], Mstar[:, Keta_idx], Mstar[:, mu_idx])
Pdep_expected = sp.Matrix(
    [
        [1 + chi0s, 0, 0],
        [-Fstar, -1, 1],
        [0, -1, 0],
    ]
)
print("P_(T,K_eta,mu) =")
sp.pprint(Pdep)
expect_zero("pivot block - expected block", Pdep - Pdep_expected)
print("det(P_(T,K_eta,mu)) =")
sp.pprint(sp.factor(Pdep.det()))
expect_zero("det(P) - (1+chi0_*)", sp.factor(Pdep.det()) - (1 + chi0s))

Pdep_inv = sp.simplify(Pdep.inv())
Pdep_inv_expected = sp.Matrix(
    [
        [1 / (1 + chi0s), 0, 0],
        [0, 0, -1],
        [Fstar / (1 + chi0s), 1, -1],
    ]
)
print("P_(T,K_eta,mu)^(-1) =")
sp.pprint(Pdep_inv)
expect_zero("pivot inverse - expected inverse", Pdep_inv - Pdep_inv_expected)
expect_zero("P_inv P - I", Pdep_inv * Pdep - sp.eye(3))
expect_zero("P P_inv - I", Pdep * Pdep_inv - sp.eye(3))

subbanner("II. Exact canonical quotient section")
Edep = sp.zeros(8, 3)
Edep[T_idx, 0] = 1
Edep[Keta_idx, 1] = 1
Edep[mu_idx, 2] = 1

Sdep = sp.simplify(Edep * Pdep_inv)
Sdep_expected = sp.Matrix(
    [
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, -1],
        [0, 0, 0],
        [Fstar / (1 + chi0s), 1, -1],
        [1 / (1 + chi0s), 0, 0],
    ]
)
print("S_(T,K_eta,mu) = E P^(-1) =")
sp.pprint(Sdep)
expect_zero("section - expected section", Sdep - Sdep_expected)
expect_zero("M_* S - I_3", sp.simplify(Mstar * Sdep) - sp.eye(3))

subbanner("III. Exact complementary projectors")
Qquot = sp.simplify(Sdep * Mstar)
Oorb = sp.simplify(sp.eye(8) - Qquot)
print("Q_quot =")
sp.pprint(Qquot)
print("O_orb =")
sp.pprint(Oorb)
expect_zero("Q^2 - Q", sp.simplify(Qquot * Qquot - Qquot))
expect_zero("O^2 - O", sp.simplify(Oorb * Oorb - Oorb))
expect_zero("Q O", sp.simplify(Qquot * Oorb))
expect_zero("O Q", sp.simplify(Oorb * Qquot))
expect_zero("M_* O", sp.simplify(Mstar * Oorb))
expect_zero("M_* Q - M_*", sp.simplify(Mstar * Qquot - Mstar))

subbanner("IV. Quotient-failure support lies only on the dependent triple")
# rows 0,1,2,3,5 must vanish in Q
for idx in [0, 1, 2, 3, 5]:
    expect_zero(f"Q row {idx}", sp.Matrix([Qquot[idx, :]]))

Dl, Dc, Dg, DU, DKeta, DW, Dmu, DT = sp.symbols(
    "Delta_lambda Delta_c Delta_gamma Delta_U Delta_Keta Delta_W Delta_mu Delta_T",
    real=True,
)
Dx = sp.Matrix([Dl, Dc, Dg, DU, DKeta, DW, Dmu, DT])
qvec = sp.simplify(Mstar * Dx)
qtr, qnt, qeta = qvec
print("q = M_* Delta x =")
sp.pprint(qvec)

Dx_fail = sp.simplify(Qquot * Dx)
Dx_fail_expected = sp.Matrix(
    [
        0,
        0,
        0,
        0,
        -qeta,
        0,
        Fstar / (1 + chi0s) * qtr + qnt - qeta,
        qtr / (1 + chi0s),
    ]
)
print("Delta x_fail = Q Delta x =")
sp.pprint(Dx_fail)
expect_zero("Q Delta x - expected dependent-triple support", Dx_fail - Dx_fail_expected)
expect_zero("M_* Delta x_fail - q", sp.simplify(Mstar * Dx_fail - qvec))

subbanner("V. Exact orbit piece preserves free coordinates and enforces the single-orbit law")
alpha = sp.simplify((1 + deltaUs) / (1 + chi0s))
Dx_orbit = sp.simplify(Oorb * Dx)
Dx_orbit_expected = sp.Matrix(
    [
        Dl,
        Dc,
        Dg,
        DU,
        2 * Dc - DU,
        DW,
        2 * Dc - DU + 2 * DW - 2 * Dl - Estar * (2 * Dg + 2 * Dl - DU - DW)
        - Fstar * alpha * (Dg + Dc - DU),
        DU - alpha * (Dg + Dc - DU),
    ]
)
print("Delta x_orbit = O Delta x =")
sp.pprint(Dx_orbit)
expect_zero("O Delta x - expected orbit law", sp.simplify(Dx_orbit - Dx_orbit_expected))
expect_zero("M_* Delta x_orbit", sp.simplify(Mstar * Dx_orbit))
expect_zero("Delta x - (orbit + fail)", sp.simplify(Dx - Dx_orbit - Dx_fail))

subbanner("VI. Exact orbit-lock equivalence")
# Canonical quotient representative carried by an independent q packet
qtr_i, qnt_i, qeta_i = sp.symbols("q_tr q_nt q_eta", real=True)
qind = sp.Matrix([qtr_i, qnt_i, qeta_i])
Dx_q = sp.simplify(Sdep * qind)
print("Canonical quotient representative S q =")
sp.pprint(Dx_q)
expect_zero("M_* (S q) - q", sp.simplify(Mstar * Dx_q - qind))
expect_zero("Q (S q) - S q", sp.simplify(Qquot * Dx_q - Dx_q))
expect_zero("O (S q)", sp.simplify(Oorb * Dx_q))

zero_q = sp.Matrix([0, 0, 0])
expect_zero("S * 0", sp.simplify(Sdep * zero_q))
expect_zero("orbit-lock packet on orbit-law branch", sp.simplify(Mstar * Dx_orbit_expected))

banner("STAGE 175 LEDGER")
print("1. The finite Packet-B quotient coordinates are exactly q = M_* Delta x.")
print("2. The dependent triple (Delta_T, Delta_Keta, Delta_mu) forms an exact pivot block")
print("   with determinant 1 + chi0_* > 0.")
print("3. The exact quotient section S = E P^(-1) is a right inverse of M_* and is supported")
print("   only on the dependent triple.")
print("4. Therefore the complementary projectors")
print("      Q_quot = S M_*,   O_orb = I - Q_quot")
print("   split any microscopic drift uniquely into orbit motion plus quotient failure.")
print("5. The quotient-failure piece lives only on (Delta_T, Delta_Keta, Delta_mu), while")
print("   the orbit piece keeps the five free similarity coordinates unchanged and enforces")
print("   the exact Stage-170 single-orbit law on the dependent triple.")
print("6. So orbit lock is exactly the zero-set condition")
print("      Q_quot Delta x = 0  <=>  M_* Delta x = 0.")
