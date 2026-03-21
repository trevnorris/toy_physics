#!/usr/bin/env python3
"""
2PN P0/P2 mouth-port operator rebuild prototype
-----------------------------------------------
Purpose
-------
Push the constructive 2PN wake story one step past the channel-kernel solve.
The goal here is to show that the solved added 2PN comparable-mass cross block
is not merely a polynomial fit: it closes exactly in a small body-local
mouth-port basis.

Main result
-----------
After carrying forward the frozen 1PN dipole wake, the genuinely new 2PN
quartic tensor sector is *exactly* a positive overlap of:
    P0   : monopole-like scalar port,
    P20  : axisymmetric quadrupole (m=0),
    P21  : quadrupole m=±1 pair,
    P22  : quadrupole m=±2 pair.
So the new conservative 2PN cross sector lives in a minimal
    P0 ⊕ P2
port family.

This is the first direct algebraic bridge from the solved 2PN comparable-mass
cross block to the inner-throat notes' mouth-port language.
"""
from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)


# ---------------------------------------------------------------------------
# Symbols and kinematics
# ---------------------------------------------------------------------------

G, r, cLight = sp.symbols("G r c", positive=True)
mA, mB = sp.symbols("mA mB", positive=True)

vAx, vAy, vAz, vBx, vBy, vBz = sp.symbols(
    "vAx vAy vAz vBx vBy vBz", real=True
)

uAx, uAy = vAx, vAy
uBx, uBy = vBx, vBy

dA = vAz
dB = vBz

vA2 = sp.expand(vAx**2 + vAy**2 + vAz**2)
vB2 = sp.expand(vBx**2 + vBy**2 + vBz**2)
vAB = sp.expand(vAx * vBx + vAy * vBy + vAz * vBz)

uAdotuB = sp.expand(uAx * uBx + uAy * uBy)

# Local potentials felt by each body from the other.
UA = sp.expand(G * mB / r)
UB = sp.expand(G * mA / r)


# ---------------------------------------------------------------------------
# Frozen 1PN wake and solved added 2PN cross target
# ---------------------------------------------------------------------------

L1WakeTarget = sp.expand(-sp.Rational(7, 2) * vAB - sp.Rational(1, 2) * dA * dB)

QuarticTarget = sp.expand(
    -sp.Rational(7, 4) * vAB * (vA2 + vB2)
    - sp.Rational(1, 4) * dA * dB * (vA2 + vB2)
    + sp.Rational(11, 8) * vA2 * vB2
    + sp.Rational(1, 4) * vAB**2
    - sp.Rational(5, 8) * (vA2 * dB**2 + vB2 * dA**2)
    + sp.Rational(3, 2) * vAB * dA * dB
    + sp.Rational(3, 8) * dA**2 * dB**2
)

QuadraticTarget = sp.expand(
    sp.Rational(11, 8) * (mA * vA2 + mB * vB2)
    - sp.Rational(15, 4) * (mA + mB) * vAB
    + sp.Rational(15, 8) * (mA * dA**2 + mB * dB**2)
)

StaticCrossTarget = sp.expand(sp.Rational(5, 4) * UA * UB)

AddedCrossTarget = sp.expand(
    (G * mA * mB / (cLight**4 * r)) * QuarticTarget
    + (G**2 * mA * mB / (cLight**4 * r**2)) * QuadraticTarget
    + 5 * G**3 * mA**2 * mB**2 / (4 * cLight**4 * r**3)
)


# ---------------------------------------------------------------------------
# 1) Frozen 1PN wake as a body-local dipole kernel
# ---------------------------------------------------------------------------

banner("1) Frozen 1PN wake as a body-local dipole kernel")

# Transverse dipole pair (m=±1) and longitudinal dipole (m=0).
Pi1pm_A = sp.Matrix([
    sp.sqrt(sp.Rational(7, 2)) * uAx,
    sp.sqrt(sp.Rational(7, 2)) * uAy,
])
Pi1pm_B = sp.Matrix([
    sp.sqrt(sp.Rational(7, 2)) * uBx,
    sp.sqrt(sp.Rational(7, 2)) * uBy,
])
Pi10_A = 2 * dA
Pi10_B = 2 * dB

L1WakeDipole = sp.expand(-(Pi1pm_A.dot(Pi1pm_B) + Pi10_A * Pi10_B))

print("L1WakeDipole - L1WakeTarget =")
print(sp.expand(L1WakeDipole - L1WakeTarget))


# ---------------------------------------------------------------------------
# 2) Universal leg dressing plus new quartic residual
# ---------------------------------------------------------------------------

banner("2) Universal leg dressing plus new quartic residual")

LegDressed1PN = sp.expand(sp.Rational(1, 2) * (vA2 + vB2) * L1WakeDipole)
QuarticResidual = sp.expand(QuarticTarget - LegDressed1PN)

print("QuarticResidual =")
print(QuarticResidual)


# ---------------------------------------------------------------------------
# 3) P0 / P2 port basis for the new quartic tensor sector
# ---------------------------------------------------------------------------

banner("3) Exact P0/P2 mouth-port factorization of the quartic residual")

# P0 monopole port
Pi0_A = sp.sqrt(5) * vA2 / 2
Pi0_B = sp.sqrt(5) * vB2 / 2

# P2, m=0 axisymmetric quadrupole port
P20_A = sp.expand(3 * dA**2 - vA2)
P20_B = sp.expand(3 * dB**2 - vB2)
Pi20_A = P20_A / 2
Pi20_B = P20_B / 2

# P2, m=±1 real pair
Pi21_A = sp.Matrix([
    sp.sqrt(2) * dA * uAx,
    sp.sqrt(2) * dA * uAy,
])
Pi21_B = sp.Matrix([
    sp.sqrt(2) * dB * uBx,
    sp.sqrt(2) * dB * uBy,
])

# P2, m=±2 real pair
Pi22_A = sp.Matrix([
    (uAx**2 - uAy**2) / (2 * sp.sqrt(2)),
    (2 * uAx * uAy) / (2 * sp.sqrt(2)),
])
Pi22_B = sp.Matrix([
    (uBx**2 - uBy**2) / (2 * sp.sqrt(2)),
    (2 * uBx * uBy) / (2 * sp.sqrt(2)),
])

QuarticPorts = sp.expand(
    Pi0_A * Pi0_B
    + Pi20_A * Pi20_B
    + Pi21_A.dot(Pi21_B)
    + Pi22_A.dot(Pi22_B)
)

print("QuarticPorts - QuarticResidual =")
print(sp.expand(QuarticPorts - QuarticResidual))


# ---------------------------------------------------------------------------
# 4) Scalar TL sector diagonalizes exactly into P0 and P20
# ---------------------------------------------------------------------------

banner("4) The old scalar TL sector diagonalizes exactly into P0 and P20")

TA, LA, TB, LB = sp.symbols("TA LA TB LB", real=True)
KTL = sp.Matrix([
    [sp.Rational(3, 2), sp.Rational(3, 4)],
    [sp.Rational(3, 4), sp.Rational(9, 4)],
])

# y = M x with x=(T,L), y=(P0,P20)=(T+L, -T+2L)
M = sp.Matrix([
    [1, 1],
    [-1, 2],
])
Minv = M.inv()
KPort = sp.simplify(Minv.T * KTL * Minv)

print("KTL in the (P0,P20) basis =")
print(KPort)
print("Expected diagonal result = diag(5/4, 1/4)")


# ---------------------------------------------------------------------------
# 5) Potential dressing only drives the scalar P0 / P20 ports
# ---------------------------------------------------------------------------

banner("5) Exact scalar-potential dressing in the P0/P20 basis")

QuadraticPorts = sp.expand(
    -sp.Rational(15, 4) * (mA + mB) * vAB
    + mA * (2 * vA2 + sp.Rational(5, 8) * P20_A)
    + mB * (2 * vB2 + sp.Rational(5, 8) * P20_B)
)

print("QuadraticPorts - QuadraticTarget =")
print(sp.expand(QuadraticPorts - QuadraticTarget))

QuadraticPortsNormalized = sp.expand(
    -sp.Rational(15, 4) * (UA + UB) * vAB
    + UA * ((4 / sp.sqrt(5)) * Pi0_B + sp.Rational(5, 4) * Pi20_B)
    + UB * ((4 / sp.sqrt(5)) * Pi0_A + sp.Rational(5, 4) * Pi20_A)
)

QuadraticPortsBlock = sp.expand((G * mA * mB / (cLight**4 * r)) * QuadraticPortsNormalized)
QuadraticTargetBlock = sp.expand((G**2 * mA * mB / (cLight**4 * r**2)) * QuadraticTarget)

print("QuadraticPortsBlock - QuadraticTargetBlock =")
print(sp.expand(QuadraticPortsBlock - QuadraticTargetBlock))


# ---------------------------------------------------------------------------
# 6) Static cross term is a pure monopole-potential overlap
# ---------------------------------------------------------------------------

banner("6) Static cross term as a pure monopole-potential overlap")

StaticCrossPorts = sp.expand(sp.Rational(5, 4) * UA * UB)
print("StaticCrossPorts - StaticCrossTarget =")
print(sp.expand(StaticCrossPorts - StaticCrossTarget))


# ---------------------------------------------------------------------------
# 7) Full added 2PN cross block in dipole + P0/P2 language
# ---------------------------------------------------------------------------

banner("7) Full added 2PN cross block in dipole + P0/P2 port language")

AddedCrossPorts = sp.expand(
    (G * mA * mB / (cLight**4 * r))
    * (
        LegDressed1PN
        + QuarticPorts
        - sp.Rational(15, 4) * (UA + UB) * vAB
        + UA * ((4 / sp.sqrt(5)) * Pi0_B + sp.Rational(5, 4) * Pi20_B)
        + UB * ((4 / sp.sqrt(5)) * Pi0_A + sp.Rational(5, 4) * Pi20_A)
        + sp.Rational(5, 4) * UA * UB
    )
)

print("AddedCrossPorts - AddedCrossTarget =")
print(sp.expand(AddedCrossPorts - AddedCrossTarget))


# ---------------------------------------------------------------------------
# 8) Interpretation
# ---------------------------------------------------------------------------

banner("8) Takeaways")
print(
    "1. The frozen 1PN wake is already a body-local dipole kernel: m=±1 transverse\n"
    "   plus m=0 longitudinal.\n"
    "2. After removing the universal leg dressing, the genuinely new quartic 2PN\n"
    "   tensor sector is exactly a positive overlap of P0 plus the full real P2\n"
    "   multiplet {m=0, ±1, ±2}.\n"
    "3. The old scalar (T,L) block diagonalizes *exactly* into the port basis\n"
    "      P0 = v^2,\n"
    "      P20 = 3(v·n)^2 - v^2,\n"
    "   with weights diag(5/4, 1/4).\n"
    "4. The scalar-potential dressing excites only the axisymmetric scalar ports\n"
    "   P0 and P20, with couplings 2 and 5/8 in the unnormalized basis.\n"
    "5. Therefore the full solved added conservative 2PN cross block already has\n"
    "   a minimal constructive mouth-port interpretation:\n"
    "      carried-forward dipole wake  (1PN),\n"
    "      plus a new P0 ⊕ P2 response layer (2PN),\n"
    "      plus scalar-potential source terms on P0/P20,\n"
    "      plus the static monopole-potential overlap (5/4) U_A U_B.\n"
)
