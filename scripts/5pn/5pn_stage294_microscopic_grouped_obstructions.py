#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage289_292_common import *

"""
Stage 294 — exact microscopic decomposition of the linear grouped outlet obstructions.
"""

banner("STAGE 294 — MICROSCOPIC GROUPED OBSTRUCTIONS")

P0 = sp.symbols("P_0", real=True)
dK, dM = sp.symbols("deltaK deltaM", real=True)
dB0, dB2, dB4 = sp.symbols("deltaB0 deltaB2 deltaB4", real=True)
dZ0, dZ2, dZ4 = sp.symbols("deltaZ0 deltaZ2 deltaZ4", real=True)
dN0 = sp.symbols("deltaN0", real=True)

subbanner("I. Full-bundle decomposition of the two outlet obstructions")
dD0 = sp.simplify(dK - dB0 - dZ0)
dD2 = sp.simplify(-(dM + dB2 + dZ2))
dD4 = sp.simplify(-(dB4 + dZ4))

Kobs = sp.simplify(dD2 + dD0 / 9)
Gobs = sp.simplify(dN0 - P0 * dD0)

W_A = sp.simplify(dK / 9 - dM)
B_A = sp.simplify(dB2 + dB0 / 9)
Z_A = sp.simplify(dZ2 + dZ0 / 9)
N_A = sp.simplify(dN0 + P0 * dZ0)

print("K_A obstruction =")
sp.pprint(Kobs)
print("G_A obstruction =")
sp.pprint(Gobs)
expect_zero("K_A = W_A - B_A - Z_A", sp.simplify(Kobs - (W_A - B_A - Z_A)))
expect_zero("G_A = -P0 dK + P0 dB0 + N_A", sp.simplify(Gobs - (-P0 * dK + P0 * dB0 + N_A)))

subbanner("II. Exact BdG contribution")
c, wv, dc, dw = sp.symbols("c varpi deltac deltavarpi", positive=True, real=True)
B0_mode = c**2 / wv**2
B2_mode = c**2 / wv**4

dB0_mode = sp.simplify(sp.diff(B0_mode, c) * dc + sp.diff(B0_mode, wv) * dw)
dB2_mode = sp.simplify(sp.diff(B2_mode, c) * dc + sp.diff(B2_mode, wv) * dw)
Bcomb_mode = sp.simplify(dB2_mode + dB0_mode / 9)
print("delta B0 (one mode) =")
sp.pprint(dB0_mode)
print("delta B2 (one mode) =")
sp.pprint(dB2_mode)
print("B_A contribution (one mode) =")
sp.pprint(Bcomb_mode)
expect_zero(
    "BdG combination formula",
    sp.simplify(
        Bcomb_mode
        - (
            2 * c * (wv**-4 + sp.Rational(1, 9) * wv**-2) * dc
            - 2 * c**2 * (2 * wv**-5 + sp.Rational(1, 9) * wv**-3) * dw
        )
    ),
)

subbanner("III. Exact Maxwell/mixed conservative and outgoing bundles")
Delta, S, Q, G, P = sp.symbols("Delta S Q G P", positive=True, real=True)
dDelta, dS, dQ, dG, dP = sp.symbols("deltaDelta deltaS deltaQ deltaG deltaP", real=True)

Z0 = Q / Delta
Z2 = (Q * S - G * Delta) / Delta**2
N0_expr = P**2 / Delta**2

dZ0_expr = sp.simplify(sp.diff(Z0, Q) * dQ + sp.diff(Z0, Delta) * dDelta)
dZ2_expr = sp.simplify(
    sp.diff(Z2, Q) * dQ + sp.diff(Z2, S) * dS + sp.diff(Z2, G) * dG + sp.diff(Z2, Delta) * dDelta
)
dN0_expr = sp.simplify(sp.diff(N0_expr, P) * dP + sp.diff(N0_expr, Delta) * dDelta)

Zcomb = sp.simplify(dZ2_expr + dZ0_expr / 9)
Ncomb = sp.simplify(dN0_expr + P0 * dZ0_expr)

print("delta Z0 =")
sp.pprint(dZ0_expr)
print("delta Z2 =")
sp.pprint(dZ2_expr)
print("Z_A contribution =")
sp.pprint(Zcomb)
print("N_A contribution =")
sp.pprint(Ncomb)

expect_zero(
    "dZ0 formula",
    sp.simplify(dZ0_expr - (Delta * dQ - Q * dDelta) / Delta**2),
)
expect_zero(
    "dN0 formula",
    sp.simplify(dN0_expr - (2 * P * dP / Delta**2 - 2 * P**2 * dDelta / Delta**3)),
)

subbanner("IV. Primitive port-variation dictionary")
OmU2, OmW2, R, gU, gW = sp.symbols("OmegaU2 OmegaW2 R gU gW", positive=True, real=True)
dOmU2, dOmW2, dR, dgU, dgW = sp.symbols("deltaOmegaU2 deltaOmegaW2 deltaR deltagU deltagW", real=True)

Delta_p = OmU2 * OmW2 - R**2
S_p = OmU2 + OmW2
G_p = gU**2 + gW**2
P_p = OmU2 * gW + R * gU
Q_p = gU**2 * OmW2 + 2 * gU * gW * R + gW**2 * OmU2

dDelta_p = sp.simplify(sp.diff(Delta_p, OmU2) * dOmU2 + sp.diff(Delta_p, OmW2) * dOmW2 + sp.diff(Delta_p, R) * dR)
dS_p = sp.simplify(sp.diff(S_p, OmU2) * dOmU2 + sp.diff(S_p, OmW2) * dOmW2)
dG_p = sp.simplify(sp.diff(G_p, gU) * dgU + sp.diff(G_p, gW) * dgW)
dP_p = sp.simplify(
    sp.diff(P_p, OmU2) * dOmU2 + sp.diff(P_p, R) * dR + sp.diff(P_p, gU) * dgU + sp.diff(P_p, gW) * dgW
)
dQ_p = sp.simplify(
    sp.diff(Q_p, OmU2) * dOmU2 + sp.diff(Q_p, OmW2) * dOmW2 + sp.diff(Q_p, R) * dR + sp.diff(Q_p, gU) * dgU + sp.diff(Q_p, gW) * dgW
)

print("delta Delta =")
sp.pprint(dDelta_p)
print("delta S =")
sp.pprint(dS_p)
print("delta G =")
sp.pprint(dG_p)
print("delta P =")
sp.pprint(dP_p)
print("delta Q =")
sp.pprint(dQ_p)

subbanner("V. Interpretation")
print("The full linear grouped outlet problem depends only on four microscopic defect")
print("bundles: a wall baseline piece, a BdG support piece, a conservative")
print("Maxwell/mixed piece, and an outgoing-transfer piece.")
print("On the weak axisymmetric branch this collapses one stage further to two scalar")
print("amplitudes, K_1 and G_1.")
