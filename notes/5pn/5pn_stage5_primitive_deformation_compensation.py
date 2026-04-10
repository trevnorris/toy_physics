#!/usr/bin/env python3
"""
5pn_stage5_primitive_deformation_compensation.py

Fifth executable SymPy audit for the 5PN grouped-real-P2 program.

What this script does
---------------------
1. Starts from the explicit Stage-3 isotropic finite-throat overlap prototype.
2. Introduces a primitive weak-axisymmetric microscopic deformation of that
   prototype through the slopes
      dK, dM, d(lambda_B), d(varpi), d(lambda_U), d(lambda_W),
      d(lambda_R), d(Omega_U), d(Omega_W).
3. Computes the induced grouped-lane slope data
      D01, D21, D41, N01
   by exact first-order differentiation at the invariant level.
4. Verifies the grouped obstruction and physical-slope combinations
      K_1 = D21 + D01/9,
      G_1 = N01 - P0 D01,
      Xi_load = N01/N0 - D01/D0.
5. Derives two exact compensation surfaces in primitive variables:
      (a) the even-preserving branch K_1 = 0 solved for dM,
      (b) the odd/normalization branch Xi_load = 0 solved for dK.
6. Records the remaining hidden-even residual
      D41 - (2/3) D21 - D01/27,
   which is the 5PN linear consistency slot still to be controlled.

Interpretation
--------------
This is the first script that turns the explicit Stage-3 overlap model into the
actual primitive deformation problem Stage 4 asked for. It does not yet choose a
unique physical anisotropy mechanism. It says exactly how any chosen primitive
mechanism will feed D01, D21, D41, and N01, and how the even-preserving and
normalization-preserving compensations are solved algebraically.
"""

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


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


# ---------------------------------------------------------------------------
# I. Baseline explicit isotropic prototype (same overlap family as Stage 3)
# ---------------------------------------------------------------------------

banner("I. BASELINE EXPLICIT ISOTROPIC PROTOTYPE")

I_mix = sp.symbols("I_mix", positive=True, real=True)
lamB, lamU, lamW, lamR = sp.symbols("lambda_B lambda_U lambda_W lambda_R", real=True)
varpi = sp.symbols("varpi", positive=True, real=True)
OmegaU, OmegaW = sp.symbols("Omega_U Omega_W", positive=True, real=True)
K, M = sp.symbols("K M", positive=True, real=True)

C = I_mix * lamB
GU = lamU
GW = I_mix * lamW
R = I_mix * lamR

B0 = C**2 / varpi**2
B2 = C**2 / varpi**4
B4 = C**2 / varpi**6

Delta = OmegaU**2 * OmegaW**2 - R**2
S2 = OmegaU**2 + OmegaW**2
Q = GU**2 * OmegaW**2 + 2 * GU * GW * R + GW**2 * OmegaU**2
H = GU**2 + GW**2
Pproto = OmegaU**2 * GW + R * GU

Z0 = Q / Delta
Z2 = (Q * S2 - H * Delta) / Delta**2
Z4 = (Q * (S2**2 - Delta) - S2 * H * Delta) / Delta**3
N0 = Pproto**2 / Delta**2

D0 = K - B0 - Z0
P0 = N0 / D0

print("B0 =", B0)
print("B2 =", B2)
print("B4 =", B4)
print("Z0 =", Z0)
print("Z2 =", Z2)
print("Z4 =", Z4)
print("N0 =", N0)
print("D0 =", D0)
print("P0 =", P0)


# ---------------------------------------------------------------------------
# II. Primitive microscopic deformation slopes
# ---------------------------------------------------------------------------

banner("II. PRIMITIVE MICROSCOPIC DEFORMATION SLOPES")

dK, dM = sp.symbols("dK dM", real=True)
dlamB, dvarpi = sp.symbols("d_lambda_B d_varpi", real=True)
dlamU, dlamW, dlamR = sp.symbols("d_lambda_U d_lambda_W d_lambda_R", real=True)
dOmegaU, dOmegaW = sp.symbols("d_Omega_U d_Omega_W", real=True)

dC = I_mix * dlamB
dGU = dlamU
dGW = I_mix * dlamW
dR = I_mix * dlamR

# BdG moments: exact first-order slopes.
dB0 = 2 * C * dC / varpi**2 - 2 * C**2 * dvarpi / varpi**3
dB2 = 2 * C * dC / varpi**4 - 4 * C**2 * dvarpi / varpi**5
dB4 = 2 * C * dC / varpi**6 - 6 * C**2 * dvarpi / varpi**7

print("dB0 =", dB0)
print("dB2 =", dB2)
print("dB4 =", dB4)

# Maxwell/mixed invariant slopes.
dDelta = 2 * OmegaU * OmegaW**2 * dOmegaU + 2 * OmegaW * OmegaU**2 * dOmegaW - 2 * R * dR
dS2 = 2 * OmegaU * dOmegaU + 2 * OmegaW * dOmegaW
dQ = (
    2 * GU * OmegaW**2 * dGU
    + 2 * GU**2 * OmegaW * dOmegaW
    + 2 * (dGU * GW * R + GU * dGW * R + GU * GW * dR)
    + 2 * GW * OmegaU**2 * dGW
    + 2 * GW**2 * OmegaU * dOmegaU
)
dH = 2 * GU * dGU + 2 * GW * dGW
dPproto = 2 * OmegaU * GW * dOmegaU + OmegaU**2 * dGW + dR * GU + R * dGU

print("\ndDelta  =", dDelta)
print("dS2     =", dS2)
print("dQ      =", dQ)
print("dH      =", dH)
print("dPproto =", dPproto)


# ---------------------------------------------------------------------------
# III. Exact invariant-level derivatives for Z0, Z2, Z4, and N0
# ---------------------------------------------------------------------------

banner("III. INVARIANT-LEVEL SLOPE FORMULAS")

t = sp.symbols("t", real=True)
Delta_t = Delta + t * dDelta
S2_t = S2 + t * dS2
Q_t = Q + t * dQ
H_t = H + t * dH
Pproto_t = Pproto + t * dPproto

Z0_t = Q_t / Delta_t
Z2_t = (Q_t * S2_t - H_t * Delta_t) / Delta_t**2
Z4_t = (Q_t * (S2_t**2 - Delta_t) - S2_t * H_t * Delta_t) / Delta_t**3
N0_t = Pproto_t**2 / Delta_t**2

dZ0 = sp.simplify(sp.diff(Z0_t, t).subs(t, 0))
dZ2 = sp.simplify(sp.diff(Z2_t, t).subs(t, 0))
dZ4 = sp.simplify(sp.diff(Z4_t, t).subs(t, 0))
dN0 = sp.simplify(sp.diff(N0_t, t).subs(t, 0))

print("dZ0 =", dZ0)
print("dZ2 =", dZ2)
print("dZ4 =", dZ4)
print("dN0 =", dN0)

subbanner("Cross-check against the compact dZ0 and dN0 formulas")
expect_zero("dZ0 - (Delta*dQ - Q*dDelta)/Delta^2", dZ0 - (Delta * dQ - Q * dDelta) / Delta**2)
expect_zero("dN0 - 2*P*(Delta*dP - P*dDelta)/Delta^3", dN0 - 2 * Pproto * (Delta * dPproto - Pproto * dDelta) / Delta**3)


# ---------------------------------------------------------------------------
# IV. Grouped operator / transfer slopes induced by the primitive deformation
# ---------------------------------------------------------------------------

banner("IV. D01, D21, D41, N01 FROM THE PRIMITIVE DEFORMATION")

D01 = dK - dB0 - dZ0
D21 = -(dM + dB2 + dZ2)
D41 = -(dB4 + dZ4)
N01 = dN0

print("D01 =", D01)
print("D21 =", D21)
print("D41 =", D41)
print("N01 =", N01)

K1 = sp.simplify(D21 + D01 / 9)
G1 = sp.simplify(N01 - P0 * D01)
Xi_load = sp.simplify(N01 / N0 - D01 / D0)

print("\nK1       =", K1)
print("G1       =", G1)
print("Xi_load  =", Xi_load)
expect_zero("Xi_load - G1/N0", Xi_load - G1 / N0)

subbanner("Hidden-even residual")
hidden_even_residual = sp.simplify(D41 - (sp.Rational(2, 3) * D21 + D01 / 27))
print("D41 - (2/3) D21 - D01/27 =", hidden_even_residual)


# ---------------------------------------------------------------------------
# V. Exact compensation surfaces in primitive variables
# ---------------------------------------------------------------------------

banner("V. EXACT PRIMITIVE COMPENSATION SURFACES")

# 1) Even-preserving compensation: K1 = 0 fixes dM linearly.
dM_even = sp.simplify(sp.solve(sp.Eq(K1, 0), dM)[0])
print("dM on the even-preserving branch =", dM_even)
expect_zero("K1 after imposing dM_even", K1.subs(dM, dM_even))

# 2) Odd/normalization compensation: Xi_load = 0 fixes dK linearly.
dK_odd = sp.simplify(sp.solve(sp.Eq(Xi_load, 0), dK)[0])
print("\ndK on the odd/normalization-preserving branch =", dK_odd)
expect_zero("Xi_load after imposing dK_odd", Xi_load.subs(dK, dK_odd))

subbanner("Useful structural reading")
print("The even-preserving branch fixes only the inertia-side static slope dM.")
print("The odd/normalization-preserving branch fixes only the static loading slope dK.")
print("So in the primitive explicit prototype the two compensations separate cleanly.")


# ---------------------------------------------------------------------------
# VI. Final ledger
# ---------------------------------------------------------------------------

banner("VI. PRIMITIVE-DEFORMATION LEDGER")
print("1. Primitive weak-axisymmetric microscopic deformations feed the explicit")
print("   finite-throat prototype through")
print("      dK, dM, d(lambda_B), d(varpi), d(lambda_U), d(lambda_W),")
print("      d(lambda_R), d(Omega_U), d(Omega_W).")
print("2. Those induce exact grouped-lane slope data")
print("      D01 = dK - dB0 - dZ0,")
print("      D21 = -(dM + dB2 + dZ2),")
print("      D41 = -(dB4 + dZ4),")
print("      N01 = dN0.")
print("3. The two live grouped obstruction channels are")
print("      K1 = D21 + D01/9,")
print("      G1 = N01 - P0 D01,")
print("   and the remaining linear normalization defect is")
print("      Xi_load = N01/N0 - D01/D0 = G1/N0.")
print("4. The even-preserving compensation is algebraic and fixes dM exactly.")
print("5. The odd/normalization-preserving compensation is algebraic and fixes dK exactly.")
print("6. The remaining 5PN linear consistency slot is")
print("      D41 - (2/3) D21 - D01/27,"
      "which is the hidden-even residual that still has to be controlled by the"
      "primitive physical deformation mechanism.")
print("7. So the next honest move is no longer abstract: choose a concrete primitive")
print("   anisotropy mechanism (wall-only, BdG-only, Maxwell/mixed-only, or a mixed")
print("   combination), substitute its slopes here, and test")
print("      K1 = 0,   Xi_load = 0,   hidden-even residual = 0.")
