#!/usr/bin/env python3
"""
Moving-Throat PDE — Stage 49 SymPy audit

Purpose
-------
Compress the explicit thin-wall parent thresholds into a single dimensionless wall
figure of merit.

Main checks
-----------
1. Define
      W_wall = 4*pi*a^2*L^2*J1*V0^2 / (T_X*ell).
2. Then the Stage-48 thin-wall thresholds become
      W_fail = Pe_req/Delta_inf,
      W_suff = Pe_req/Delta_0.
3. Therefore:
      W_wall <= W_fail  -> fail,
      W_wall >= W_suff  -> succeed,
   and only the intermediate band needs the full fixed-point solve.
4. In the constant-compressibility wall layer,
      W_H = 4*pi*a^2*L^2*I_f*V0^2 / (H_w*T_X*ell),
   with the same thresholds.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


banner("STAGE 66 — DIMENSIONLESS WALL FIGURE OF MERIT")

V0, ell, a, L, J1, TX = sp.symbols("V0 ell a L J1 T_X", positive=True, real=True)
Pe_req, Delta0, Deltainf = sp.symbols("Pe_req Delta_0 Delta_inf", positive=True, real=True)
Hw, If = sp.symbols("H_w I_f", positive=True, real=True)

W_wall = sp.simplify(4 * sp.pi * a**2 * L**2 * J1 * V0**2 / (TX * ell))
W_fail = sp.simplify(Pe_req / Deltainf)
W_suff = sp.simplify(Pe_req / Delta0)

print("W_wall =", W_wall)
print("W_fail =", W_fail)
print("W_suff =", W_suff)

# Stage-48 thresholds
V0_fail_sq = sp.simplify(TX * ell * Pe_req / (4 * sp.pi * a**2 * L**2 * J1 * Deltainf))
V0_suff_sq = sp.simplify(TX * ell * Pe_req / (4 * sp.pi * a**2 * L**2 * J1 * Delta0))

expect_zero("W_wall(V0_fail)^2 - W_fail", sp.simplify(W_wall.subs(V0**2, V0_fail_sq) - W_fail))
expect_zero("W_wall(V0_suff)^2 - W_suff", sp.simplify(W_wall.subs(V0**2, V0_suff_sq) - W_suff))

banner("MONOTONICITY")

# W_wall is manifestly monotone in V0^2, a^2, L^2, J1 and inverse ell, T_X
# (notes section 3: six signed directions).
Vp = sp.symbols("Vp", positive=True, real=True)
W_Vp = sp.simplify(W_wall.subs(V0**2, Vp))

dW_dV0sq = sp.simplify(sp.diff(W_Vp, Vp))
dW_da    = sp.simplify(sp.diff(W_wall, a))
dW_dL    = sp.simplify(sp.diff(W_wall, L))
dW_dell  = sp.simplify(sp.diff(W_wall, ell))
dW_dJ1   = sp.simplify(sp.diff(W_wall, J1))
dW_dTX   = sp.simplify(sp.diff(W_wall, TX))

print("dW/d(V0^2) =", dW_dV0sq)
print("dW/da =", dW_da)
print("dW/dL =", dW_dL)
print("dW/dell =", dW_dell)
print("dW/dJ1 =", dW_dJ1)
print("dW/dT_X =", dW_dTX)

assert sp.simplify(dW_dV0sq > 0) is sp.true, "dW/d(V0^2) should be positive"
assert sp.simplify(dW_da    > 0) is sp.true, "dW/da should be positive"
assert sp.simplify(dW_dL    > 0) is sp.true, "dW/dL should be positive"
assert sp.simplify(dW_dell  < 0) is sp.true, "dW/dell should be negative"
assert sp.simplify(dW_dJ1   > 0) is sp.true, "dW/dJ1 should be positive"
assert sp.simplify(dW_dTX   < 0) is sp.true, "dW/dT_X should be negative"

banner("CONSTANT-COMPRESSIBILITY FORM")
W_H = sp.simplify(4 * sp.pi * a**2 * L**2 * If * V0**2 / (Hw * TX * ell))
print("W_H =", W_H)
expect_zero("J1 = I_f/H_w reduction", sp.simplify(W_wall.subs(J1, If / Hw) - W_H))
expect_zero("W_H(V0_fail)^2 - W_fail", sp.simplify(W_H.subs(V0**2, Hw * TX * ell * Pe_req / (4 * sp.pi * a**2 * L**2 * If * Deltainf)) - W_fail))
expect_zero("W_H(V0_suff)^2 - W_suff", sp.simplify(W_H.subs(V0**2, Hw * TX * ell * Pe_req / (4 * sp.pi * a**2 * L**2 * If * Delta0)) - W_suff))

banner("STAGE 66 AUDIT PASSED")
print("1. The explicit thin-wall parent branch is controlled by one dimensionless wall figure of merit W_wall.")
print("2. Fail/succeed are equivalent to comparing W_wall against the exact bounds Pe_req/Delta_inf and Pe_req/Delta_0.")
print("3. In the constant-compressibility wall layer, the same statement holds for W_H.")
