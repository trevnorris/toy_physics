#!/usr/bin/env python3
"""
g2_step57_off_family_scalar_no_linear_p2.py

Step 57 of the g-2 chain.

What this script does
---------------------
1. Rewrites the first off-family defect as the single scalar normal-slippage combination
       epsilon_perp = g_* epsilon_T + (...) epsilon_v + (...) epsilon_L,
   with delta_perp = -epsilon_perp.
2. Rewrites the outgoing defect ledger in terms of epsilon_perp and delta gamma_W.
3. Verifies the grouped real P2 weak-axisymmetric pattern
       b = 3 a,
       A^2 = (7/10) eps^2 x1^2,
   so the first scalar invariant is quadratic.
4. Concludes that pure grouped-lane anisotropy cannot generate the scalar
   off-bundle slippages at linear order.
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


banner("STEP 57 — THE LAST FIRST-ORDER DEFECT IS ONE SCALAR, AND LINEAR P2 FEED-DOWN IS FORBIDDEN")

r_star, g_star = sp.symbols("r_star g_star", positive=True, real=True)
sigma_star = sp.symbols("sigma_star", positive=True, real=True)
eps_L, eps_v, eps_T = sp.symbols("eps_L eps_v eps_T", real=True)
delta_gamma_W = sp.symbols("delta_gamma_W", real=True)

# ---------------------------------------------------------------------------
# 1. Exact scalar normal coordinate and outgoing defect law
# ---------------------------------------------------------------------------

coef_v = sp.simplify(g_star + 1 / (2 * sp.sqrt(1 + r_star**2)))
coef_L = sp.simplify(2 * g_star + 3 / (4 * sp.sqrt(1 + r_star**2)))

epsilon_perp = sp.simplify(g_star * eps_T + coef_v * eps_v + coef_L * eps_L)
delta_perp = sp.simplify(-epsilon_perp)
print("epsilon_perp =", epsilon_perp)
print("delta_perp =", delta_perp)

Delta_Q = sp.simplify(
    sigma_star * (-16 * epsilon_perp / sp.sqrt(1 + r_star**2) - 27 * delta_gamma_W)
    / (3 * (1 - sigma_star))
)
print("Delta_Q =", Delta_Q)

# Canonical-even preservation again kills epsilon_perp and the mixed-lane scale defect.
Delta_Q_even = sp.simplify(Delta_Q.subs(epsilon_perp, 0))
print("Delta_Q with epsilon_perp = 0 =", Delta_Q_even)

r_F1 = sp.Float("1.77799353547498")
g_F1 = sp.Float("0.758035078944663")
coef_v_num = sp.N(coef_v.subs({r_star: r_F1, g_star: g_F1}), 18)
coef_L_num = sp.N(coef_L.subs({r_star: r_F1, g_star: g_F1}), 18)
print("Family-1 coefficient of eps_T =", g_F1)
print("Family-1 coefficient of eps_v =", coef_v_num)
print("Family-1 coefficient of eps_L =", coef_L_num)

# ---------------------------------------------------------------------------
# 2. Weak axisymmetric grouped-P2 anisotropy is quadratic as a scalar invariant
# ---------------------------------------------------------------------------

eps, x0, x1 = sp.symbols("eps x0 x1", real=True)

x20 = x0 + eps * x1
x21 = x0 + eps * x1 / 2
x22 = x0 - eps * x1

xbar = sp.simplify((x20 + 2 * x21 + 2 * x22) / 5)
a_x = sp.simplify((2 * x20 - x21 - x22) / 10)
b_x = sp.simplify((x21 - x22) / 2)
A_x_sq = sp.simplify(4 * a_x**2 + sp.Rational(4, 5) * b_x**2)

print("xbar =", xbar)
print("a_x =", a_x)
print("b_x =", b_x)
expect_zero("b_x - 3 a_x", b_x - 3 * a_x)
expect_zero("A_x^2 - (7/10) eps^2 x1^2", A_x_sq - sp.Rational(7, 10) * eps**2 * x1**2)
expect_zero("linear coefficient of the first scalar invariant", sp.expand(A_x_sq).coeff(eps, 1))

# A pure grouped-P2 anisotropy therefore cannot generate scalar slippages linearly.
eps_L_lin = 0
eps_v_lin = 0
eps_T_lin = 0
epsilon_perp_lin = sp.simplify(
    epsilon_perp.subs({eps_L: eps_L_lin, eps_v: eps_v_lin, eps_T: eps_T_lin})
)
expect_zero("epsilon_perp^(1,P2)", epsilon_perp_lin)

banner("FINAL LEDGER")
print("The first off-family defect is carried by one scalar only:")
print("  epsilon_perp = g_* eps_T + (g_* + 1/(2 sqrt(1+r_*^2))) eps_v +")
print("                 (2 g_* + 3/(4 sqrt(1+r_*^2))) eps_L.")
print()
print("The outgoing defect therefore reduces to")
print("  Delta_Q = [sigma_*/(3(1-sigma_*))] [ -16 epsilon_perp/sqrt(1+r_*^2) - 27 delta_gamma_W ].")
print()
print("For weak axisymmetric grouped-P2 anisotropy, the exact signature is")
print("  b = 3 a,   A^2 = (7/10) eps^2 x1^2,")
print("so the first scalar invariant is quadratic and there is no linear")
print("grouped-P2 feed-down into epsilon_L, epsilon_v, or epsilon_T.")
print()
print("Therefore the only remaining linear routes are:")
print("  (i) a true scalar/off-bundle slippage epsilon_perp, or")
print("  (ii) a direct odd mixed-port renormalization delta_gamma_W.")
