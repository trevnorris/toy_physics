#!/usr/bin/env python3
"""Independent S11c-c spec review: B0c composition algebra (Grok decision leg).

Prints the linear system that closes (delta_p, J, v_bulk) against (V, mu_theta)
and the coefficient of the bulk amplitude. No conclusions.
"""
from sympy import Eq, Matrix, collect, simplify, solve, symbols

print("=" * 72)
print("S11CC_GROK_LEG: PERMEABLE CLOSURE LINEAR SYSTEM")
print("=" * 72)

# Symbols. Z is the (possibly operator-valued) bulk DtN; treat as a scalar
# symbol here to see the algebraic structure of B0c, which S11c-c §3b copies.
rho_m, Z = symbols("rho_m Z")
LamA, LamV, LamX = symbols("Lambda_A Lambda_V Lambda_X")
V, mu_s, mu_theta, rho_br = symbols("V mu_s mu_theta rho_br_bg")
delta_p, J, v_bulk = symbols("delta_p J v_bulk")

# Three supplied relations (S11b B0c / S11c-c §1d / §3b):
#   (B)  delta_p = Z * v_bulk
#   (K)  v_bulk  = V + J / rho_m
#   (C)  J       = LamA * A + LamV * V
#   (A)  A       = mu_s - delta_p / rho_m
#   (T)  t       = -(delta_p + LamX * A) n     (not needed for the (delta_p, J) solve)
A = mu_s - delta_p / rho_m
eqB = Eq(delta_p, Z * v_bulk)
eqK = Eq(v_bulk, V + J / rho_m)
eqC = Eq(J, LamA * A + LamV * V)

print("\n--- supplied system ---")
print("B:", eqB)
print("K:", eqK)
print("C:", eqC)
print("A:", A)

# Eliminate v_bulk using K into B: delta_p = Z (V + J/rho_m)
eqB2 = Eq(delta_p, Z * (V + J / rho_m))
print("\nB after K:", eqB2)

sols = solve([eqB2, eqC], [delta_p, J], dict=True)
print("\n--- solve for (delta_p, J) in (V, mu_s) ---")
print("n_solutions =", len(sols))
for sol in sols:
    dp = simplify(sol[delta_p])
    Jj = simplify(sol[J])
    print("delta_p =", dp)
    print("J       =", Jj)
    print("delta_p collected in V, mu_s:")
    print("  ", collect(dp.expand(), [V, mu_s]))
    print("coeff_V(delta_p)    =", simplify(dp.expand().coeff(V)))
    print("coeff_mu_s(delta_p) =", simplify(dp.expand().coeff(mu_s)))
    print("coeff_V(J)          =", simplify(Jj.expand().coeff(V)))
    print("coeff_mu_s(J)       =", simplify(Jj.expand().coeff(mu_s)))

# Locus: face equation loses dependence on the bulk amplitude.
# Bulk amplitude here is the free half-space mode coefficient. Equivalently,
# after using B, the map from v_bulk to the face residuals.
# Drive the face with prescribed (V, mu_s) and ask when the homogeneous
# (V=mu_s=0) problem admits a free (delta_p, J, v_bulk) != 0, or when the
# driven map drops rank.
print("\n--- homogeneous (V=0, mu_s=0) kernel ---")
eqB2h = Eq(delta_p, Z * (J / rho_m))
eqCh = Eq(J, LamA * (-delta_p / rho_m))
solh = solve([eqB2h, eqCh], [delta_p, J], dict=True)
print("homogeneous solutions =", solh)
# Matrix form on (delta_p, J):
#   delta_p - (Z/rho_m) J           = Z V
#   (LamA/rho_m) delta_p + J        = LamA mu_s + LamV V
M = Matrix([
    [1, -Z / rho_m],
    [LamA / rho_m, 1],
])
rhs = Matrix([Z * V, LamA * mu_s + LamV * V])
print("system matrix M =")
print(M)
print("det M =", simplify(M.det()))
print("rhs =", rhs)

print("\n--- degenerate locus equations (det M = 0) ---")
print("det_eq:", Eq(simplify(M.det()), 0))
print("det expanded:", simplify(M.det()))

print("\n--- T-i vs traction: which law carries Lambda_X ---")
# T-i (S11c-a §4): shape-differentiate J - LamA A - LamV V = 0. No Lambda_X.
# T-d traction: t = -(delta_p + LamX A) n. Lambda_X lives here.
print("T_i_closure_contains_Lambda_X =", False)
print("T_d_traction_contains_Lambda_X =", True)
print("theta_row_S11c_b = evolution_mass_balance - sum closure_shape_deriv")
print("Lambda_X_enters_slab_via = traction / virtual_work (thickness and in-plane rows)")

print("\nDONE")



