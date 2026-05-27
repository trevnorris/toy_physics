#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

def banner(t: str):
    print("\n" + "="*88)
    print(t)
    print("="*88)

z, L = sp.symbols('z L', positive=True, real=True)
k = sp.pi/(2*L)

banner("STAGE 126 — EXPLICIT POSITIVE SOURCE FAMILIES")

# Self-matched derivative profile
sigma_match = k*sp.cos(k*z)
norm_match = sp.simplify(sp.integrate(sigma_match, (z, 0, L)))
g_match = sp.simplify(sp.integrate(sigma_match*sp.cos(k*z), (z, 0, L)))

print("sigma_match normalization =", norm_match)
print("g_match =", g_match)
if sp.simplify(norm_match - 1) != 0 or sp.simplify(g_match - sp.pi/4) != 0:
    raise AssertionError("Self-matched profile formulas failed.")

R = sp.sqrt(4107 - 100*sp.pi**2)
gminus = sp.simplify((2*R - 37*sp.sqrt(3))/(20*sp.pi))
delta_g = sp.simplify(sp.pi/4 - gminus)
traction_ratio = sp.simplify((sp.pi/4)/gminus)

print("\ng_-^F1 =", gminus)
print("delta_g_match =", delta_g)
print("T_- / T_match =", traction_ratio)
print("g_match numeric =", sp.N(sp.pi/4, 20))
print("g_- numeric     =", sp.N(gminus, 20))
print("traction ratio  =", sp.N(traction_ratio, 20))

# Convex positive family
xi = sp.symbols('xi', real=True)
sigma_xi = (1-xi)*k*sp.cos(k*z) + xi/L
norm_xi = sp.simplify(sp.integrate(sigma_xi, (z, 0, L)))
g_xi = sp.simplify(sp.integrate(sigma_xi*sp.cos(k*z), (z, 0, L)))

print("\nsigma_xi normalization =", norm_xi)
print("g_xi =", g_xi)
if sp.simplify(norm_xi - 1) != 0:
    raise AssertionError("Convex family normalization failed.")

# Positivity: sigma_xi(z) >= 0 on z in [0, L], xi in [0, 1].
# Decomposition: sigma_xi = (1 - xi)*k*cos(k*z) + xi/L. Each term is >= 0 on
# the stated domain (cos(k*z) >= 0 since k*z in [0, pi/2]; xi/L > 0 for xi >= 0).
min_match = sp.calculus.util.minimum(k*sp.cos(k*z), z, sp.Interval(0, L))
print("min sigma_match on [0,L] =", min_match)
if sp.simplify(min_match) != 0:
    raise AssertionError("sigma_match minimum on [0,L] should be 0 (at z=L).")
min_xi0 = sp.calculus.util.minimum(sigma_xi.subs(xi, 0), z, sp.Interval(0, L))
val_xi1 = sp.simplify(sigma_xi.subs(xi, 1))
print("min sigma_xi(xi=0) on [0,L] =", min_xi0)
print("sigma_xi(xi=1) =", val_xi1)
if sp.simplify(min_xi0) != 0:
    raise AssertionError("sigma_xi(xi=0) minimum on [0,L] should be 0.")
if sp.simplify(val_xi1 - 1/L) != 0:
    raise AssertionError("sigma_xi(xi=1) should equal 1/L.")
sigma_corner = sp.simplify(sigma_xi.subs([(z, L), (xi, 0)]))
print("sigma_xi(z=L, xi=0) =", sigma_corner)
if sp.simplify(sigma_corner) != 0:
    raise AssertionError("sigma_xi(z=L, xi=0) should equal 0.")

xi_star = sp.simplify(sp.solve(sp.Eq(g_xi, gminus), xi)[0])
print("xi_* =", xi_star)
print("xi_* numeric =", sp.N(xi_star, 20))

g_xi_star = sp.simplify(g_xi.subs(xi, xi_star) - gminus)
print("g_xi(xi_*) - g_- =", g_xi_star)
if sp.simplify(g_xi_star) != 0:
    raise AssertionError("Exact convex-family compensation failed.")

# Interval check
print("\n2/pi numeric =", sp.N(2/sp.pi, 20))
print("pi/4 numeric =", sp.N(sp.pi/4, 20))
interval_check = bool(sp.N(2/sp.pi) < sp.N(gminus) < sp.N(sp.pi/4))
print("Check 2/pi < g_- < pi/4 ->", interval_check)
if not interval_check:
    raise AssertionError("g_-^F1 does not lie strictly between 2/pi and pi/4.")
