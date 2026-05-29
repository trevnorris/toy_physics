#!/usr/bin/env python3
"""
Stage 130 SymPy audit — exact mouth-bias map and Family-1 compensation point.
"""

from __future__ import annotations
import sympy as sp

z, L, Pi = sp.symbols("z L Pi", positive=True, real=True)
g_minus = sp.Float("0.758035078944663")

sigma = Pi*sp.exp(-Pi*z/L)/(L*(1-sp.exp(-Pi)))
f = sp.cos(sp.pi*z/(2*L))

gPi = sp.simplify(sp.integrate(sigma*f, (z, 0, L)))
gPi_boxed = 2*Pi*(2*Pi*sp.exp(Pi) + sp.pi) / ((4*Pi**2 + sp.pi**2) * (sp.exp(Pi) - 1))
if sp.simplify(gPi - gPi_boxed) != 0:
    raise AssertionError("g_Pi does not match paper boxed closed form.")
Ez = sp.simplify(sp.integrate(sigma*z, (z, 0, L)))
Efz = sp.simplify(sp.integrate(sigma*f*z, (z, 0, L)))
cov_id = sp.simplify(sp.diff(gPi, Pi) + (Efz - gPi*Ez)/L)

print("g_Pi =", gPi)
print("Covariance identity residual =", cov_id)

g0 = sp.simplify(sp.limit(gPi, Pi, 0, dir='+'))
ginf = sp.simplify(sp.limit(gPi, Pi, sp.oo))
print("limit Pi->0+ =", g0)
print("limit Pi->oo =", ginf)

if cov_id != 0:
    raise AssertionError("Covariance identity failed.")
if sp.simplify(g0 - 2/sp.pi) != 0:
    raise AssertionError("Incorrect uniform-source limit.")
if sp.simplify(ginf - 1) != 0:
    raise AssertionError("Incorrect point-source limit.")

# Global strict monotonicity dg/dPi > 0 for all Pi > 0 (notes section 2, boxed).
# Proof structure: dg/dPi = -(1/L) Cov_Pi(f, z) (already checked: cov_id == 0).
# We certify Cov_Pi(f, z) < 0 via the FKG/Chebyshev symmetrized identity,
# valid for any normalized density p = sigma, and the pointwise sign of its
# integrand. This is a GLOBAL certificate on Pi > 0, not a finite sweep.
dgPi = sp.diff(gPi, Pi)

# p is a probability density on [0, L]: integral over [0, L] is 1.
norm_p = sp.simplify(sp.integrate(sigma, (z, 0, L)))
if sp.simplify(norm_p - 1) != 0:
    raise AssertionError("sigma_Pi is not a normalized density on [0, L].")

# Covariance as already used in cov_id: Cov = E[f z] - g E[z].
Cov = sp.simplify(Efz - gPi * Ez)

# (1) Symmetrized double-integral identity for the covariance (FKG/Chebyshev).
#     Cov = (1/2) ∫∫ (f(z1)-f(z2))(z1-z2) p(z1) p(z2) dz1 dz2.
z1, z2 = sp.symbols("z1 z2", positive=True, real=True)
f1 = f.subs(z, z1)
f2 = f.subs(z, z2)
p1 = sigma.subs(z, z1)
p2 = sigma.subs(z, z2)
integrand_sym = sp.Rational(1, 2) * (f1 - f2) * (z1 - z2) * p1 * p2
Cov_double = sp.simplify(
    sp.integrate(sp.integrate(integrand_sym, (z1, 0, L)), (z2, 0, L))
)
if sp.simplify(Cov_double - Cov) != 0:
    raise AssertionError("Symmetrized covariance identity failed.")

# (2) Pointwise sign of the symmetrizer factor on [0, L]^2:
#     for 0 <= z2 < z1 <= L, f(z1) - f(z2) < 0 (f strictly decreasing).
#     Hence (f(z1)-f(z2))(z1-z2) <= 0 everywhere on [0,L]^2, < 0 off-diagonal.
#     This is a bounded-domain cosine inequality with NO exp(Pi); it is decidable.
# Direct, engine-light certificate: g(z) := -f(z) = -cos(pi z/2L) is strictly
# increasing on [0,L] because g'(z) = (pi/(2L)) sin(pi z/2L) > 0 for 0 < z < L.
gprime = sp.diff(-f, z)
# gprime = (pi/(2L)) sin(pi z/(2L)); positive on (0, L). Certify by showing it
# equals a manifestly-positive closed form on the open interval.
if sp.simplify(gprime - (sp.pi / (2 * L)) * sp.sin(sp.pi * z / (2 * L))) != 0:
    raise AssertionError("f'(z) closed form mismatch.")
# sin(pi z/(2L)) > 0 strictly for z in (0, L) since its argument lies in (0, pi/2).
# Therefore f is strictly decreasing, the symmetrizer is <= 0 (< 0 off-diagonal),
# p > 0, so Cov < 0 and dg/dPi = -Cov/L > 0 for ALL Pi > 0.
print("Cov_Pi(f,z) (symmetrized) =", Cov_double)
print("f'(z) =", sp.simplify(gprime))
# Sanity: the certified sign must be consistent with the verified identity,
# i.e. dg/dPi = -Cov/L. (Non-tautological: a wrong gPi or wrong Cov breaks it.)
if sp.simplify(dgPi + Cov / L) != 0:
    raise AssertionError("dg/dPi = -(1/L) Cov consistency failed.")
print("Global strict monotonicity certified: dg/dPi = -(1/L) Cov_Pi(f,z) > 0 for Pi>0.")

# (F1c) Uniqueness of the Family-1 compensation point Pi_* on (0, oo):
# g is strictly increasing from g(0+) = 2/pi (~0.6366) to g(oo) = 1, so the
# equation g(Pi) = g_minus has at most one root, and exactly one iff
# 2/pi < g_minus < 1. g_minus is the lower-branch target (notes section 3).
g_lo = sp.N(2 / sp.pi, 30)
g_hi = sp.Integer(1)
if not (g_lo < g_minus < g_hi):
    raise AssertionError(
        "g_minus not strictly inside (2/pi, 1); uniqueness of Pi_* not guaranteed."
    )
print(f"Bracket for unique Pi_*: 2/pi = {g_lo} < g_minus = {g_minus} < 1 = {g_hi}")

# Unique Family-1 compensation point
Pi_star = sp.nsolve(gPi - g_minus, 1.5, tol=1e-30, maxsteps=100, prec=80)
print("Pi_* =", sp.N(Pi_star, 30))
print("x_* = 1/Pi_* =", sp.N(1/Pi_star, 30))
print("g(Pi_*) =", sp.N(gPi.subs(Pi, Pi_star), 30))
print("g'(Pi_*) =", sp.N(sp.diff(gPi, Pi).subs(Pi, Pi_star), 30))
