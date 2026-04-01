#!/usr/bin/env python3
"""
Stage 113 SymPy audit — exact mouth-bias map and Family-1 compensation point.
"""

from __future__ import annotations
import sympy as sp

z, L, Pi = sp.symbols("z L Pi", positive=True, real=True)
g_minus = sp.Float("0.758035078944663")

sigma = Pi*sp.exp(-Pi*z/L)/(L*(1-sp.exp(-Pi)))
f = sp.cos(sp.pi*z/(2*L))

gPi = sp.simplify(sp.integrate(sigma*f, (z, 0, L)))
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

# Unique Family-1 compensation point
Pi_star = sp.nsolve(gPi - g_minus, 1.5, tol=1e-30, maxsteps=100, prec=80)
print("Pi_* =", sp.N(Pi_star, 30))
print("x_* = 1/Pi_* =", sp.N(1/Pi_star, 30))
print("g(Pi_*) =", sp.N(gPi.subs(Pi, Pi_star), 30))
print("g'(Pi_*) =", sp.N(sp.diff(gPi, Pi).subs(Pi, Pi_star), 30))
