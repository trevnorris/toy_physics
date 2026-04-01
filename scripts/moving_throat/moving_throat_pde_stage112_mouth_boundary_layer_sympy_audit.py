#!/usr/bin/env python3
"""
Stage 112 SymPy audit — explicit GNLS + localized-Maxwell mouth boundary layer.
"""

from __future__ import annotations
import sympy as sp

z, L, Pi = sp.symbols("z L Pi", positive=True, real=True)
Theta, V1, M = sp.symbols("Theta V1 M", positive=True, real=True)
sigma_star = sp.symbols("sigma_star", positive=True, real=True)

sigma = Pi*sp.exp(-Pi*z/L)/(L*(1-sp.exp(-Pi)))
J = -M*(Theta*sp.diff(sigma, z) + V1*sigma)

print("sigma_Pi(z) =", sp.simplify(sigma))
print("Normalization =", sp.simplify(sp.integrate(sigma, (z, 0, L))))

J_sub = sp.simplify(J.subs({V1: Pi*Theta/L}))
print("Zero-flux current J_sigma =", J_sub)

mu = Theta*sp.log(sigma/sigma_star) + V1*z
res = sp.simplify(Theta*sp.diff(sigma, z) + V1*sigma)
print("Stationary zero-flux ODE residual =", res.subs({V1: Pi*Theta/L}))

if sp.simplify(sp.integrate(sigma, (z, 0, L)) - 1) != 0:
    raise AssertionError("Profile not normalized.")
if sp.simplify(J_sub) != 0:
    raise AssertionError("Current is not zero on the stationary branch.")
if sp.simplify(res.subs({V1: Pi*Theta/L})) != 0:
    raise AssertionError("Boundary-layer ODE not satisfied.")

print("\nDerived electrochemical bias:")
print("Pi_m = V1*L/Theta")
