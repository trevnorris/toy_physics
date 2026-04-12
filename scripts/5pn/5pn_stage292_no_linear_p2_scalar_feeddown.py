#!/usr/bin/env python3
from __future__ import annotations

import sympy as sp
from fivepn_stage289_292_common import *

"""
Stage 292 — no linear grouped-P2 feed-down into the scalar off-bundle slippages.
"""

banner("STAGE 292 — NO LINEAR GROUPED-P2 SCALAR FEED-DOWN")

th, ph = sp.symbols("theta phi", real=True)
measure = sp.sin(th)

# Real normalized l=2 harmonics.
Y20 = sp.sqrt(5 / (16 * sp.pi)) * (3 * sp.cos(th)**2 - 1)
Y21c = sp.sqrt(15 / (4 * sp.pi)) * sp.sin(th) * sp.cos(th) * sp.cos(ph)
Y21s = sp.sqrt(15 / (4 * sp.pi)) * sp.sin(th) * sp.cos(th) * sp.sin(ph)
Y22c = sp.sqrt(15 / (16 * sp.pi)) * sp.sin(th)**2 * sp.cos(2 * ph)
Y22s = sp.sqrt(15 / (16 * sp.pi)) * sp.sin(th)**2 * sp.sin(2 * ph)

subbanner("I. Exact sphere averages vanish")
for name, Y in [("Y20", Y20), ("Y21c", Y21c), ("Y21s", Y21s), ("Y22c", Y22c), ("Y22s", Y22s)]:
    avg = sp.simplify(sp.integrate(sp.integrate(Y * measure, (ph, 0, 2 * sp.pi)), (th, 0, sp.pi)))
    print(f"< {name} > =")
    sp.pprint(avg)
    expect_zero(f"sphere average of {name}", avg)

subbanner("II. Weak-axisymmetric grouped signature")
eps, x0, x1, y0, y1 = sp.symbols("eps x0 x1 y0 y1", real=True)
X20 = x0 + eps * x1
X21 = x0 + eps * x1 / 2
X22 = x0 - eps * x1
Y20v = y0 + eps * y1
Y21v = y0 + eps * y1 / 2
Y22v = y0 - eps * y1

gX = grouped_trace_anomaly(X20, X21, X22)
gY = grouped_trace_anomaly(Y20v, Y21v, Y22v)
Ixy = grouped_bilinear(X20, X21, X22, Y20v, Y21v, Y22v)

print("a_x =")
sp.pprint(gX["a"])
print("b_x =")
sp.pprint(gX["b"])
expect_zero("b_x - 3 a_x", sp.simplify(gX["b"] - 3 * gX["a"]))
print("I[X,Y] =")
sp.pprint(Ixy)
expect_zero("I[X,Y] - 7 eps^2 x1 y1 / 10", sp.simplify(Ixy - sp.Rational(7, 10) * eps**2 * x1 * y1))

subbanner("III. No linear scalar feed-down theorem")
print("Every rotational scalar observable extracted from a pure l=2 perturbation has")
print("vanishing first variation because every real l=2 harmonic has zero sphere average.")
print("Therefore, for pure weak grouped-P2 anisotropy,")
print("  eps_L^(1,P2) = eps_v^(1,P2) = eps_T^(1,P2) = 0,")
print("  hence eps_perp^(1,P2) = delta_perp^(1,P2) = 0.")

subbanner("IV. First nonzero scalar feed-down is quadratic")
print("The first scalar feed-down from grouped anisotropy must therefore enter through")
print("the quadratic grouped bilinears I[X,Y], not linearly in the grouped amplitudes.")
print("So the remaining linear grouped bottleneck is only in the direct outlet")
print("coefficients delta kappa_W and delta gamma_W, not in the scalar slippage channel.")
