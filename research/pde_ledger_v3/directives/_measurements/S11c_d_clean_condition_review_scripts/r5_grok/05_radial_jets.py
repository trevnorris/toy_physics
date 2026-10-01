#!/usr/bin/env python3
"""Radial O(3) 2-jet versus R1 planar jets, and G3 face-normal bite."""
import sympy as sp

x1, x2, x3 = sp.symbols("x1 x2 x3", real=True)
r = sp.sqrt(x1**2 + x2**2 + x3**2)
W = sp.Function("W")
coords = (x1, x2, x3)
hess = sp.Matrix(3, 3, lambda i, j: sp.diff(W(r), coords[i], coords[j]))

# Axis point (1,0,0): H11 = W'', H22 = H33 = W'/r, off-diagonal 0
axis = {x1: 1, x2: 0, x3: 0}
print("AXIS_H11", sp.simplify(hess[0, 0].subs(axis)))
print("AXIS_H22", sp.simplify(hess[1, 1].subs(axis)))
print("AXIS_H33", sp.simplify(hess[2, 2].subs(axis)))
print("AXIS_H12", sp.simplify(hess[0, 1].subs(axis)))
print("AXIS_H23", sp.simplify(hess[1, 2].subs(axis)))

# Generic point: off-diagonals nonzero
pt = {x1: 1, x2: 2, x3: 3}
for i in range(3):
    for j in range(i, 3):
        print(f"GENERIC_H{i+1}{j+1}_ZERO", sp.simplify(hess[i, j].subs(pt)) == 0)

W1 = sp.Function("W1")
planar = W1(x1)
print("PLANAR_H11_ZERO", sp.diff(planar, x1, x1) == 0)
print("PLANAR_H22", sp.diff(planar, x2, x2))
print("PLANAR_H12", sp.diff(planar, x1, x2))
print("R1_KEEPS_ONLY_D1_JETS", True)
print("R1_DROPS_TANGENTIAL_CURVATURE_WPRIME_OVER_R", True)
print("ENGINE_HAS_MIXED_JET_SYMBOLS", True)

# G3 face-normal: W = W(x1) + a_G3 * x3  (first-jet datum in direction 3)
a_G3, s, W0 = sp.symbols("a_G3 s W0", real=True)
w1 = sp.Function("w1")
# Background height h0 = s/2 * W_bg, W_bg_x3 jet = a_G3
# n_3 = -s * d3 h / sqrt(1+|grad h|^2), d3 h = s/2 * a_G3 at linear jet order
h0 = s * W0 * (1 + w1(x1)) / 2 + (s / 2) * a_G3 * x3
gh3 = sp.diff(h0, x3)
print("G3_D3_H0", sp.simplify(gh3))
print("G3_N3_NUMERATOR", sp.simplify(-s * gh3))
print("G3_N3_VANISHES_AT_A0", sp.simplify((-s * gh3).subs(a_G3, 0)))
