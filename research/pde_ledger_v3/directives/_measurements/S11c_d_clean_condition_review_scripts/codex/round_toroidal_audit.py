#!/usr/bin/env python3
"""A small Cartesian symbolic model of a round-background toroidal displacement."""

import sympy as sp

x, y, z, a = sp.symbols("x y z a", real=True)
R2 = x**2 + y**2 + z**2
Q = R2 + a * R2**2
u = sp.Matrix([-y * (1 + R2), x * (1 + R2), 0])
grad_Q = sp.Matrix([sp.diff(Q, c) for c in (x, y, z)])
div_u = sum(sp.diff(u[i], c) for i, c in enumerate((x, y, z)))

print("ROUND_BACKGROUND_GRADIENT=", grad_Q)
print("TOROIDAL_U_DOT_GRAD_Q=", sp.expand(u.dot(grad_Q)))
print("TOROIDAL_DIVERGENCE=", sp.expand(div_u))
print("TOROIDAL_FACE_NORMAL_VELOCITY_CARRIER=", sp.expand(u.dot(grad_Q)))
print("MATERIAL_ANCHOR_LINEAR_CARRIER_MINUS_U_DOT_GRAD_Q=", sp.expand(-u.dot(grad_Q)))
print("QUADRATIC_TOROIDAL_SCALAR_U_SQUARED=", sp.expand(u.dot(u)))
