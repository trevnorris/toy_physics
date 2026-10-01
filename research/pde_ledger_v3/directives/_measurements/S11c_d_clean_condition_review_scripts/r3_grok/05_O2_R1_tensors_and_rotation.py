#!/usr/bin/env python3
"""Does R1 as written (polar vectors along dir 1; tensors rotation-invariant
about dir 1) imply the full O(2) A uses for the oblique-incidence claim?

Also: which rotation takes an oblique (k2, k3) Fourier mode to class P
(independent of direction 3)? About ŵ (interface normal) or about ê1
(background-gradient axis)?
"""
from __future__ import annotations

import sympy as sp

# SO(2) about axis 1, acting on components (1,2,3):
#  v1' = v1
#  v2' =  c v2 - s v3
#  v3' =  s v2 + c v3
phi = sp.symbols("phi")
c, s = sp.cos(phi), sp.sin(phi)
R = sp.Matrix([[1, 0, 0], [0, c, -s], [0, s, c]])
print("SO2_ABOUT_DIR1", R)

# Polar vector along dir 1: V = (V1, 0, 0)
V = sp.Matrix([sp.symbols("V1"), 0, 0])
print("POLAR_ALONG_1_ROTATED", sp.simplify(R * V))
# Reflection of axis 3: M = diag(1,1,-1)
M = sp.diag(1, 1, -1)
print("POLAR_ALONG_1_REFLECTED", M * V)

# Polar vector along dir 2 (Codex R2-2 counterexample)
V2 = sp.Matrix([0, sp.symbols("V2"), 0])
print("POLAR_ALONG_2_ROTATED", sp.simplify(R * V2))
print("POLAR_ALONG_2_REFLECTED", M * V2)
print("POLAR_ALONG_2_BREAKS_SO2", sp.simplify(R * V2 - V2) != sp.Matrix([0, 0, 0]))
print("POLAR_ALONG_2_BREAKS_MIRROR", M * V2 != V2)

# Axial vector along dir 1 (or equivalently the area tensor in the 2-3 plane)
A = sp.Matrix([sp.symbols("A1"), 0, 0])
# Axial vectors transform with det(R): A' = det(R) R A. SO(2) has det=+1
print("AXIAL_ALONG_1_SO2", sp.simplify(R * A))
print("AXIAL_ALONG_1_MIRROR_DET", M.det())
print("AXIAL_ALONG_1_AFTER_MIRROR", M.det() * M * A)
print("AXIAL_ALONG_1_BREAKS_MIRROR", M.det() * M * A != A)

# Rank-2 tensor SO(2)-invariant about axis 1: T11, T22=T33, T23=-T32, T12=T13=0
a, b, d = sp.symbols("T11 T_perp T23")
T = sp.Matrix([[a, 0, 0], [0, b, d], [0, -d, b]])
print("T_SO2_INVARIANT_ANSATZ", T)
Trot = sp.simplify(R * T * R.T)
print("T_ROTATED_MINUS_T", sp.simplify(Trot - T))
Tref = M * T * M.T  # polar 2-tensor
print("T_REFLECTED", Tref)
print("T_MIRROR_RESIDUAL", sp.simplify(Tref - T))
print("ANTISYMMETRIC_T23_BREAKS_MIRROR", sp.simplify((Tref - T)[1, 2]) != 0)

# Energy mixing from T23: (1/2) u_i T_ij u_j contains T23 (u2 u3 - u3 u2)=0
# because the antisymmetric part drops from a symmetric bilinear in u.
# Derivative coupling: T_ij (∂_i θ) u_j  or T_ij (curl-type).
# The invariant axial coupling is A · (∇×u) = A1 (∂2 u3 - ∂3 u2).
u2, u3, d2u3, d3u2, theta = sp.symbols("u2 u3 d2u3 d3u2 theta")
mix = A[0] * (d2u3 - d3u2) * theta
print("AXIAL_MIXING_TERM", mix)
# On class P, d3u2=0: A1 * d2u3 * theta  -- odd-even mixing
print("AXIAL_MIXING_ON_P", mix.subs(d3u2, 0))

# Symmetric T22=T33 does not mix u2 with u3 under a mass term T_perp (u2^2+u3^2).
# It is O(2) even.

# Rotation that maps k = (0, k2, k3) to (0, k, 0):
k2, k3 = sp.symbols("k2 k3", real=True)
kvec = sp.Matrix([0, k2, k3])
# Want R k = (0, k_perp, 0). Choose phi so that
# c k2 - s k3 = k_perp, s k2 + c k3 = 0  => tan phi = -k3/k2
# This R is SO(2) about direction 1, NOT a rotation about ŵ.
# Rotation about ŵ (the slab normal, orthogonal to all of 1,2,3) is an
# in-plane SO(3) action that can mix direction 1 with 2 and 3, which would
# rotate the background gradient off axis 1.
print("K_INPLANE", kvec)
print("ROTATION_THAT_KILLS_K3_IS_ABOUT_DIR1", True)
print("ROTATION_ABOUT_W_MIXES_DIR1_INTO_23", True)

# Background W(x1) after a rotation about ŵ that rotates axis 1 toward axis 2:
# W(x1) becomes W(c x1 + s x2), which depends on x2, breaking R1 and the
# x3-mirror of the original coordinates.
print("W_OF_X1_NOT_INVARIANT_UNDER_ROTATION_ABOUT_W", True)
print("W_OF_X1_INVARIANT_UNDER_SO2_ABOUT_DIR1", True)
