#!/usr/bin/env python3
"""Inversion parities of toroidal vs spheroidal polar displacements vs scalars.

Active inversion on a polar vector: (P u)(r) = -u(-r).
On a scalar: (P φ)(r) = φ(-r).
Toroidal: u_T = r × ∇ Y_ℓm  (tangent to r=const, div-free).
Spheroidal radial: u_r = r̂ Y_ℓm.
Spheroidal longitudinal-horizontal: u_L = r ∇_⊥ Y_ℓm.
"""
from __future__ import annotations

import sympy as sp

x, y, z = sp.symbols("x y z", real=True)
r2 = x * x + y * y + z * z

def inversion_scalar(expr):
    return expr.subs({x: -x, y: -y, z: -z}, simultaneous=True)

def inversion_polar_vector(vec):
    """(P u)(r) = - u(-r)."""
    uminus = [inversion_scalar(c) for c in vec]
    return [-c for c in uminus]


def check_eigen(name, obj, expected, is_vector=False):
    if is_vector:
        Pu = inversion_polar_vector(obj)
        residual = [sp.expand(sp.simplify(Pu[i] - expected * obj[i])) for i in range(3)]
        ok = all(r == 0 for r in residual)
        print(f"{name}_EXPECTED", expected)
        print(f"{name}_RESIDUAL", residual)
        print(f"{name}_OK", ok)
    else:
        Pphi = inversion_scalar(obj)
        residual = sp.expand(sp.simplify(Pphi - expected * obj))
        print(f"{name}_EXPECTED", expected)
        print(f"{name}_RESIDUAL", residual)
        print(f"{name}_OK", residual == 0)


# ℓ=1, m=0: Y10 ∝ z / r  (use Cartesian harmonic Φ = z, harmonic homogeneous)
# For inversion tests the radial factor r^ℓ Y_ℓm is a harmonic polynomial.
# ℓ=1: Φ1 = z. Scalar parity (-1)^1 = -1.
Phi1 = z
check_eigen("SCALAR_L1", Phi1, -1, is_vector=False)

# Toroidal ℓ=1: r × ∇Φ1. ∇Φ1 = (0,0,1), r×∇Φ1 = (y, -x, 0)
uT1 = [y, -x, sp.Integer(0)]
# div: ∂x y + ∂y (-x) + 0 = 0
# r·u = x y - y x = 0
print("TOROIDAL_L1_DIV", sp.diff(uT1[0], x) + sp.diff(uT1[1], y) + sp.diff(uT1[2], z))
print("TOROIDAL_L1_RADIAL", sp.expand(x * uT1[0] + y * uT1[1] + z * uT1[2]))
check_eigen("TOROIDAL_L1", uT1, (-1) ** (1 + 1), is_vector=True)  # (-1)^{ℓ+1} = +1

# Spheroidal radial ℓ=1: r̂ Φ1, use u = r Φ1 / something; Cartesian r̂~r, take u = (x,y,z)*z? 
# Cleaner: u_r ∝ Φ1 * (x,y,z) / r^2 * r = Φ1 r̂. Use u = Φ1 * (x,y,z) = z (x,y,z)
uS1 = [z * x, z * y, z * z]
check_eigen("SPHEROIDAL_RADIAL_L1", uS1, (-1) ** 1, is_vector=True)  # (-1)^ℓ = -1

# Spheroidal horizontal: r^2 ∇Φ1 - r (r·∇Φ1) r̂-direction... ∇_⊥ Φ1.
# For Φ1=z, ∇Φ1=(0,0,1), radial projection (z/r^2)(x,y,z), 
# ∇_⊥ Φ1 = (0,0,1) - z (x,y,z)/r^2. Multiply by r^2: u = (0,0,r2) - z(x,y,z)
uH1 = [sp.Integer(0) - z * x, 0 - z * y, r2 - z * z]
check_eigen("SPHEROIDAL_HORIZ_L1", uH1, (-1) ** 1, is_vector=True)

# ℓ=2, m=0: Φ2 = 3z^2 - r^2  (∝ r^2 Y20), scalar parity (+1)
Phi2 = 3 * z * z - r2
check_eigen("SCALAR_L2", Phi2, 1, is_vector=False)

# Toroidal ℓ=2: r × ∇Φ2
gPhi2 = [sp.diff(Phi2, x), sp.diff(Phi2, y), sp.diff(Phi2, z)]
# r × g = (y gz - z gy, z gx - x gz, x gy - y gx)
uT2 = [
    y * gPhi2[2] - z * gPhi2[1],
    z * gPhi2[0] - x * gPhi2[2],
    x * gPhi2[1] - y * gPhi2[0],
]
uT2 = [sp.expand(c) for c in uT2]
print("TOROIDAL_L2_DIV", sp.expand(sp.diff(uT2[0], x) + sp.diff(uT2[1], y) + sp.diff(uT2[2], z)))
print("TOROIDAL_L2_RADIAL", sp.expand(x * uT2[0] + y * uT2[1] + z * uT2[2]))
check_eigen("TOROIDAL_L2", uT2, (-1) ** (2 + 1), is_vector=True)  # -1

uS2 = [sp.expand(Phi2 * x), sp.expand(Phi2 * y), sp.expand(Phi2 * z)]
check_eigen("SPHEROIDAL_RADIAL_L2", uS2, (-1) ** 2, is_vector=True)

# Mixing: a polar-vector inner product u_T · u_S is a scalar. For it to be
# inversion-even as an energy density, the two vectors must have the same π.
# If π_T ≠ π_S they cannot appear bilinearly in a true scalar without an extra
# pseudoscalar (Levi-Civita / axial background).
print("L1_PI_TOROIDAL", 1)
print("L1_PI_SPHEROIDAL", -1)
print("L1_PI_SCALAR", -1)
print("L2_PI_TOROIDAL", -1)
print("L2_PI_SPHEROIDAL", 1)
print("L2_PI_SCALAR", 1)
print("SAME_L_TOROIDAL_NE_SCALAR", True)
print("SAME_L_SPHEROIDAL_EQ_SCALAR", True)
