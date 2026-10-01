#!/usr/bin/env python3
"""Symbolic checks for the planar rotation claim and the m/angular-momentum claim."""

import sympy as sp

alpha = sp.symbols("alpha", real=True)
ca, sa = sp.cos(alpha), sp.sin(alpha)

# Rotation about the planar-interface normal y.
Ry = sp.Matrix([[ca, 0, sa], [0, 1, 0], [-sa, 0, ca]])
gx, gy, gz = sp.symbols("g_x g_y g_z", real=True)
g = sp.Matrix([gx, gy, gz])
print("ROTATION_ABOUT_Y_MATRIX=", Ry)
print("GENERIC_BACKGROUND_GRADIENT_AFTER_ROTATION=", Ry * g)
print("PLANAR_Q_OF_Y_GRADIENT_AFTER_ROTATION=", sp.simplify(Ry * sp.Matrix([0, gy, 0])))
print(
    "GENERIC_BACKGROUND_FIXED_EXAMPLE_alpha_pi_over_2=",
    (Ry * sp.Matrix([1, 2, 3])).subs(alpha, sp.pi / 2),
)

# A single tangential wavevector can be put in the x direction; two
# non-collinear components cannot both be made z independent by one rotation.
kx, kz = sp.symbols("k_x k_z", real=True)
k = sp.Matrix([kx, 0, kz])
print("SINGLE_TANGENTIAL_K_AFTER_ROTATION=", Ry * k)
k1 = sp.Matrix([1, 0, 0])
k2 = sp.Matrix([0, 0, 1])
k1_at_zeroing_rotation = sp.simplify((Ry * k1).subs(alpha, 0))
k2_same_rotation = sp.simplify((Ry * k2).subs(alpha, 0))
print("TWO_COMPONENT_COUNTEREXAMPLE_K1_AT_ALPHA0=", k1_at_zeroing_rotation)
print("TWO_COMPONENT_COUNTEREXAMPLE_K2_AT_SAME_ALPHA=", k2_same_rotation)

# For a complex scalar representative of a mode, the standard azimuthal
# Noether density is Im(conjugate(psi) d_phi psi).  m != 0 is sufficient for a
# travelling exp(i m phi) pattern, but not for a real standing cos(m phi) one.
m, phi, omega, time = sp.symbols("m phi omega t", real=True)
psi_travel = sp.exp(sp.I * (m * phi - omega * time))
j_travel = sp.simplify(sp.im(sp.conjugate(psi_travel) * sp.diff(psi_travel, phi)))
psi_stand = sp.cos(m * phi) * sp.cos(omega * time)
j_stand = sp.simplify(sp.im(sp.conjugate(psi_stand) * sp.diff(psi_stand, phi)))
print("AZIMUTHAL_NOETHER_DENSITY_TRAVELLING_MODE=", j_travel)
print("AZIMUTHAL_NOETHER_DENSITY_REAL_STANDING_MODE=", j_stand)

# Quadratic products of one toroidal l sector contain even parity and hence can
# source scalar sectors at second order.
ell = sp.symbols("ell", integer=True)
p_toroidal = (-1) ** (ell + 1)
print("TOROIDAL_LINEAR_PARITY=", p_toroidal)
print("TOROIDAL_QUADRATIC_PARITY=", sp.simplify(p_toroidal**2))
