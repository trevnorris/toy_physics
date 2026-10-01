#!/usr/bin/env python3
"""Kinematic counterexample to an exact 'fast bulk is safe' claim at a smooth defect."""

import sympy as sp

omega = sp.Rational(1)
c_light = sp.Rational(1)
c_bulk = sp.Rational(10)
k_incident = omega / c_light

# Uniform interface: tangential momentum is conserved.
q2_uniform = sp.simplify((omega / c_bulk) ** 2 - k_incident**2)

# A static smooth Gaussian profile has a nonzero Fourier component at every
# finite momentum transfer.  Choose G=-k_incident so the scattered tangential
# momentum is zero.
L = sp.symbols("L", positive=True)
G = -k_incident
gaussian_fourier_carrier = sp.exp(-L**2 * G**2 / 2)
k_scattered = sp.simplify(k_incident + G)
q2_scattered = sp.simplify((omega / c_bulk) ** 2 - k_scattered**2)

print("C_BULK_OVER_C_LIGHT=", c_bulk / c_light)
print("UNIFORM_INTERFACE_Q_SQUARED=", q2_uniform)
print("SMOOTH_GAUSSIAN_TRANSFER_CARRIER_AT_G_MINUS_K=", gaussian_fourier_carrier)
print("SCATTERED_TANGENTIAL_K=", k_scattered)
print("SCATTERED_BULK_Q_SQUARED=", q2_scattered)
print("GAUSSIAN_TRANSFER_IS_IDENTICALLY_ZERO=", sp.simplify(gaussian_fourier_carrier) == 0)
