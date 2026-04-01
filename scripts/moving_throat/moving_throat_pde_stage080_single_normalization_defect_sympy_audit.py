#!/usr/bin/env python3
"""
Stage 80 SymPy audit: the actual isotropic passive/outgoing quadrupole branch
collapses to a single normalization defect.
"""
from __future__ import annotations
import sympy as sp

G, c, c_s, a = sp.symbols('G c c_s a', positive=True, real=True)
K0, Omega = sp.symbols('K0 Omega', positive=True, real=True)

# Actual branch conservative module coefficients.
K2 = sp.simplify(K0 / (4 * Omega**2))
K4 = sp.simplify(K0 / (4 * Omega**4))
Gamma5 = sp.simplify(9 * K2**sp.Rational(5, 2) / K0**sp.Rational(3, 2))
Gamma5_expected = sp.simplify(9 * K0 / (32 * Omega**5))
assert sp.simplify(Gamma5 - Gamma5_expected) == 0

# GR target normalization.
K0_target = sp.simplify(64 * G * Omega**5 / (45 * c**5))
Omega_geom = sp.simplify(3 * c_s / (2 * a))
K0_target_geom = sp.simplify(K0_target.subs(Omega, Omega_geom))
assert sp.simplify(K0_target_geom - 54 * G * c_s**5 / (5 * a**5 * c**5)) == 0

K2_target = sp.simplify(K0_target / (4 * Omega**2))
K4_target = sp.simplify(K0_target / (4 * Omega**4))
Gamma5_target = sp.simplify(9 * K2_target**sp.Rational(5, 2) / K0_target**sp.Rational(3, 2))
assert sp.simplify(Gamma5_target - 2 * G / (5 * c**5)) == 0

# Single normalization defect.
NQ = sp.symbols('N_Q', positive=True, real=True)
subs_actual = {K0: sp.simplify(NQ * K0_target)}
R0 = sp.simplify(K0.subs(subs_actual) / K0_target - 1)
R2 = sp.simplify(K2.subs(subs_actual) / K2_target - 1)
R4 = sp.simplify(K4.subs(subs_actual) / K4_target - 1)
R5 = sp.simplify(Gamma5.subs(subs_actual) / Gamma5_target - 1)
assert sp.simplify(R0 - (NQ - 1)) == 0
assert sp.simplify(R2 - (NQ - 1)) == 0
assert sp.simplify(R4 - (NQ - 1)) == 0
assert sp.simplify(R5 - (NQ - 1)) == 0

print('K2 =', K2)
print('K4 =', K4)
print('Gamma5 =', Gamma5_expected)
print('K0_target =', K0_target)
print('K0_target (Omega=3 c_s/(2a)) =', K0_target_geom)
print('Gamma5_target =', Gamma5_target)
print('R0 =', sp.factor(R0))
print('R2 =', sp.factor(R2))
print('R4 =', sp.factor(R4))
print('R5 =', sp.factor(R5))
print('\nSTAGE 80 AUDIT PASSED')
