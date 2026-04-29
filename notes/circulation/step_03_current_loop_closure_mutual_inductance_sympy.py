#!/usr/bin/env python3
"""
Step 3 SymPy audit: fixed-current coaxial loop closure and finite-mouth
far-field mutual inductance expansion.

Run:
    python step_03_current_loop_closure_mutual_inductance_sympy.py
"""
import sympy as sp

# Positive geometric and physical quantities.
mu0, R1, R2, d, I1, I2 = sp.symbols('mu0 R1 R2 d I1 I2', positive=True)
t1, t2, u = sp.symbols('t1 t2 u', real=True)

# Derive the reduced Neumann integral from the double line integral.
loop1 = sp.Matrix([R1 * sp.cos(t1), R1 * sp.sin(t1), 0])
loop2 = sp.Matrix([R2 * sp.cos(t2), R2 * sp.sin(t2), d])
dl1_dt = sp.diff(loop1, t1)
dl2_dt = sp.diff(loop2, t2)

dot_dl = sp.trigsimp(dl1_dt.dot(dl2_dt))
distance_sq = sp.trigsimp((loop1 - loop2).dot(loop1 - loop2))

dot_reduced = sp.trigsimp(dot_dl.subs(t1, u + t2))
distance_sq_reduced = sp.trigsimp(distance_sq.subs(t1, u + t2))
reduced_kernel = sp.cos(u) / sp.sqrt(d**2 + R1**2 + R2**2 - 2 * R1 * R2 * sp.cos(u))
derived_kernel = sp.simplify(dot_reduced / (R1 * R2 * sp.sqrt(distance_sq_reduced)))

A = R1**2 + R2**2
B = 2 * R1 * R2
# Distance denominator: sqrt(d^2 + A - B cos u)
# Large-d expansion: 1/sqrt(d^2 + A - B cos u)
# = 1/d [1 - eps/2 + 3 eps^2/8 - 5 eps^3/16 + O(eps^4)].
eps = (A - B * sp.cos(u)) / d**2
inv_dist_series = (1/d) * (
    1 - eps/2 + sp.Rational(3, 8) * eps**2 - sp.Rational(5, 16) * eps**3
)

integrand_series = sp.cos(u) * inv_dist_series
integral_series = sp.simplify(sp.integrate(integrand_series, (u, 0, 2*sp.pi)))
M_series = sp.simplify(mu0 * R1 * R2 / 2 * integral_series)
M_expected = mu0 * sp.pi * R1**2 * R2**2 / (2*d**3) * (
    1
    - sp.Rational(3, 2)*(R1**2 + R2**2)/d**2
    + sp.Rational(15, 8)*(R1**4 + 3*R1**2*R2**2 + R2**4)/d**4
)

# Fixed-current mechanical potential and radial force. F_d < 0 means attraction.
U_fixed_current = -I1 * I2 * M_series
F_d = sp.simplify(-sp.diff(U_fixed_current, d))
F_expected = sp.simplify(I1 * I2 * sp.diff(M_expected, d))

# Factor the force into the leading sign and finite-size correction.
leading_force = -3 * mu0 * sp.pi * I1 * I2 * R1**2 * R2**2 / (2*d**4)
force_ratio = sp.simplify(F_d / leading_force)
force_ratio_expected = sp.simplify(
    1
    - sp.Rational(5, 2)*(R1**2 + R2**2)/d**2
    + sp.Rational(35, 8)*(R1**4 + 3*R1**2*R2**2 + R2**4)/d**4
)

assert dot_reduced == R1 * R2 * sp.cos(u)
assert sp.simplify(distance_sq_reduced - (d**2 + R1**2 + R2**2 - 2*R1*R2*sp.cos(u))) == 0
assert sp.simplify(derived_kernel - reduced_kernel) == 0
assert sp.simplify(M_series - M_expected) == 0
assert sp.simplify(F_d - F_expected) == 0
assert sp.simplify(force_ratio - force_ratio_expected) == 0

print("STEP 3: FIXED-CURRENT CURRENT-LOOP CLOSURE")
print("Double-line Neumann pieces:")
print("dl1.dt dot dl2.dt =", dot_dl)
print("|r1-r2|^2 =", distance_sq)
print("After u=t1-t2, kernel =", derived_kernel)
print("Translation symmetry turns double integral into 2*pi times the u integral.")
print()
print("Neumann mutual inductance expansion for coaxial loops:")
print("M(d) =", sp.factor(M_series))
print()
print("Fixed-current mechanical potential:")
print("U_I(d) =", sp.factor(U_fixed_current))
print()
print("Radial force F_d = -dU/dd, with F_d < 0 attraction:")
print("F_d =", sp.factor(F_d))
print()
print("Leading force =", leading_force)
print("Finite-size force ratio =", force_ratio)
print("Verdict: for d >> R_A, I1*I2 > 0 attracts; I1*I2 < 0 repels.")
