#!/usr/bin/env python3
"""
Step 5 SymPy audit: full 3D dipole-orientation law and finite-size expansion
in the coaxial far-field limit.

The dipole law is built from actual vectors m_A = pi R_A^2 I_A s_A. The
finite-size expansion is treated as an asymptotic far-field expansion, not as a
near-contact convergence proof.

Run:
    python step_05_3d_orientation_and_finite_size_sympy.py
"""
import sympy as sp

mu0, R1, R2, d, I1, I2 = sp.symbols("mu0 R1 R2 d I1 I2", positive=True)
u = sp.symbols("u", real=True)

ex = sp.Matrix([1, 0, 0])
ey = sp.Matrix([0, 1, 0])
ez = sp.Matrix([0, 0, 1])


def dipole_force(s1_hat, s2_hat, dhat):
    """Fixed-current radial force from vector dipoles."""
    m1_vec = sp.pi * R1**2 * I1 * s1_hat
    m2_vec = sp.pi * R2**2 * I2 * s2_hat
    orientation = sp.simplify(
        3 * m1_vec.dot(dhat) * m2_vec.dot(dhat) - m1_vec.dot(m2_vec)
    )
    potential = sp.simplify(-mu0 * orientation / (4 * sp.pi * d**3))
    force = sp.simplify(-sp.diff(potential, d))
    return m1_vec, m2_vec, orientation, potential, force


cases = {
    "coaxial aligned": (ez, ez, ez),
    "coaxial anti-aligned": (ez, -ez, ez),
    "side-by-side parallel": (ex, ex, ez),
    "side-by-side anti-parallel": (ex, -ex, ez),
    "one axial one transverse": (ez, ex, ez),
}
case_data = {name: dipole_force(*vectors) for name, vectors in cases.items()}
case_forces = {name: data[-1] for name, data in case_data.items()}

step3_leading_force = -3 * mu0 * sp.pi * I1 * I2 * R1**2 * R2**2 / (2 * d**4)

# Coaxial finite-size expansion through the next asymptotic term. This repeats
# the Step-3 Neumann expansion locally so the finite-size statement is audited
# rather than imported.
A = R1**2 + R2**2
B = 2 * R1 * R2
eps = (A - B * sp.cos(u)) / d**2
inv_dist_series_next = (1/d) * (
    1 - eps/2 + sp.Rational(3, 8) * eps**2 - sp.Rational(5, 16) * eps**3
)
integral_series_next = sp.integrate(sp.cos(u) * inv_dist_series_next, (u, 0, 2 * sp.pi))
M_series_next = sp.simplify(mu0 * R1 * R2 / 2 * integral_series_next)
F_series_next = sp.simplify(I1 * I2 * sp.diff(M_series_next, d))
finite_ratio_next = sp.simplify(F_series_next / step3_leading_force)
finite_ratio_expected = sp.simplify(
    1
    - sp.Rational(5, 2) * (R1**2 + R2**2) / d**2
    + sp.Rational(35, 8) * (R1**4 + 3 * R1**2 * R2**2 + R2**4) / d**4
)
first_small_parameter = (R1**2 + R2**2) / d**2

assert sp.simplify(case_forces["coaxial aligned"] - step3_leading_force) == 0
assert sp.simplify(case_forces["coaxial anti-aligned"] + step3_leading_force) == 0
assert sp.simplify(case_forces["side-by-side parallel"] - 3*mu0*sp.pi*I1*I2*R1**2*R2**2/(4*d**4)) == 0
assert sp.simplify(case_forces["side-by-side anti-parallel"] + 3*mu0*sp.pi*I1*I2*R1**2*R2**2/(4*d**4)) == 0
assert case_forces["one axial one transverse"] == 0
assert sp.simplify(finite_ratio_next - finite_ratio_expected) == 0

print("STEP 5: 3D VECTOR DIPOLE LAW AND FINITE-SIZE EXPANSION")
print("Dipoles are m_A = pi R_A^2 I_A s_A.")
print("Fixed-current potential:")
print("U = -mu0/(4*pi*d^3) * [3(m1.dhat)(m2.dhat) - m1.m2]")
print("F_d = -dU/dd, with F_d < 0 attraction")
print()
for name, (_, _, orientation, potential, force) in case_data.items():
    if force == 0:
        sign = "zero at leading dipole order"
    elif force.could_extract_minus_sign():
        sign = "attraction"
    else:
        sign = "repulsion"
    print(f"{name:28s}:")
    print("  orientation numerator =", orientation)
    print("  U =", potential)
    print("  F_d =", force, " ->", sign)
print()
print("Coaxial aligned force exactly matches Step-3 leading term:")
print("F_leading =", step3_leading_force)
print()
print("Coaxial finite-size force ratio through the next asymptotic term:")
print("F/F_leading =", finite_ratio_next)
print("small far-field parameter =", first_small_parameter)
print("This is an asymptotic d >> R_A check, not a near-contact convergence proof.")
