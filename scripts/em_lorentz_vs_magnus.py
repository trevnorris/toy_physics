"""
em_lorentz_vs_magnus.py

Purpose:
    Compare the classical Magnus force on a vortex/flux tube in a
    superfluid with a Lorentz force q(E + u x B) in a simple setup,
    to make the mapping "charge = circulation × area" explicit.

Setup:
    - Straight vortex line along the z-axis with circulation Γ.
    - Background flow v_inf in the x-direction.
    - Defect core moving with velocity u relative to the lab frame.
    - Mass density rho0, core radius a (so mouth area ~ pi a^2).

We show:
    F_Magnus = rho0 * Gamma * k × (u - v_inf)
    F_Lorentz = q (E + u × B)

In this script we focus on the *magnetic-like* part:

    F_Magnus,mag ~ rho0 Γ k × u
    F_Lorentz,mag = q u × B

and choose an effective B such that these match when
    q = rho0 * pi * a**2 * Gamma.

The "electric-like" part from pressure/enthalpy gradients is left in
symbolic form; the script is mainly to support the algebraic mapping
between Magnus and the q u × B term.
"""

from sympy import symbols, pi, simplify, Eq, pprint
from sympy.vector import CoordSys3D

# ----------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------

rho0, a, Gamma = symbols('rho0 a Gamma', positive=True)
# Effective charge from toy model:
q = symbols('q', real=True)

# Effective magnetic field amplitude along z
B0 = symbols('B0', real=True)

# Coordinate system and basis
N = CoordSys3D('N')
i, j, k = N.i, N.j, N.k

# Velocities
ux, uy, uz = symbols('ux uy uz', real=True)
u_vec = ux*i + uy*j + uz*k

v0 = symbols('v0', real=True)
v_inf = v0 * i  # background flow in +x

# ----------------------------------------------------------------------
# 1. Magnus force on a straight vortex line (per unit length)
# ----------------------------------------------------------------------

# Classical Magnus force: F_M = rho0 * Gamma * k × (u - v_inf)
F_M = rho0 * Gamma * k.cross(u_vec - v_inf)

print("Magnus force per unit length, F_M = rho0 * Gamma * k × (u - v_inf):")
pprint(F_M)
print()

# Split into u-dependent (magnetic-like) and v_inf-dependent (background) parts
F_M_u_part = simplify(rho0 * Gamma * k.cross(u_vec))
F_M_vinf_part = simplify(-rho0 * Gamma * k.cross(v_inf))

print("u-dependent (magnetic-like) part of Magnus force:")
pprint(F_M_u_part)
print()

print("Background-flow (v_inf) contribution to Magnus force:")
pprint(F_M_vinf_part)
print()


# ----------------------------------------------------------------------
# 2. Lorentz force with q = rho0 * pi * a^2 * Gamma
# ----------------------------------------------------------------------

q_def = rho0 * pi * a**2 * Gamma
print("Toy-model effective charge: q_defect = rho0 * pi * a^2 * Gamma")
pprint(Eq(q, q_def))
print()

# Effective magnetic field B along z-axis:
B_vec = B0 * k

# Magnetic part of Lorentz force: F_L,mag = q * u × B
F_L_mag = simplify(q * u_vec.cross(B_vec))

print("Magnetic part of Lorentz force, F_L,mag = q * u × B:")
pprint(F_L_mag)
print()

# Substitute q_def into Lorentz force
F_L_mag_qdef = simplify(F_L_mag.subs(q, q_def))
print("F_L,mag with q = rho0 * pi * a^2 * Gamma:")
pprint(F_L_mag_qdef)
print()


# ----------------------------------------------------------------------
# 3. Compare with Magnus u-part
# ----------------------------------------------------------------------

F_M_u_simplified = simplify(F_M_u_part)
F_L_u_simplified = simplify(F_L_mag_qdef)

print("u-dependent Magnus term:")
pprint(F_M_u_simplified)
print()

print("u-dependent Lorentz term (with q_def):")
pprint(F_L_u_simplified)
print()

print(r"""
We can read off the proportionality by inspection:
    F_M_u = rho0 * Gamma * (k × u)
           = -rho0 * Gamma * (u × k)

    F_L_u = q_def * B0 * (u × k)
          = rho0 * pi * a^2 * Gamma * B0 * (u × k)

So matching F_L_u = F_M_u gives
    -rho0 * Gamma = rho0 * pi * a^2 * Gamma * B0
    => B0 = -1 / (pi * a^2).

Thus, with the identification
    q = rho0 * pi * a^2 * Gamma,
    B = -(1 / (pi * a^2)) k,
the u-dependent Magnus force matches the magnetic part of the Lorentz force.

The remaining background-flow term F_M_vinf can be grouped with
pressure-gradient forces and interpreted as part of an effective
electric field E in the Lorentz force qE.
""")

if __name__ == "__main__":
    print("em_lorentz_vs_magnus.py: algebraic comparison complete.")
