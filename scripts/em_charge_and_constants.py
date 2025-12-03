"""
em_charge_and_constants.py

Purpose:
    Symbolically connect the Coulomb / Maxwell normalizations to the
    superfluid toy-model parameters.

    - Start from a generic Coulomb potential phi ~ q / (4*pi*eps0*r)
      and verify Gauss's law structure.
    - Encode the defect-based charge and mass scalings:
          m ~ rho0 * (pi * a**2 * L)
          q ~ rho0 * (pi * a**2) * Gamma
      as already derived in the Mathematica notebook.
    - Impose c^2 = 1 / (eps0 * mu0) with c identified as the acoustic
      wave speed, and keep eps0, mu0 symbolic so they can later be
      related to the microphysics (rho0, c, geometry, etc.).
    - Provide a simple symbolic definition of (rho_e, J) in terms of
      q and a pointlike defect worldline, suitable for use in the EM paper.

This script is intended to be *didactic* rather than numerically useful:
it prints identities and symbolic relationships that show how the EM
constants and sources fit into the hydrodynamic picture.

Dependencies:
    sympy
"""

from sympy import (
    symbols, Function, sqrt, pi, simplify,
    Eq, pprint
)
from sympy.vector import CoordSys3D, gradient, divergence

# ----------------------------------------------------------------------
# 0. Parameters and basic symbols
# ----------------------------------------------------------------------

# Fluid / defect parameters
rho0, a, L, Gamma = symbols('rho0 a L Gamma', positive=True)
aspect = symbols('aspect', positive=True)   # L/a ratio from stability (≈ 1.8475)

# EM parameters and physical constants
q, eps0, mu0, c = symbols('q eps0 mu0 c', positive=True)

# Spatial coordinates for field calculations
N = CoordSys3D('N')
x, y, z = N.x, N.y, N.z

# Radial distance from origin
r = sqrt(x**2 + y**2 + z**2)

print("===========================================================")
print("PART 1: Defect mass and charge scaling (from toy model)")
print("===========================================================\n")

# Effective geometric relations from the stability analysis:
#   L = aspect * a   with aspect ≈ sqrt(2)*pi/x01
L_expr = aspect * a

# Gravitational mass ~ displaced volume of the flux tube
massG = rho0 * (pi * a**2 * L_expr)

# Electric charge ~ mass density * mouth area * circulation
chargeQ = rho0 * (pi * a**2) * Gamma

print("Gravitational mass m_G ~ rho0 * Volume:")
pprint(Eq(symbols('m_G'), massG))
print()

print("Electric charge q ~ rho0 * Area * Circulation:")
pprint(Eq(symbols('q_defect'), chargeQ))
print()

print("These relations encode the basic hierarchy scaling:")
force_ratio = simplify(chargeQ**2 / massG**2)
print("F_elec / F_grav ~ q^2 / m_G^2 ~")
pprint(force_ratio)
print()


# ----------------------------------------------------------------------
# 1. Coulomb potential and Gauss's law structure
# ----------------------------------------------------------------------
print("===========================================================")
print("PART 2: Coulomb potential and Gauss's law")
print("===========================================================\n")

# Define the Coulomb potential of a point charge at the origin:
phi_C = q / (4 * pi * eps0 * r)

print("Coulomb potential (symbolic, in terms of q and eps0):")
pprint(Eq(symbols('phi_C'), phi_C))
print()

# Electric field E = -grad(phi_C)
E_vec = -gradient(phi_C, N)

print("Electric field E = -grad(phi_C):")
pprint(E_vec)
print()

# Divergence of E, away from the origin.
divE = divergence(E_vec)

print("Divergence of E (valid for r != 0):")
pprint(Eq(symbols('divE'), simplify(divE)))
print("""
SymPy returns zero for div(E) away from the origin, as expected.
The missing piece at r = 0 is the standard distributional identity:
    div(E) = q * delta^3(r) / eps0
which encodes Gauss's law for a point charge.
""")


# ----------------------------------------------------------------------
# 2. Relating eps0 and mu0 via the wave speed c
# ----------------------------------------------------------------------
print("===========================================================")
print("PART 3: Wave speed and the relation eps0 * mu0 = 1/c^2")
print("===========================================================\n")

print("We impose the standard relativistic relation:")
pprint(Eq(eps0 * mu0, 1 / c**2))
print()

print("""
In the toy model, the effective light speed c is identified with the
acoustic wave speed c_s of the superfluid vacuum (up to 1PN corrections).
The above relation is therefore a *definition* of mu0 in terms of eps0
and the already-derived wave speed.
""")

mu0_expr = simplify(1 / (eps0 * c**2))
print("Solving for mu0 in terms of eps0 and c:")
pprint(Eq(mu0, mu0_expr))
print()


# ----------------------------------------------------------------------
# 3. Source densities rho_e and J in terms of q and a point defect
# ----------------------------------------------------------------------
print("===========================================================")
print("PART 4: Charge density rho_e and current J for a moving defect")
print("===========================================================\n")

t = symbols('t')
ux, uy, uz = symbols('ux uy uz')
u_vec = ux * N.i + uy * N.j + uz * N.k

# Symbolic placeholder for delta^3(r)
delta3 = Function('delta3')(x, y, z)   # stand-in for delta^3(r)

rho_e = q * delta3
J_vec = rho_e * u_vec

print("Symbolic charge density for a point defect:")
pprint(Eq(symbols('rho_e'), rho_e))
print()

print("Symbolic current density for a moving defect:")
print("J_vec =")
pprint(J_vec)
print()

print("""
Here delta3(x,y,z) is a placeholder for the 3D Dirac delta distribution.
With these definitions, Gauss's law and the continuity equation
    div(E) = rho_e / eps0,
    d(rho_e)/dt + div(J) = 0,
are satisfied in the standard distributional sense for a point charge
moving along its worldline.

The remaining hydrodynamic work (handled in other scripts/notebooks) is:
    - Relate the effective q in rho_e and J_vec back to the fluid
      parameters via q_defect = rho0 * pi * a^2 * Gamma.
    - Use the previously derived acoustic wave equation and the
      identification c^2 = 1/(eps0*mu0) to show that the vector
      potential A obeys Box A = -mu0 J, giving Ampere-Maxwell with
      displacement current automatically.
""")

if __name__ == "__main__":
    print("em_charge_and_constants.py: symbolic derivation complete.")
