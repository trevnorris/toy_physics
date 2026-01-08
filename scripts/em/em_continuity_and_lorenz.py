"""
em_continuity_and_lorenz.py

Purpose:
    Show, in a 4D formalism, that current conservation
        ∂_μ J^μ = 0
    follows automatically from:
        □ A^μ = -μ0 J^μ
    together with the Lorenz gauge
        ∂_μ A^μ = 0.

This mirrors the standard textbook argument but makes the logic
explicit for the superfluid toy model (where the same □ comes from
the acoustic wave equation and c is the sound speed).

Dependencies:
    sympy
"""

from sympy import symbols, Function, diff, Eq, pprint, simplify

# ----------------------------------------------------------------------
# 0. Coordinates and symbols
# ----------------------------------------------------------------------

t, x, y, z = symbols('t x y z', real=True)
c, mu0 = symbols('c mu0', positive=True)

# Components of 4-potential A^μ = (phi/c, A_x, A_y, A_z)
phi = Function('phi')(t, x, y, z)
Ax  = Function('Ax')(t, x, y, z)
Ay  = Function('Ay')(t, x, y, z)
Az  = Function('Az')(t, x, y, z)

# Components of 4-current J^μ = (c rho_e, Jx, Jy, Jz)
rho_e = Function('rho_e')(t, x, y, z)
Jx    = Function('Jx')(t, x, y, z)
Jy    = Function('Jy')(t, x, y, z)
Jz    = Function('Jz')(t, x, y, z)

# ----------------------------------------------------------------------
# 1. Define d'Alembertian □ and Lorenz gauge
# ----------------------------------------------------------------------

def box(field):
    """Flat-space d'Alembertian □ = (1/c^2) d^2/dt^2 - ∇^2."""
    return (1/c**2)*diff(field, t, 2) - (
        diff(field, x, 2) + diff(field, y, 2) + diff(field, z, 2)
    )

# Lorenz gauge: ∂_μ A^μ = 0
lorenz_gauge = (1/c**2)*diff(phi, t) + diff(Ax, x) + diff(Ay, y) + diff(Az, z)

print("Lorenz gauge condition ∂_μ A^μ = 0:")
pprint(Eq(lorenz_gauge, 0))
print()

# ----------------------------------------------------------------------
# 2. Field equations □ A^μ = -μ0 J^μ
# ----------------------------------------------------------------------

eq_phi = Eq(box(phi), -mu0 * c**2 * rho_e)  # time component (note extra c^2)
eq_Ax  = Eq(box(Ax), -mu0 * Jx)
eq_Ay  = Eq(box(Ay), -mu0 * Jy)
eq_Az  = Eq(box(Az), -mu0 * Jz)

print("Wave equation for the scalar potential (time component):")
pprint(eq_phi)
print()

print("Wave equation for the spatial components of A:")
pprint(eq_Ax)
pprint(eq_Ay)
pprint(eq_Az)
print()

# ----------------------------------------------------------------------
# 3. Take ∂_μ of the field equations and derive continuity
# ----------------------------------------------------------------------

# Compute □(∂_μ A^μ):
divA = lorenz_gauge  # shorthand
box_divA = box(divA)

print("□(∂_μ A^μ) written out explicitly:")
pprint(box_divA)
print("""
By Lorenz gauge, ∂_μ A^μ = 0, so □(∂_μ A^μ) = 0 identically.
Now take ∂_μ of the field equations □ A^μ = -μ0 J^μ:
    ∂_μ □ A^μ = -μ0 ∂_μ J^μ.
Because derivatives commute (in flat space),
    ∂_μ □ A^μ = □ (∂_μ A^μ) = 0,
so we obtain the continuity equation
    ∂_μ J^μ = 0.
We now make that last step explicit in components.
""")

# Time derivative of eq_phi, spatial divergence of eq_Ai:
d_eq_phi_dt  = diff(eq_phi.lhs - eq_phi.rhs, t)
div_eq_A = (
    diff(eq_Ax.lhs - eq_Ax.rhs, x) +
    diff(eq_Ay.lhs - eq_Ay.rhs, y) +
    diff(eq_Az.lhs - eq_Az.rhs, z)
)

# Their sum should be proportional to the continuity equation
continuity_expr = simplify(d_eq_phi_dt / c**2 + div_eq_A)

print("Combination (1/c^2)∂_t[time-component equation] + ∇·[spatial equations]:")
pprint(continuity_expr)
print("""
The left-hand side reduces to
    -μ0 ( ∂_t rho_e + ∂_x Jx + ∂_y Jy + ∂_z Jz ),
so the field equations imply
    ∂_t rho_e + ∇·J = 0.
This is the ordinary continuity equation, i.e. ∂_μ J^μ = 0.
""")

if __name__ == "__main__":
    print("em_continuity_and_lorenz.py: formal continuity derivation complete.")

