"""
em_continuity_and_lorenz.py

Purpose:
    Current-supportive symbolic check for the 4d_em_fields paper.

    In flat 3+1 form, verify that
        Box(phi) = -mu0 c^2 rho_e,
        Box(A_i) = -mu0 J_i,
    together with the Lorenz gauge condition
        div(A) + (1/c^2) dphi/dt = 0,
    imply charge continuity
        d(rho_e)/dt + div(J) = 0.

Dependencies:
    sympy
"""

from sympy import Eq, Function, diff, pprint, simplify, symbols


def box(field, time_coord, spatial_coords, light_speed):
    x, y, z = spatial_coords
    return diff(field, time_coord, 2) / light_speed**2 - (
        diff(field, x, 2) + diff(field, y, 2) + diff(field, z, 2)
    )


t, x, y, z = symbols("t x y z", real=True)
c, mu0 = symbols("c mu0", positive=True)

phi = Function("phi")(t, x, y, z)
Ax = Function("Ax")(t, x, y, z)
Ay = Function("Ay")(t, x, y, z)
Az = Function("Az")(t, x, y, z)

rho_e = Function("rho_e")(t, x, y, z)
Jx = Function("Jx")(t, x, y, z)
Jy = Function("Jy")(t, x, y, z)
Jz = Function("Jz")(t, x, y, z)

lorenz_gauge = diff(phi, t) / c**2 + diff(Ax, x) + diff(Ay, y) + diff(Az, z)

eq_phi = Eq(box(phi, t, (x, y, z), c), -mu0 * c**2 * rho_e)
eq_Ax = Eq(box(Ax, t, (x, y, z), c), -mu0 * Jx)
eq_Ay = Eq(box(Ay, t, (x, y, z), c), -mu0 * Jy)
eq_Az = Eq(box(Az, t, (x, y, z), c), -mu0 * Jz)

continuity = diff(rho_e, t) + diff(Jx, x) + diff(Jy, y) + diff(Jz, z)
box_divA = box(lorenz_gauge, t, (x, y, z), c)
continuity_combo = simplify(diff(eq_phi.lhs - eq_phi.rhs, t) / c**2 + (
    diff(eq_Ax.lhs - eq_Ax.rhs, x)
    + diff(eq_Ay.lhs - eq_Ay.rhs, y)
    + diff(eq_Az.lhs - eq_Az.rhs, z)
))
identity_residual = simplify(continuity_combo - box_divA - mu0 * continuity)

print("Lorenz-gauge continuity check")
print("=============================")
print()

print("Lorenz gauge condition:")
pprint(Eq(lorenz_gauge, 0))
print()

print("Wave equations:")
pprint(eq_phi)
pprint(eq_Ax)
pprint(eq_Ay)
pprint(eq_Az)
print()

print("Check exact identity:")
print("div(field equations) - Box(Lorenz gauge) - mu0 * continuity =")
pprint(identity_residual)
print("PASS" if identity_residual == 0 else "FAIL")
print()

print("Conclusion:")
print(
    "The exact identity is div(field equations) = Box(Lorenz gauge) + mu0 * continuity. "
    "Under the Lorenz gauge, continuity therefore follows."
)
