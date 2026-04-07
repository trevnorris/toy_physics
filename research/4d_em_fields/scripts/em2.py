"""
em2.py

Purpose:
    Current-supportive symbolic check for the 4d_em_fields paper.

    Starting from the standard potential definitions
        B = curl(A),
        E = -grad(phi) - dA/dt,
    this script verifies three identities:

        div(B) = 0,
        curl(E) + dB/dt = 0,
        curl(B) - (1/c^2) dE/dt
          = -Box(A) + grad(div(A) + (1/c^2) dphi/dt).

    The last line makes the paper's logic explicit: if the brane gauge field
    satisfies the wave equation Box(A) = -mu0 J and the Lorenz gauge condition
    div(A) + (1/c^2) dphi/dt = 0, then the Ampere-Maxwell law follows
    automatically.

Dependencies:
    sympy
"""

from sympy import Function, diff, simplify, symbols
from sympy.vector import CoordSys3D, curl, divergence, gradient


def scalar_laplacian(field, coords):
    x, y, z = coords
    return diff(field, x, 2) + diff(field, y, 2) + diff(field, z, 2)


N = CoordSys3D("N")
x, y, z = N.x, N.y, N.z
t = symbols("t", real=True)
c = symbols("c", positive=True)

Ax = Function("Ax")(x, y, z, t)
Ay = Function("Ay")(x, y, z, t)
Az = Function("Az")(x, y, z, t)
phi = Function("phi")(x, y, z, t)

A = Ax * N.i + Ay * N.j + Az * N.k
B = curl(A)
E = -gradient(phi) - diff(A, t)

div_B = simplify(divergence(B))
faraday_residual = curl(E) + diff(B, t)

laplacian_A = (
    scalar_laplacian(Ax, (x, y, z)) * N.i
    + scalar_laplacian(Ay, (x, y, z)) * N.j
    + scalar_laplacian(Az, (x, y, z)) * N.k
)
box_A = laplacian_A - diff(A, t, 2) / c**2

gauge_scalar = divergence(A) + diff(phi, t) / c**2
grad_gauge = gradient(gauge_scalar)
ampere_lhs = curl(B) - diff(E, t) / c**2
ampere_residual = ampere_lhs - (-box_A + grad_gauge)

component_names = {"x": N.i, "y": N.j, "z": N.k}
faraday_checks = {
    name: simplify(faraday_residual.dot(direction))
    for name, direction in component_names.items()
}
ampere_checks = {
    name: simplify(ampere_residual.dot(direction))
    for name, direction in component_names.items()
}

print("Potential-based Maxwell identity checks")
print("=======================================")
print()

print("Gauss-for-B residual div(B):")
print(div_B)
print("PASS" if div_B == 0 else "FAIL")
print()

print("Faraday residuals curl(E) + dB/dt:")
for name in ("x", "y", "z"):
    print(f"  {name}: {faraday_checks[name]}")
print(
    "PASS"
    if all(expr == 0 for expr in faraday_checks.values())
    else "FAIL"
)
print()

print("Ampere-Maxwell identity residuals")
print("curl(B) - (1/c^2) dE/dt - [-Box(A) + grad(gauge)] :")
for name in ("x", "y", "z"):
    print(f"  {name}: {ampere_checks[name]}")
print(
    "PASS"
    if all(expr == 0 for expr in ampere_checks.values())
    else "FAIL"
)
print()

print("Conclusion:")
print(
    "If Box(A) = -mu0 J and div(A) + (1/c^2) dphi/dt = 0, "
    "then curl(B) - (1/c^2) dE/dt = mu0 J."
)
