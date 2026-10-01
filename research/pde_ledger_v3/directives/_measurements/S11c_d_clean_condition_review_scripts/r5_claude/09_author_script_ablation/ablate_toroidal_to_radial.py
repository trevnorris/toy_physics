#!/usr/bin/env python3
"""Compute the round-sector carrier checks added to packet v5.

The script uses explicit Gaussian-weighted l=1 vector fields.  It computes
helicity, inversion parity, angular momentum, a Coriolis image's scalar
projections, and the O(3)-fixed subspace of a vector carrier.
"""

import sympy as sp

x, y, z, t = sp.symbols("x y z t", real=True)
omega, mixing, Omega = sp.symbols("omega mixing Omega", real=True)
coordinates = (x, y, z)
r_vector = [x, y, z]
weight = sp.exp(-(x**2 + y**2 + z**2))


def curl(vector):
    return [
        sp.diff(vector[2], y) - sp.diff(vector[1], z),
        sp.diff(vector[0], z) - sp.diff(vector[2], x),
        sp.diff(vector[1], x) - sp.diff(vector[0], y),
    ]


def divergence(vector):
    return sum(sp.diff(vector[index], coordinates[index]) for index in range(3))


def cross(left, right):
    return [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]


def integrate_all(expression):
    return sp.simplify(
        sp.integrate(
            sp.expand(expression),
            (x, -sp.oo, sp.oo),
            (y, -sp.oo, sp.oo),
            (z, -sp.oo, sp.oo),
        )
    )


def inversion_parity(vector):
    image = [
        sp.expand(-component.subs({x: -x, y: -y, z: -z}, simultaneous=True))
        for component in vector
    ]
    if all(sp.simplify(image[index] - vector[index]) == 0 for index in range(3)):
        return 1
    if all(sp.simplify(image[index] + vector[index]) == 0 for index in range(3)):
        return -1
    return "MIXED"


toroidal_x = [weight * component for component in [x, y, z]]
toroidal_y = [weight * component for component in cross([0, 1, 0], r_vector)]
poloidal_x = curl(toroidal_x)
quadrature = [
    sp.cos(omega * t) * toroidal_x[index]
    + sp.sin(omega * t) * toroidal_y[index]
    for index in range(3)
]
mixed = [
    toroidal_x[index] + mixing * poloidal_x[index] for index in range(3)
]


def helicity(vector):
    vector_curl = curl(vector)
    return integrate_all(
        sum(vector[index] * vector_curl[index] for index in range(3))
    )


print("HELICITY_TOROIDAL", helicity(toroidal_x))
print("HELICITY_TOROIDAL_QUADRATURE", helicity(quadrature))
print("HELICITY_TOROIDAL_POLOIDAL", helicity(mixed))
print("PARITY_TOROIDAL", inversion_parity(toroidal_x))
print("PARITY_POLOIDAL", inversion_parity(poloidal_x))

quadrature_time_derivative = [sp.diff(component, t) for component in quadrature]
angular_momentum = [
    integrate_all(component)
    for component in cross(quadrature, quadrature_time_derivative)
]
print("ANGULAR_MOMENTUM_TOROIDAL_QUADRATURE", angular_momentum)

swirl_axis = [0, 0, Omega]
toroidal_z = [weight * component for component in cross([0, 0, 1], r_vector)]
coriolis_image = cross(swirl_axis, toroidal_z)
coriolis_divergence = sp.simplify(divergence(coriolis_image))
print(
    "CORIOLIS_SCALAR_L0_PROJECTION",
    integrate_all(coriolis_divergence * weight),
)
print(
    "CORIOLIS_SCALAR_L2_PROJECTION",
    integrate_all(
        coriolis_divergence * (2 * z**2 - x**2 - y**2) * weight
    ),
)

# A vector fixed by all infinitesimal rotations must lie in the joint nullspace
# of the three rotation generators.
Lx, Ly, Lz = sp.symbols("Lx Ly Lz")
carrier = sp.Matrix([Lx, Ly, Lz])
generators = (
    sp.Matrix([[0, 0, 0], [0, 0, -1], [0, 1, 0]]),
    sp.Matrix([[0, 0, 1], [0, 0, 0], [-1, 0, 0]]),
    sp.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 0]]),
)
fixed_equations = []
for generator in generators:
    fixed_equations.extend(list(generator * carrier))
print(
    "O3_FIXED_VECTOR_SOLUTIONS",
    sp.solve(fixed_equations, (Lx, Ly, Lz), dict=True),
)
