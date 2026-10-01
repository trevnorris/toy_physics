#!/usr/bin/env python3
"""H-round inversion parities for polar vector spherical harmonics, ell=1 and ell=2.

Polar-vector inversion image: u'( -x ) = - u(x) composed with the point map,
i.e. image_i(x) = - u_i(-x). Eigenvalue +1 means image = +u (A's toroidal
(-1)^{ell+1} for ell=1 is +1).
"""
import sympy as sp

x, y, z = sp.symbols("x y z", real=True)
r2 = x**2 + y**2 + z**2
weight = sp.exp(-r2)
r = sp.Matrix([x, y, z])


def curl(v):
    return sp.Matrix(
        [
            sp.diff(v[2], y) - sp.diff(v[1], z),
            sp.diff(v[0], z) - sp.diff(v[2], x),
            sp.diff(v[1], x) - sp.diff(v[0], y),
        ]
    )


def polar_image(v):
    flipped = v.subs({x: -x, y: -y, z: -z}, simultaneous=True)
    return sp.simplify(-flipped)


def parity_eigenvalue(v):
    img = polar_image(v)
    if all(sp.simplify(img[i] - v[i]) == 0 for i in range(3)):
        return 1
    if all(sp.simplify(img[i] + v[i]) == 0 for i in range(3)):
        return -1
    return "MIXED"


def helicity(v):
    c = curl(v)
    integrand = sp.expand(sum(v[i] * c[i] for i in range(3)))
    return sp.simplify(
        sp.integrate(integrand, (x, -sp.oo, sp.oo), (y, -sp.oo, sp.oo), (z, -sp.oo, sp.oo))
    )


# ell=1 toroidal: r × e_x * weight  (Y_1m dipole)
t1 = weight * r.cross(sp.Matrix([1, 0, 0]))
p1 = curl(t1)  # poloidal / spheroidal image of the same harmonic family

print("ELL1_TOROIDAL_PARITY", parity_eigenvalue(t1))
print("ELL1_POLOIDAL_PARITY", parity_eigenvalue(p1))
print("ELL1_A_CLAIMS_TOROIDAL", (-1) ** (1 + 1))
print("ELL1_A_CLAIMS_SPHEROIDAL", (-1) ** 1)
print("ELL1_TOROIDAL_HELICITY", helicity(t1))
print("ELL1_MIXED_HELICITY", helicity(t1 + p1))

# ell=2 toroidal: r × ∇(quadratic harmonic Q = x y), times weight.
Q2 = x * y
gradQ = sp.Matrix([sp.diff(Q2, x), sp.diff(Q2, y), sp.diff(Q2, z)])
t2 = weight * r.cross(gradQ)
p2 = curl(t2)
print("ELL2_TOROIDAL_PARITY", parity_eigenvalue(t2))
print("ELL2_POLOIDAL_PARITY", parity_eigenvalue(p2))
print("ELL2_A_CLAIMS_TOROIDAL", (-1) ** (2 + 1))
print("ELL2_A_CLAIMS_SPHEROIDAL", (-1) ** 2)

# Scalar Y_ellm parity is (-1)^ell without the polar minus.
def scalar_parity(s):
    img = sp.simplify(s.subs({x: -x, y: -y, z: -z}, simultaneous=True))
    if sp.simplify(img - s) == 0:
        return 1
    if sp.simplify(img + s) == 0:
        return -1
    return "MIXED"


print("SCALAR_XY_PARITY_ELL2", scalar_parity(weight * Q2))
print("SCALAR_Z_PARITY_ELL1", scalar_parity(weight * z))

# A net in-plane vector is not O(3) fixed.
Lx, Ly, Lz = sp.symbols("Lx Ly Lz")
carrier = sp.Matrix([Lx, Ly, Lz])
gens = (
    sp.Matrix([[0, 0, 0], [0, 0, -1], [0, 1, 0]]),
    sp.Matrix([[0, 0, 1], [0, 0, 0], [-1, 0, 0]]),
    sp.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 0]]),
)
eqs = []
for g in gens:
    eqs.extend(list(g * carrier))
print("O3_FIXED_VECTOR", sp.solve(eqs, [Lx, Ly, Lz], dict=True))
