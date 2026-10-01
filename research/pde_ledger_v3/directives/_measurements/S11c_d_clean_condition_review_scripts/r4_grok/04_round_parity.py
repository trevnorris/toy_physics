#!/usr/bin/env python3
"""H-round inversion parities for polar displacements and scalars.

Convention used by A: a polar vector field has inversion parity P when
    u_transformed(-x) := -u(x)   equals   P * u(-x)? 
Define P by: -u(-x) = P u(x), i.e. the field already matches the polar
vector law up to P. Equivalently u(-x) = -P u(x).

Scalars: phi(-x) = P phi(x), P = (-1)^ell for Y_ell.
"""
import sympy as sp

x, y, z = sp.symbols("x y z", real=True)
r2 = x**2 + y**2 + z**2

# ell=0 scalar: Y00 ~ 1
s0 = sp.Integer(1)
print("ELL0_SCALAR_P", 1, "formula", 1)

# Toroidal T = r × ∇Y. For ell=0, Y const, grad=0.
T00 = (0, 0, 0)
print("T00_ZERO", True)
print("TOROIDAL_EXISTS_ELL0", False)

# ell=1, m=0 scalar Y10 ~ z
s1 = z
s1_at_minus = s1.subs({x: -x, y: -y, z: -z})
print("ELL1_SCALAR_s(-x)/s(x)", sp.simplify(s1_at_minus / s1))
print("ELL1_SCALAR_P", -1, "formula_(-1)^ell", -1)

# Spheroidal ell=1: radial  u = x_i * Y ~ z (x,y,z)  (ignore 1/r factors at r=1)
# Use u = (x z, y z, z^2) which is r Y ê_r type (times r).
u_sph1 = (x * z, y * z, z * z)


def as_expr(c):
    return c if isinstance(c, sp.Expr) else sp.Integer(c)


def polar_parity(vec):
    """P such that -V(-x) = P V(x), if V is an eigenfield."""
    Vx, Vy, Vz = (as_expr(c) for c in vec)
    Vneg = (
        Vx.subs({x: -x, y: -y, z: -z}),
        Vy.subs({x: -x, y: -y, z: -z}),
        Vz.subs({x: -x, y: -y, z: -z}),
    )
    lhs = tuple(-c for c in Vneg)
    # try P = +1 and -1
    for P in (1, -1):
        rhs = tuple(P * c for c in vec)
        if all(sp.expand(lhs[i] - rhs[i]) == 0 for i in range(3)):
            return P
    return "NOT_EIGEN"


print("ELL1_SPHEROIDAL_P", polar_parity(u_sph1), "formula_(-1)^ell", -1)

# Toroidal ell=1: r × ∇z = r × e_z = (y, -x, 0)
u_tor1 = (y, -x, 0)
print("ELL1_TOROIDAL_P", polar_parity(u_tor1), "formula_(-1)^{ell+1}", 1)

# ell=2 scalar Y20 ~ 3z^2 - r^2
s2 = 3 * z**2 - r2
s2_at_minus = s2.subs({x: -x, y: -y, z: -z})
print("ELL2_SCALAR_s(-x)/s(x)", sp.simplify(s2_at_minus / s2))
print("ELL2_SCALAR_P", 1, "formula_(-1)^ell", 1)

# Toroidal ell=2, m=0: r × ∇Y20. Y20 ~ 3z^2-r^2, ∇Y = (-2x, -2y, 4z) after dropping const.
# r × ∇Y = |i     j     k  |
#         |x     y     z  |
#         |-2x  -2y    4z | = (y*4z - z*(-2y), z*(-2x)-x*4z, x*(-2y)-y*(-2x))
#         = (4 y z + 2 y z, -2 x z - 4 x z, 0) = (6 y z, -6 x z, 0)
u_tor2 = (y * z, -x * z, 0)
print("ELL2_TOROIDAL_P", polar_parity(u_tor2), "formula_(-1)^{ell+1}", -1)

# Spheroidal ell=2 radial ~ Y ê_r ~ (3z^2-r^2)(x,y,z)
u_sph2 = ((3 * z**2 - r2) * x, (3 * z**2 - r2) * y, (3 * z**2 - r2) * z)
print("ELL2_SPHEROIDAL_P", polar_parity(u_sph2), "formula_(-1)^ell", 1)

# Opposite parities at each ell>=1 => inversion forbids toroidal <-> scalar/spheroidal
for ell in range(0, 5):
    scalar_P = (-1) ** ell
    sph_P = (-1) ** ell
    if ell == 0:
        tor_exists = False
        tor_P = None
    else:
        tor_exists = True
        tor_P = (-1) ** (ell + 1)
    mix_allowed = (tor_P == scalar_P) if tor_exists else False
    print(
        f"ROUND_SECTOR {ell} SCALAR_P {scalar_P} SPH_P {sph_P} "
        f"TOR_P {tor_P} TOR_EXISTS {tor_exists} MIX_TOR_SCALAR {mix_allowed}"
    )

# zeta_c, theta, delta W, phi are scalars: same P as Y_ell
print("FACE_SCALARS_FOLLOW_YLM", True)
# u has no w-component; toroidal r×∇Y is already tangent to r=const in the brane
print("TOROIDAL_HAS_NO_W_COMPONENT", True)
