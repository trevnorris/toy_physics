#!/usr/bin/env python3
"""Small symbolic audit of the S11c-b reflection structures.

This independently reconstructs the delta-contraction family used by the
S11c-b SymPy engine and checks the z-reflection restriction.  It also checks
the exact tensor types in the anchoring, constraint, and tilted-face laws.
"""

from itertools import product
import sympy as sp


def perfect_matchings(items):
    if not items:
        return ((),)
    first = items[0]
    out = []
    for n in range(1, len(items)):
        rest = items[1:n] + items[n + 1 :]
        for matching in perfect_matchings(rest):
            out.append(((first, items[n]),) + matching)
    return tuple(out)


def component(data, indices):
    if not indices:
        return sp.sympify(data)
    if len(indices) == 1:
        return data[indices[0]]
    return data[indices[0]][indices[1]]


def delta_contractions(factors):
    slots = tuple(
        (factor, index)
        for factor, (_, rank, _) in enumerate(factors)
        for index in range(rank)
    )
    if len(slots) % 2:
        return ()
    expressions = []
    for matching in perfect_matchings(slots):
        total = 0
        for assignment in product(range(3), repeat=len(matching)):
            index_map = {
                slot: assignment[pair_index]
                for pair_index, pair in enumerate(matching)
                for slot in pair
            }
            term = 1
            for factor_index, (_, rank, data) in enumerate(factors):
                indices = tuple(index_map[(factor_index, slot)] for slot in range(rank))
                term *= component(data, indices)
            total += term
        expanded = sp.expand(total)
        if expanded not in expressions:
            expressions.append(expanded)
    return tuple(expressions)


def unique(expressions):
    out = []
    for expression in expressions:
        expression = sp.expand(expression)
        if expression != 0 and expression not in out:
            out.append(expression)
    return tuple(out)


P = sp.diag(1, 1, -1)
P4 = sp.diag(1, 1, -1, 1)

u = sp.Matrix(sp.symbols("u_x u_y u_z"))
g = sp.Matrix(sp.symbols("g_x g_y g_z"))
q = sp.Matrix(sp.symbols("q_x q_y q_z"))
r = sp.Matrix(sp.symbols("r_x r_y r_z"))
G_symbols = sp.symbols("G_xx G_xy G_xz G_yx G_yy G_yz G_zx G_zy G_zz")
Gm = sp.Matrix(3, 3, G_symbols)
G = tuple(tuple(Gm[a, i] for i in range(3)) for a in range(3))
theta, e = sp.symbols("theta e")

print("REFLECTION_MATRIX=", P)
print("DOT_VECTOR_RESIDUAL=", sp.expand((P * u).dot(P * g) - u.dot(g)))
print("TRACE_GRADIENT_RESIDUAL=", sp.expand(sp.trace(P * Gm * P) - sp.trace(Gm)))
constraint = theta + e + sp.trace(Gm) + u.dot(g)
constraint_reflected = theta + e + sp.trace(P * Gm * P) + (P * u).dot(P * g)
print("CONSTRAINT_FOLD_SCALAR_RESIDUAL=", sp.expand(constraint_reflected - constraint))
print("MATERIAL_ANCHOR_MINUS_U_DOT_G_RESIDUAL=", sp.expand(-(P * u).dot(P * g) + u.dot(g)))

slope = sp.Matrix(sp.symbols("h_x h_y h_z"))
normal = sp.Matrix([-slope[0], -slope[1], -slope[2], sp.Symbol("n_w")])
vface = sp.Matrix(sp.symbols("v_x v_y v_z v_w"))
normal_reflected = sp.Matrix([-(P * slope)[0], -(P * slope)[1], -(P * slope)[2], normal[3]])
print("TILTED_NORMAL_COVARIANCE_RESIDUAL=", sp.simplify(normal_reflected - P4 * normal))
print(
    "FACE_NORMAL_VELOCITY_SCALAR_RESIDUAL=",
    sp.expand(normal_reflected.dot(P4 * vface) - normal.dot(vface)),
)

z_restriction = {
    g[2]: 0,
    q[2]: 0,
    r[2]: 0,
    slope[2]: 0,
    **{Gm[a, 2]: 0 for a in range(3)},
}
print(
    "Z_INVARIANT_FACE_SCALAR_D_DVZ=",
    sp.diff(normal.dot(vface).subs({slope[2]: 0}), vface[2]),
)
print(
    "Z_INVARIANT_ANCHOR_D_DUZ=",
    sp.diff(u.dot(g).subs({g[2]: 0}), u[2]),
)

# Reconstruct every raw Kronecker candidate, a superset of the quotient-selected
# 10 uniform + 15 per first-jet source used by the actual engine.
data_uniform = (
    ("GRAD_U", 2, G),
    ("THETA", 0, theta),
    ("GRAD_THETA", 1, tuple(q)),
    ("E", 0, e),
    ("GRAD_E", 1, tuple(r)),
)
uniform_raw = []
for i, left in enumerate(data_uniform):
    for right in data_uniform[i:]:
        uniform_raw.extend(delta_contractions((left, right)))
uniform_raw = unique(uniform_raw)

data_new = (
    ("U", 1, tuple(u)),
    ("GRAD_U", 2, G),
    ("THETA", 0, theta),
    ("GRAD_THETA", 1, tuple(q)),
    ("E", 0, e),
    ("GRAD_E", 1, tuple(r)),
)
new_raw = []
spurion = ("BACKGROUND_FIRST_JET", 1, tuple(g))
for i, left in enumerate(data_new):
    for right in data_new[i:]:
        new_raw.extend(delta_contractions((spurion, left, right)))
new_raw = unique(new_raw)

odd = (u[2], Gm[2, 0], Gm[2, 1])
even = (
    theta,
    e,
    q[0],
    q[1],
    r[0],
    r[1],
    Gm[0, 0],
    Gm[0, 1],
    Gm[1, 0],
    Gm[1, 1],
)


def mixed_failures(expressions):
    failures = []
    for index, expression in enumerate(expressions):
        reduced = sp.expand(expression.subs(z_restriction))
        for odd_var in odd:
            for even_var in even:
                residual = sp.expand(sp.diff(reduced, odd_var, even_var))
                if residual != 0:
                    failures.append((index, odd_var, even_var, residual))
    return failures


uniform_failures = mixed_failures(uniform_raw)
new_failures = mixed_failures(new_raw)
print("UNIFORM_RAW_KRONECKER_CANDIDATE_COUNT=", len(uniform_raw))
print("UNIFORM_ODD_EVEN_MIXED_FAILURE_COUNT=", len(uniform_failures))
print("FIRST_JET_RAW_KRONECKER_CANDIDATE_COUNT=", len(new_raw))
print("FIRST_JET_ODD_EVEN_MIXED_FAILURE_COUNT=", len(new_failures))

# A single Levi-Civita contraction is parity odd.  It can make the protected
# block live, but not every syntactically valid one survives class P.
epsilon_active = sp.expand(g.dot(u.cross(q)))
epsilon_reflected = sp.expand((P * g).dot((P * u).cross(P * q)))
epsilon_active_z = sp.expand(epsilon_active.subs(z_restriction))
epsilon_inert = sp.expand(g.dot(q.cross(r)))
epsilon_inert_z = sp.expand(epsilon_inert.subs(z_restriction))
print("LEVI_CIVITA_PARITY_SUM_EPRIME_PLUS_E=", sp.expand(epsilon_reflected + epsilon_active))
print("K2_ACTIVE_Z_RESTRICTION=", epsilon_active_z)
print("K2_ACTIVE_MIXED_DUZ_DQX=", sp.diff(epsilon_active_z, u[2], q[0]))
print("K2_ACTIVE_MIXED_DUZ_DQY=", sp.diff(epsilon_active_z, u[2], q[1]))
print("K2_INERT_Z_RESTRICTION=", epsilon_inert_z)

# Two terms both satisfying B's prose description for K1 have different
# behavior in the target block.
beta = sp.symbols("beta", real=True)
axis = sp.Matrix([0, sp.sin(beta), sp.cos(beta)])
k1_inert = sp.expand(axis.dot(q) * sp.trace(Gm))
k1_active = sp.expand(axis.dot(q) * (axis.T * Gm * axis)[0])
k1_inert_z = sp.expand(k1_inert.subs(z_restriction))
k1_active_z = sp.expand(k1_active.subs(z_restriction))
print("K1_INERT_MIXED_DGZY_DQY=", sp.diff(k1_inert_z, Gm[2, 1], q[1]))
print("K1_ACTIVE_MIXED_DGZY_DQY=", sp.factor(sp.diff(k1_active_z, Gm[2, 1], q[1])))

# Inversion parity of scalar, spheroidal polar-vector, and toroidal
# polar-displacement harmonics.  The output is independent of m.
print("ROUND_INVERSION_PARITY_TABLE_l_scalar_spheroidal_toroidal=")
for ell in range(5):
    scalar = (-1) ** ell
    spheroidal = (-1) ** ell
    toroidal = (-1) ** (ell + 1)
    print((ell, scalar, spheroidal, toroidal))

odd_amplitude = sp.symbols("a_odd")
print("ODD_LINEAR_REFLECTION=", -odd_amplitude)
print("ODD_SQUARED_REFLECTION_RESIDUAL=", sp.expand((-odd_amplitude) ** 2 - odd_amplitude**2))
