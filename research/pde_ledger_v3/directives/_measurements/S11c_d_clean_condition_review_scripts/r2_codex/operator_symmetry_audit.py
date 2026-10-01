#!/usr/bin/env python3
"""Independent small symbolic audit for the S11c-d clean-condition review.

This does not import or run either production engine.  It independently
reconstructs the Kronecker-contraction candidate families used at
S11c_b_brane_operator_sympy_audit.py:1313-1529, applies the proposed R1/P
restriction before testing reflection blocks, and models the exact tensor
structure of the supplied face/constraint laws.
"""

from itertools import product
import hashlib
from pathlib import Path
import sympy as sp

REPO = Path("/var/projects/toy_physics")
ENGINE = REPO / "research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py"


def perfect_matchings(items):
    if not items:
        return ((),)
    first = items[0]
    out = []
    for p in range(1, len(items)):
        rest = items[1:p] + items[p + 1 :]
        for matching in perfect_matchings(rest):
            out.append(((first, items[p]),) + matching)
    return tuple(out)


def component(data, indices):
    if not indices:
        return sp.sympify(data)
    if len(indices) == 1:
        return data[indices[0]]
    return data[indices[0]][indices[1]]


def delta_contractions(factors):
    slots = tuple((f, i) for f, (_, rank, _) in enumerate(factors) for i in range(rank))
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
            for f, (_, rank, data) in enumerate(factors):
                term *= component(data, tuple(index_map[(f, j)] for j in range(rank)))
            total += term
        expanded = sp.expand(total)
        if expanded != 0 and expanded not in expressions:
            expressions.append(expanded)
    return tuple(expressions)


def unique(expressions):
    out = []
    for expression in expressions:
        expression = sp.expand(expression)
        if expression != 0 and expression not in out:
            out.append(expression)
    return tuple(out)


def sidx(*indices):
    return tuple(sorted(indices))


u = sp.symbols("u1:4")
G = tuple(tuple(sp.symbols(f"G{a+1}{i+1}") for i in range(3)) for a in range(3))
theta = sp.Symbol("theta")
q = sp.symbols("q1:4")
e = sp.Symbol("e")
r = sp.symbols("r1:4")
g = sp.symbols("g1:4")


def second(prefix, component_index=None):
    tag = "" if component_index is None else str(component_index + 1)
    return {
        (i, j): sp.Symbol(f"{prefix}{tag}_{i+1}{j+1}")
        for i in range(3)
        for j in range(i, 3)
    }


fields = tuple(
    [(u[a], G[a], second("U2", a)) for a in range(3)]
    + [(theta, q, second("T2")), (e, r, second("E2"))]
)
g2 = second("BG2")


def signatures(candidates, live_background=False):
    derivative_maps = {i: {} for i in range(3)}
    for field, first, second_jets in fields:
        for i in range(3):
            derivative_maps[i][field] = first[i]
            for j in range(3):
                derivative_maps[i][first[j]] = second_jets[sidx(i, j)]
    if live_background:
        for i in range(3):
            for j in range(3):
                derivative_maps[i][g[j]] = g2[sidx(i, j)]

    def dx(expr, direction):
        return sp.expand(sum(sp.diff(expr, atom) * deriv for atom, deriv in derivative_maps[direction].items()))

    out = []
    for candidate in candidates:
        row = []
        for field, first, _ in fields:
            row.append(sp.expand(sp.diff(candidate, field) - sum(dx(sp.diff(candidate, first[i]), i) for i in range(3))))
        out.append(tuple(row))
    return tuple(out)


def independent_indices(candidates, rows):
    variables = sorted(
        {atom for row in rows for expr in row for atom in expr.free_symbols},
        key=sp.default_sort_key,
    )
    monomials = sorted(
        {mon for row in rows for expr in row for mon in sp.Poly(expr, *variables).monoms()}
    )
    matrix = sp.zeros(len(rows[0]) * len(monomials), 0)
    rank = 0
    selected = []
    for index, row in enumerate(rows):
        column = []
        for expr in row:
            poly = sp.Poly(expr, *variables)
            column.extend(poly.coeff_monomial(mon) for mon in monomials)
        trial = matrix.row_join(sp.Matrix(column))
        trial_rank = trial.rank()
        if trial_rank > rank:
            selected.append(index)
            matrix = trial
            rank = trial_rank
    return tuple(selected)


uniform_data = (
    ("GRAD_U", 2, G),
    ("THETA", 0, theta),
    ("GRAD_THETA", 1, q),
    ("E_LOCAL", 0, e),
    ("GRAD_E_LOCAL", 1, r),
)
uniform_raw = []
for left_index, left in enumerate(uniform_data):
    for right in uniform_data[left_index:]:
        uniform_raw.extend(delta_contractions((left, right)))
antisym = sp.expand(sum((G[i][j] - G[j][i]) ** 2 for i in range(3) for j in range(i + 1, 3)))
symtf = sp.expand(
    sp.Rational(1, 2) * sum((G[i][j] + G[j][i]) ** 2 for i in range(3) for j in range(3))
    - sp.Rational(2, 3) * sum(G[i][i] for i in range(3)) ** 2
)
uniform_candidates = unique([antisym, symtf, *uniform_raw])
uniform_selected = independent_indices(uniform_candidates, signatures(uniform_candidates))

new_data = (
    ("U", 1, u),
    ("GRAD_U", 2, G),
    ("THETA", 0, theta),
    ("GRAD_THETA", 1, q),
    ("E_LOCAL", 0, e),
    ("GRAD_E_LOCAL", 1, r),
)
new_raw = []
spurion = ("BACKGROUND_FIRST_JET", 1, g)
for left_index, left in enumerate(new_data):
    for right in new_data[left_index:]:
        new_raw.extend(delta_contractions((spurion, left, right)))
new_candidates = unique(new_raw)
new_selected = independent_indices(new_candidates, signatures(new_candidates, live_background=True))

# R1: only background direction 1 survives.  P: no direction-3 derivative.
r1p = {g[1]: 0, g[2]: 0, q[2]: 0, r[2]: 0}
r1p.update({G[a][2]: 0 for a in range(3)})
odd = (u[2], G[2][0], G[2][1])
even = (
    u[0], u[1],
    G[0][0], G[0][1], G[1][0], G[1][1],
    theta, q[0], q[1], e, r[0], r[1],
)
reflection = {atom: -atom for atom in odd}


def parity_and_block_failures(expressions):
    parity_failures = []
    block_failures = []
    for index, original in enumerate(expressions):
        expr = sp.expand(original.subs(r1p, simultaneous=True))
        parity_residual = sp.expand(expr.xreplace(reflection) - expr)
        if parity_residual != 0:
            parity_failures.append((index, parity_residual))
        for odd_atom in odd:
            for even_atom in even:
                mixed = sp.expand(sp.diff(expr, odd_atom, even_atom))
                if mixed != 0:
                    block_failures.append((index, odd_atom, even_atom, mixed))
    return parity_failures, block_failures


uniform_exprs = tuple(uniform_candidates[i] for i in uniform_selected)
new_exprs = tuple(new_candidates[i] for i in new_selected)
uniform_parity, uniform_blocks = parity_and_block_failures(uniform_exprs)
new_parity, new_blocks = parity_and_block_failures(new_exprs)


def restricted_unique_nonzero(expressions):
    restricted = [sp.expand(expr.subs(r1p, simultaneous=True)) for expr in expressions]
    return tuple(dict.fromkeys(expr for expr in restricted if expr != 0))

print("ENGINE_SHA256", hashlib.sha256(ENGINE.read_bytes()).hexdigest())
print("UNIFORM_SELECTED_COUNT", len(uniform_selected))
print("FIRST_JET_SELECTED_COUNT_PER_SOURCE", len(new_selected))
print("TOTAL_ENGINE_FAMILY_COUNT", len(uniform_selected) + 2 * len(new_selected))
print("R1P_DISTINCT_NONZERO_UNIFORM_CARRIERS", len(restricted_unique_nonzero(uniform_exprs)))
print("R1P_DISTINCT_NONZERO_FIRST_JET_CARRIERS_PER_SOURCE", len(restricted_unique_nonzero(new_exprs)))
print("R1P_UNIFORM_PARITY_FAILURES", uniform_parity)
print("R1P_UNIFORM_ODD_EVEN_HESSIAN_FAILURES", uniform_blocks)
print("R1P_FIRST_JET_PARITY_FAILURES", new_parity)
print("R1P_FIRST_JET_ODD_EVEN_HESSIAN_FAILURES", new_blocks)

# Exact supplied-law tensor structure on an R1/P graph.
s, h1, h2 = sp.symbols("s h1 h2")
ut1, ut2, ut3, ht = sp.symbols("ut1 ut2 ut3 ht")
vb1, vb2, vb3, vbw = sp.symbols("vb1 vb2 vb3 vbw")
mu, dp, LA, LV, LX, rhom = sp.symbols("mu dp LA LV LX rhom", nonzero=True)
n = sp.Matrix([-h1, -h2, 0, s])
vface = sp.Matrix([ut1, ut2, ut3, ht])
vbulk = sp.Matrix([vb1, vb2, vb3, vbw])
V = sp.expand(n.dot(vface))
normal_bulk = sp.expand(n.dot(vbulk))
affinity = mu - dp / rhom
J = sp.expand(LA * affinity + LV * V)
traction = sp.expand(-(dp + LX * affinity)) * n

print("FACE_NORMAL_R1P", tuple(n))
print("FACE_VELOCITY_DEPENDENCE_D_UT3", sp.diff(V, ut3))
print("BULK_NORMAL_DEPENDENCE_D_VB3", sp.diff(normal_bulk, vb3))
print("FLUX_DEPENDENCE_D_UT3", sp.diff(J, ut3))
print("TRACTION_COMPONENT_3", traction[2])

# Both anchoring pullbacks and pin-B constraint see only the direction-1 jet.
u1, u2, u3, Q1, rho1 = sp.symbols("u1 u2 u3 Q1 rho1")
dvt, dve, dvu1, dvu2, dvu3, dvG11, dvG22 = sp.symbols(
    "dvt dve dvu1 dvu2 dvu3 dvG11 dvG22"
)
lab_pullback = u1 * Q1
material_pullback = u1 * rho1
virtual_constraint = dvt + dve + dvG11 + dvG22 + rho1 * dvu1
print("LAB_ANCHOR_D_U3", sp.diff(lab_pullback, u3))
print("MATERIAL_ANCHOR_D_U3", sp.diff(material_pullback, u3))
print("PIN_B_CONSTRAINT_D_DVU3", sp.diff(virtual_constraint, dvu3))

# R1 plus one mirror does not by itself imply rotations about direction 1.
# A constant background/support vector along direction 2 is independent of
# directions 2/3 and is fixed by the 3-reflection, but not by a generic O(2)
# rotation in the 2-3 plane.
alpha, b = sp.symbols("alpha b", real=True)
R_about_1 = sp.Matrix([
    [1, 0, 0],
    [0, sp.cos(alpha), -sp.sin(alpha)],
    [0, sp.sin(alpha), sp.cos(alpha)],
])
mirror_3 = sp.diag(1, 1, -1)
background_vector = sp.Matrix([0, b, 0])
print("R1_VECTOR_MIRROR3_RESIDUAL", tuple(mirror_3 * background_vector - background_vector))
print("R1_VECTOR_ROTATION_ABOUT1_RESIDUAL", tuple(sp.simplify(x) for x in (R_about_1 * background_vector - background_vector)))

# Pinned controls, already reduced to R1/P.
a1, a2, beta, g1 = sp.symbols("a_K1 a_K2 beta g1")
u11, u12, u31, u32, th1, th2 = sp.symbols("u11 u12 u31 u32 th1 th2")
K1 = sp.expand(a1 * (sp.sin(beta) * (u11 * th1 + u12 * th2) + sp.cos(beta) * (u31 * th1 + u32 * th2)))
K2 = sp.expand(a2 * theta * g1 * u32)
K1_reflected = sp.expand(K1.subs({u31: -u31, u32: -u32}, simultaneous=True))
K2_reflected = sp.expand(K2.subs({u32: -u32}, simultaneous=True))
print("K1_REFLECTION_RESIDUAL", sp.factor(K1_reflected - K1))
print("K1_ODD_EVEN_MIXED_DERIVATIVES", sp.diff(K1, u31, th1), sp.diff(K1, u32, th2))
print("K2_R1P", K2)
print("K2_REFLECTION_RESIDUAL", sp.factor(K2_reflected - K2))
print("K2_ODD_EVEN_MIXED_DERIVATIVE", sp.diff(K2, u32, theta))

# Round vector-harmonic parity under the polar-vector parity operator P V(x)=-V(-x).
ell = sp.Symbol("ell", integer=True)
scalar_eval = (-1) ** ell
radial_or_poloidal_eval = (-1) ** (ell + 1)
toroidal_eval = (-1) ** ell
print("ROUND_SCALAR_PARITY", scalar_eval)
print("ROUND_SPHEROIDAL_POLAR_VECTOR_PARITY", sp.simplify(-radial_or_poloidal_eval))
print("ROUND_TOROIDAL_POLAR_VECTOR_PARITY", sp.simplify(-toroidal_eval))
print("ROUND_SAME_ELL_PARITY_RATIO_TOROIDAL_OVER_SCALAR", sp.simplify((-toroidal_eval) / scalar_eval))

# Linear selection does not imply nonlinear isolation: an odd-square/even-scalar term is allowed.
odd_amp, scalar_amp, lam = sp.symbols("odd_amp scalar_amp lambda")
nonlinear = lam * scalar_amp * odd_amp**2
print("NONLINEAR_TERM_REFLECTION_RESIDUAL", sp.expand(nonlinear.subs(odd_amp, -odd_amp) - nonlinear))
print("NONLINEAR_SCALAR_SOURCE", sp.diff(nonlinear, scalar_amp))
print("LINEAR_ODD_SCALAR_HESSIAN_AT_BACKGROUND", sp.diff(nonlinear, odd_amp, scalar_amp).subs(odd_amp, 0))

# A real physical pair of degenerate standing coordinates can carry angular momentum.
A, omega, t = sp.symbols("A omega t", real=True)
qx = A * sp.cos(omega * t)
qy_linear = sp.Integer(0)
qy_circular = A * sp.sin(omega * t)
L_linear = sp.simplify(qx * sp.diff(qy_linear, t) - qy_linear * sp.diff(qx, t))
L_circular = sp.simplify(qx * sp.diff(qy_circular, t) - qy_circular * sp.diff(qx, t))
print("REAL_LINEAR_STANDING_PAIR_ANGULAR_MOMENTUM", L_linear)
print("REAL_QUADRATURE_STANDING_PAIR_ANGULAR_MOMENTUM", L_circular)
