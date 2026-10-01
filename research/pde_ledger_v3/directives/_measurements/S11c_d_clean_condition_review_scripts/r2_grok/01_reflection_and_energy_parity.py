#!/usr/bin/env python3
"""Independent parity audit of the S11c-b field content on class P.

Background restriction R1: profile jet g = (g1, 0, 0).
Class P: all fields independent of direction 3 (index 2 in 0-based loops).
Reflection R: x3 -> -x3, polar-vector rule (Ru)_i(x) = R_ij u_j(R x).
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


# --- reflection action on z-independent fields ---
print("=== 1. Reflection representation on class P ===")
# For x3-independent fields, Ru = +u (even sector) requires:
#   u1, u2 even functions of x3 -> allowed as constants in x3
#   u3 odd function of x3 -> must vanish
# Ru = -u (odd sector) requires:
#   u1, u2 odd in x3 -> must vanish
#   u3 even in x3 -> allowed
# Scalars: even sector allowed; odd sector must vanish.
print("even_sector_fields = u1, u2, theta, e_W, zeta_plus, zeta_minus, phi, delta_p")
print("odd_sector_fields  = u3")
print("odd_scalar_on_P    = 0  (cannot be x3-independent and odd unless identically 0)")

# Polar vector: odd part is the 3-component.
R = sp.diag(1, 1, -1)
u = sp.Matrix(sp.symbols("u1 u2 u3"))
print("R_on_polar_u =", tuple(R * u))
print("odd_part_of_polar_u = u3")

n = sp.Matrix(sp.symbols("n1 n2 n3"))
print("R_on_polar_n =", tuple(R * n))

# --- Kronecker energy family, R1 + P ---
print("\n=== 2. Kronecker bilinear family on R1+P ===")
g1 = sp.symbols("g1")
g = (g1, 0, 0)
u1, u2, u3 = sp.symbols("u1 u2 u3")
th, ew = sp.symbols("theta e_W")
# z-independent first jets (index order: field component, derivative direction)
du = (
    tuple(sp.symbols(f"d{k+1}_u1") for k in range(2)) + (0,),
    tuple(sp.symbols(f"d{k+1}_u2") for k in range(2)) + (0,),
    tuple(sp.symbols(f"d{k+1}_u3") for k in range(2)) + (0,),
)
dth = tuple(sp.symbols(f"d{k+1}_th") for k in range(2)) + (0,)
dew = tuple(sp.symbols(f"d{k+1}_ew") for k in range(2)) + (0,)
bu = (u1, u2, u3)
bG = du
bq = dth
br = dew

odd_atoms = {u3, du[2][0], du[2][1]}
even_dyn = {
    u1, u2, th, ew,
    du[0][0], du[0][1], du[1][0], du[1][1],
    dth[0], dth[1], dew[0], dew[1],
}


def mixed_odd_even(expr):
    expr = sp.expand(expr)
    mixed = []
    for term in sp.Add.make_args(expr):
        odd_deg = 0
        even_deg = 0
        for s in odd_atoms:
            if term.has(s):
                odd_deg += sp.degree(sp.expand(term), s)
        for s in even_dyn:
            if term.has(s):
                even_deg += sp.degree(sp.expand(term), s)
        if odd_deg % 2 == 1 and even_deg >= 1:
            mixed.append(term)
    return sp.expand(sum(mixed)) if mixed else sp.Integer(0)


def unique(expressions):
    out = []
    for expression in expressions:
        expression = sp.expand(expression)
        if expression != 0 and expression not in out:
            out.append(expression)
    return tuple(out)


uniform_data = (
    ("GRAD_U", 2, bG),
    ("THETA", 0, th),
    ("GRAD_THETA", 1, bq),
    ("E_LOCAL", 0, ew),
    ("GRAD_E_LOCAL", 1, br),
)
raw_uniform = []
for i, left in enumerate(uniform_data):
    for right in uniform_data[i:]:
        raw_uniform.extend(delta_contractions((left, right)))
antisym = sp.expand(
    sum((bG[i][j] - bG[j][i]) ** 2 for i in range(3) for j in range(i + 1, 3))
)
raw_uniform = unique([antisym, *raw_uniform])

new_data = (
    ("U", 1, bu),
    ("GRAD_U", 2, bG),
    ("THETA", 0, th),
    ("GRAD_THETA", 1, bq),
    ("E_LOCAL", 0, ew),
    ("GRAD_E_LOCAL", 1, br),
)
spurion = ("BACKGROUND_FIRST_JET", 1, g)
raw_new = []
for i, left in enumerate(new_data):
    for right in new_data[i:]:
        raw_new.extend(delta_contractions((spurion, left, right)))
raw_new = unique(raw_new)

n_uniform_mixed = 0
n_new_mixed = 0
print(f"uniform_candidate_count = {len(raw_uniform)}")
for n, expr in enumerate(raw_uniform, 1):
    mix = mixed_odd_even(expr)
    flag = mix != 0
    n_uniform_mixed += int(flag)
    if flag:
        print(f"UNIFORM_MIXED[{n}] {expr} -> {mix}")
print(f"uniform_mixed_count = {n_uniform_mixed}")

print(f"new_spurion_candidate_count = {len(raw_new)}")
for n, expr in enumerate(raw_new, 1):
    mix = mixed_odd_even(expr)
    flag = mix != 0
    n_new_mixed += int(flag)
    if flag:
        print(f"NEW_MIXED[{n}] {expr} -> {mix}")
print(f"new_mixed_count = {n_new_mixed}")

# Explicit surviving u3-only pieces (odd-odd, allowed)
u3_only = []
for expr in list(raw_uniform) + list(raw_new):
    e2 = expr.subs({u1: 0, u2: 0, th: 0, ew: 0,
                    du[0][0]: 0, du[0][1]: 0, du[1][0]: 0, du[1][1]: 0,
                    dth[0]: 0, dth[1]: 0, dew[0]: 0, dew[1]: 0})
    e2 = sp.expand(e2)
    if e2 != 0:
        u3_only.append(e2)
print("u3_self_energy_on_P_nonzero =", unique(u3_only))

print("\n=== 3. Curl-squared cross terms on P ===")
curl = (
    du[2][1] - du[1][2],  # d2 u3 - d3 u2 = d2 u3
    du[0][2] - du[2][0],  # d3 u1 - d1 u3 = -d1 u3
    du[1][0] - du[0][1],  # d1 u2 - d2 u1
)
curl_sq = sp.expand(sum(c**2 for c in curl))
print("curl =", curl)
print("curl_sq =", curl_sq)
print("curl_sq_mixed =", mixed_odd_even(curl_sq))
