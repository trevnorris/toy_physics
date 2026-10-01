#!/usr/bin/env python3
"""Independent reconstruction of z-reflection parity and Kronecker invariants
on the z-independent planar sector W=W(y).

Coordinates: 0=x (interface, k live), 1=y (background), 2=z (ignorable).
Does not import the S11c-b engine.
"""
from itertools import combinations_with_replacement, product
import sympy as sp

# --- z-reflection on z-independent fields ---------------------------------
x, y, z = sp.symbols("x y z", real=True)
# Independent jet symbols for a z-independent configuration
u1, u2, u3 = sp.symbols("u_x u_y u_z", real=True)
ux_x, ux_y, uy_x, uy_y, uz_x, uz_y = sp.symbols(
    "u_x_dx u_x_dy u_y_dx u_y_dy u_z_dx u_z_dy", real=True
)
th, th_x, th_y = sp.symbols("theta theta_dx theta_dy", real=True)
ew, ew_x, ew_y = sp.symbols("e_W e_W_dx e_W_dy", real=True)
zp, zm = sp.symbols("zeta_plus zeta_minus", real=True)
phi, dp = sp.symbols("phi delta_p", real=True)
gy = sp.symbols("g_y", real=True)
# background gradient along y only
g = (sp.Integer(0), gy, sp.Integer(0))

# Parity under z -> -z of a polar vector / true scalar, z-independent.
# u_z flips; everything else in this list does not.
parity = {
    "u_x": +1,
    "u_y": +1,
    "u_z": -1,
    "theta": +1,
    "e_W": +1,
    "zeta_plus": +1,
    "zeta_minus": +1,
    "phi": +1,
    "delta_p": +1,
    "g_y": +1,
    "dx": +1,
    "dy": +1,
    "dz_absent": +1,
}

print("=== PARITY TABLE (z -> -z, z-independent polar fields) ===")
for name, p in parity.items():
    print(f"PARITY[{name}] = {p}")

# --- Kronecker contractions, same algorithm as S11c-b delta_contractions ---
def perfect_matchings(items):
    items = tuple(items)
    if not items:
        yield ()
        return
    a = items[0]
    for i, b in enumerate(items[1:], 1):
        rest = items[1:i] + items[i + 1 :]
        for m in perfect_matchings(rest):
            yield ((a, b),) + m


def tensor_component(data, indices):
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
        total = sp.Integer(0)
        for assignment in product(range(3), repeat=len(matching)):
            index_map = {
                slot: assignment[pair_index]
                for pair_index, pair in enumerate(matching)
                for slot in pair
            }
            term = sp.Integer(1)
            for factor_index, (_, rank, data) in enumerate(factors):
                indices = tuple(index_map[(factor_index, slot)] for slot in range(rank))
                term *= tensor_component(data, indices)
            total += term
        expressions.append(sp.expand(total))
    # unique
    out = []
    seen = set()
    for e in expressions:
        key = sp.srepr(e)
        if key not in seen:
            seen.add(key)
            out.append(e)
    return tuple(out)


# Abstract jets, z-independent: G[a][i] = 0 if i==2 (no z derivative)
G = [[sp.zeros(1)[0] for i in range(3)] for a in range(3)]
G[0][0], G[0][1] = ux_x, ux_y
G[1][0], G[1][1] = uy_x, uy_y
G[2][0], G[2][1] = uz_x, uz_y
# G[*][2] remain 0
bu = (u1, u2, u3)
bq = (th_x, th_y, sp.Integer(0))
br = (ew_x, ew_y, sp.Integer(0))
btheta, be = th, ew

uniform_data = (
    ("GRAD_U", 2, G),
    ("THETA", 0, btheta),
    ("GRAD_THETA", 1, bq),
    ("E_LOCAL", 0, be),
    ("GRAD_E_LOCAL", 1, br),
)
raw_uniform = []
for li, left in enumerate(uniform_data):
    for right in uniform_data[li:]:
        raw_uniform.extend(delta_contractions((left, right)))

antisym = sp.expand(
    sum((G[i][j] - G[j][i]) ** 2 for i in range(3) for j in range(i + 1, 3))
)
symtf = sp.expand(
    sp.Rational(1, 2) * sum((G[i][j] + G[j][i]) ** 2 for i in range(3) for j in range(3))
    - sp.Rational(2, 3) * (G[0][0] + G[1][1] + G[2][2]) ** 2
)

def unique_expr(seq):
    out, seen = [], set()
    for e in seq:
        e = sp.expand(e)
        k = sp.srepr(e)
        if k not in seen and e != 0:
            seen.add(k)
            out.append(e)
    return out

uniform_invariants = unique_expr([antisym, symtf, *raw_uniform])
print("\n=== UNIFORM INVARIANTS restricted to z-independent jets ===")
print(f"COUNT_UNIFORM_NONEZERO_ON_P = {len(uniform_invariants)}")

def u_z_even_cross(expr):
    """True if expr has a monomial mixing an odd (u_z jet) with an even non-u_z field."""
    odd_syms = {u3, uz_x, uz_y}
    even_fields = {u1, u2, ux_x, ux_y, uy_x, uy_y, th, th_x, th_y, ew, ew_x, ew_y}
    expr = sp.expand(expr)
    if expr == 0:
        return False, sp.Integer(0)
    mixed = sp.Integer(0)
    for term in sp.Add.make_args(expr):
        odd_deg = sum(term.as_poly(list(odd_syms)).degree(s) if term.has(s) else 0 for s in odd_syms)
        # simpler: total degree in odd_syms
        t = term
        deg_odd = 0
        for s in odd_syms:
            d = t.as_expr().count_ops()  # dummy
        # use exponents
        factors = {s: 0 for s in list(odd_syms) + list(even_fields)}
        tt = sp.Integer(1) * term
        for s in list(odd_syms) + list(even_fields):
            if tt.has(s):
                factors[s] = sp.degree(sp.Poly(sp.expand(tt), s), s) if tt.has(s) else 0
        # Poly degree on each
        odd_total = 0
        even_total = 0
        for s in odd_syms:
            if term.has(s):
                odd_total += sp.Poly(sp.expand(term), s).degree()
        for s in even_fields:
            if term.has(s):
                even_total += sp.Poly(sp.expand(term), s).degree()
        if odd_total % 2 == 1 and even_total >= 1:
            mixed += term
    return mixed != 0, sp.expand(mixed)

print("\n=== UNIFORM INVARIANT even/odd MIXING (u_z with even fields) ===")
n_mix_u = 0
for i, inv in enumerate(uniform_invariants, 1):
    mixes, mixed = u_z_even_cross(inv)
    # also report whether invariant contains u_z at all
    has_uz = inv.has(u3, uz_x, uz_y)
    has_even = any(inv.has(s) for s in (u1, u2, ux_x, ux_y, uy_x, uy_y, th, th_x, th_y, ew, ew_x, ew_y))
    print(f"U{i:02d} has_uz={has_uz} mixes_odd_even={mixes}")
    print(f"    expr = {inv}")
    if mixes:
        n_mix_u += 1
        print(f"    MIXED = {mixed}")
print(f"UNIFORM_ODD_EVEN_MIX_COUNT = {n_mix_u}")

# --- spurion first-jet invariants ---
new_data = (
    ("U", 1, bu),
    ("GRAD_U", 2, G),
    ("THETA", 0, btheta),
    ("GRAD_THETA", 1, bq),
    ("E_LOCAL", 0, be),
    ("GRAD_E_LOCAL", 1, br),
)
spurion = ("BACKGROUND_FIRST_JET", 1, g)
raw_new = []
for li, left in enumerate(new_data):
    for right in new_data[li:]:
        raw_new.extend(delta_contractions((spurion, left, right)))
new_invariants = unique_expr(raw_new)
print("\n=== FIRST-JET SPURION INVARIANTS on g=(0,g_y,0), z-independent ===")
print(f"COUNT_SPURION_NONEZERO_ON_P = {len(new_invariants)}")
n_mix_s = 0
for i, inv in enumerate(new_invariants, 1):
    mixes, mixed = u_z_even_cross(inv)
    has_uz = inv.has(u3, uz_x, uz_y)
    print(f"S{i:02d} has_uz={has_uz} mixes_odd_even={mixes}")
    print(f"    expr = {inv}")
    if mixes:
        n_mix_s += 1
        print(f"    MIXED = {mixed}")
print(f"SPURION_ODD_EVEN_MIX_COUNT = {n_mix_s}")

# Explicit dangerous candidates
print("\n=== EXPLICIT DANGEROUS DENSITIES ===")
u_cross = u3 * (g[1] * th_x)  # u_z g_y d_x theta
print(f"u_z * g_y * d_x theta = {u_cross}")
print("This is z-odd (u_z odd, rest even). Kronecker family mix count above is the test.")

# curl^2 cross check
curl = (
    G[2][1] - G[1][2],  # d_y u_z - d_z u_y = uz_y
    G[0][2] - G[2][0],  # d_z u_x - d_x u_z = -uz_x
    G[1][0] - G[0][1],  # d_x u_y - d_y u_x
)
curl2 = sp.expand(sum(c**2 for c in curl))
print(f"\nCURL_SQUARED_ON_P = {curl2}")
mixes, mixed = u_z_even_cross(curl2)
print(f"CURL_SQUARED mixes_odd_even={mixes} mixed={mixed}")
