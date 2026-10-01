#!/usr/bin/env python3
"""Independent Kronecker energy-parity and K1/K2/G3 mixed-block audit.

Rebuilds the S11c-b first-jet enumerator (Kronecker contractions of
{u, grad u, theta, grad theta, e, grad e} with a spurion g) and tests
x3-reflection parity on class P (perturbation d3 = 0) with R1 (g2=g3=0).
Does not import the ledger engines.
"""
from itertools import combinations, product
import sympy as sp

# Abstract jets: polar vector u, scalars theta/e, polar spurion g.
u = sp.Matrix(sp.symbols("u1 u2 u3"))
G = sp.Matrix(3, 3, lambda i, j: sp.symbols(f"G{i+1}{j+1}"))  # G_ij = d_j u_i
theta, e = sp.symbols("theta e")
q = sp.Matrix(sp.symbols("q1 q2 q3"))  # d theta
r = sp.Matrix(sp.symbols("r1 r2 r3"))  # d e
g = sp.Matrix(sp.symbols("g1 g2 g3"))


def perfect_matchings(slots):
    slots = tuple(slots)
    if not slots:
        yield ()
        return
    a = slots[0]
    for i, b in enumerate(slots[1:], start=1):
        rest = slots[1:i] + slots[i + 1 :]
        for matching in perfect_matchings(rest):
            yield ((a, b),) + matching


def delta_contractions(factors):
    # factors: list of (rank, data) with data indexed by 0..2
    slots = [(f, k) for f, (rank, _) in enumerate(factors) for k in range(rank)]
    if len(slots) % 2:
        return []
    out = []
    for matching in perfect_matchings(slots):
        total = 0
        for assignment in product(range(3), repeat=len(matching)):
            index_map = {}
            for pair_i, pair in enumerate(matching):
                for slot in pair:
                    index_map[slot] = assignment[pair_i]
            term = 1
            for f, (rank, data) in enumerate(factors):
                idx = tuple(index_map[(f, k)] for k in range(rank))
                term *= data[idx]
            total += term
        out.append(sp.expand(total))
    # unique
    uniq = []
    seen = set()
    for expr in out:
        key = sp.srepr(expr)
        if key not in seen:
            seen.add(key)
            uniq.append(expr)
    return uniq


def tensor_get(name, rank, data):
    def getter(idx):
        if rank == 0:
            return data
        if rank == 1:
            return data[idx[0]]
        return data[idx[0], idx[1]]

    return getter


# Build factor tables as callables via dicts
def as_data(rank, obj):
    if rank == 0:
        return lambda idx: obj
    if rank == 1:
        return lambda idx: obj[idx[0]]
    return lambda idx: obj[idx[0], idx[1]]


field_data = [
    ("U", 1, as_data(1, u)),
    ("GRAD_U", 2, as_data(2, G)),
    ("THETA", 0, as_data(0, theta)),
    ("GRAD_THETA", 1, as_data(1, q)),
    ("E", 0, as_data(0, e)),
    ("GRAD_E", 1, as_data(1, r)),
]


def enumerate_new(gvec):
    spurion = ("G", 1, as_data(1, gvec))
    raw = []
    for i, left in enumerate(field_data):
        for right in field_data[i:]:
            factors = [
                (1, spurion[2]),
                (left[1], left[2]),
                (right[1], right[2]),
            ]
            # rebuild using indexable arrays for delta_contractions
            # Use symbol tensors instead
            pass
    return raw


# Simpler explicit contraction using symbols packed as nested tuples
def pack_vec(vec):
    return tuple(vec[i] for i in range(3))


def pack_mat(mat):
    return tuple(tuple(mat[i, j] for j in range(3)) for i in range(3))


bu = pack_vec(u)
bG = pack_mat(G)
bq = pack_vec(q)
br = pack_vec(r)
bg = pack_vec(g)
btheta = theta
be = e


def tcomp(data, indices):
    if not indices:
        return data
    if len(indices) == 1:
        return data[indices[0]]
    return data[indices[0]][indices[1]]


def contractions(factors):
    """factors: list of (rank, data). Kronecker pairings of all indices."""
    slots = [(f, k) for f, (rank, _) in enumerate(factors) for k in range(rank)]
    if len(slots) % 2:
        return []
    exprs = []
    for matching in perfect_matchings(slots):
        total = 0
        for assignment in product(range(3), repeat=len(matching)):
            index_map = {}
            for pair_i, pair in enumerate(matching):
                for slot in pair:
                    index_map[slot] = assignment[pair_i]
            term = 1
            for f, (rank, data) in enumerate(factors):
                idx = tuple(index_map[(f, k)] for k in range(rank))
                term *= tcomp(data, idx)
            total += term
        exprs.append(sp.expand(total))
    uniq = []
    seen = set()
    for expr in exprs:
        key = sp.srepr(sp.expand(expr))
        if key not in seen:
            seen.add(key)
            uniq.append(sp.expand(expr))
    return uniq


def new_candidates(g_tuple):
    data = [
        ("U", 1, bu),
        ("GRAD_U", 2, bG),
        ("THETA", 0, btheta),
        ("GRAD_THETA", 1, bq),
        ("E", 0, be),
        ("GRAD_E", 1, br),
    ]
    spurion = (1, g_tuple)
    raw = []
    for i, left in enumerate(data):
        for right in data[i:]:
            factors = [spurion, (left[1], left[2]), (right[1], right[2])]
            raw.extend(contractions(factors))
    uniq = []
    seen = set()
    for expr in raw:
        key = sp.srepr(expr)
        if key not in seen:
            seen.add(key)
            uniq.append(expr)
    return uniq


def uniform_candidates():
    data = [
        ("GRAD_U", 2, bG),
        ("THETA", 0, btheta),
        ("GRAD_THETA", 1, bq),
        ("E", 0, be),
        ("GRAD_E", 1, br),
    ]
    raw = []
    for i, left in enumerate(data):
        for right in data[i:]:
            raw.extend(contractions([(left[1], left[2]), (right[1], right[2])]))
    antisym = sp.expand(sum((bG[i][j] - bG[j][i]) ** 2 for i in range(3) for j in range(i + 1, 3)))
    return [antisym] + raw


# Reflection x3 -> -x3 on polar tensors; class P: perturbation d3 = 0.
odd_syms = {
    u[2],
    G[2, 0],
    G[2, 1],  # d1 u3, d2 u3
    G[0, 2],
    G[1, 2],
    G[2, 2],  # d3 u_*  (set to 0 on P, listed for completeness)
    q[2],
    r[2],
    g[2],
}
# On class P, G_i3 = d3 u_i = 0 and q3=r3=0.
P_subs = {
    G[0, 2]: 0,
    G[1, 2]: 0,
    G[2, 2]: 0,
    q[2]: 0,
    r[2]: 0,
}
R1_subs = {g[1]: 0, g[2]: 0}


def is_odd_monomial(monomial):
    """A monomial is odd if it contains an odd number of odd symbols."""
    deg = 0
    if monomial == 0:
        return False
    for sym in odd_syms:
        deg += sp.degree(monomial, sym)
    return (deg % 2) == 1


def split_parity(expr):
    expr = sp.expand(expr)
    even = 0
    odd = 0
    if expr == 0:
        return sp.Integer(0), sp.Integer(0)
    for term in sp.Add.make_args(expr):
        if is_odd_monomial(term):
            odd += term
        else:
            even += term
    return sp.expand(even), sp.expand(odd)


def mixed_u3_even(expr):
    """Terms that are bilinear in {u3-family} and {even scalars/u1/u2-family}."""
    expr = sp.expand(expr)
    even, odd = split_parity(expr)
    return odd  # energy density odd under the reflection = mixed coupling


new_inv = new_candidates(bg)
print("NEW_INVARIANT_COUNT", len(new_inv))

mixed_r1_p = []
for n, inv in enumerate(new_inv, start=1):
    restricted = sp.expand(inv.subs(P_subs).subs(R1_subs))
    even, odd = split_parity(restricted)
    if odd != 0:
        mixed_r1_p.append((n, odd))

print("KRONECKER_SPURION_ODD_ENERGY_ON_R1_P_COUNT", len(mixed_r1_p))
for item in mixed_r1_p[:12]:
    print("KRONECKER_ODD", item[0], item[1])

uniform = uniform_candidates()
print("UNIFORM_RAW_COUNT", len(uniform))
mixed_u = []
for n, inv in enumerate(uniform, start=1):
    restricted = sp.expand(inv.subs(P_subs))
    even, odd = split_parity(restricted)
    if odd != 0:
        mixed_u.append((n, odd))
print("UNIFORM_ODD_ENERGY_ON_P_COUNT", len(mixed_u))
for item in mixed_u[:8]:
    print("UNIFORM_ODD", item[0], item[1])

# Curl-squared on P (the photon kinetic/potential piece)
curl_sq = sp.expand(sum((bG[i][j] - bG[j][i]) ** 2 for i in range(3) for j in range(i + 1, 3)))
curl_p = sp.expand(curl_sq.subs(P_subs))
print("CURL_SQUARED_ON_P", curl_p)
print("CURL_SQUARED_ODD_PART", split_parity(curl_p)[1])

# K1: a_K1 ehat_i (d_k u_i)(d_k theta)
a_K1, beta = sp.symbols("a_K1 beta", real=True)
ehat = sp.Matrix([sp.sin(beta), 0, sp.cos(beta)])
K1 = 0
for i in range(3):
    for k in range(3):
        K1 += a_K1 * ehat[i] * G[i, k] * q[k]
K1 = sp.expand(K1.subs(P_subs))
print("K1_ON_P", K1)
print("K1_ODD_ENERGY", split_parity(K1)[1])

# If R1 were wrongly applied to ehat, ehat_3 -> 0
K1_killed = sp.expand(K1.subs(sp.cos(beta), 0))
print("K1_ODD_IF_E3_ZEROED", split_parity(K1_killed)[1])

# K2: a_K2 theta epsilon_ijk g_i d_j u_k   (d_j u_k = G_kj)
a_K2 = sp.symbols("a_K2", real=True)
eps = lambda i, j, k: int(sp.LeviCivita(i, j, k))
K2 = 0
for i, j, k in product(range(3), repeat=3):
    K2 += a_K2 * theta * eps(i, j, k) * g[i] * G[k, j]
K2 = sp.expand(K2)
K2_r1p = sp.expand(K2.subs(P_subs).subs(R1_subs))
print("K2_ON_R1_P", K2_r1p)
print("K2_ODD_ENERGY", split_parity(K2_r1p)[1])

# G3: g -> g + a_G3 e3, held. Spurion g·u * theta is the simplest N15 channel.
a_G3 = sp.symbols("a_G3", real=True)
gG3 = sp.Matrix([g[0], g[1], g[2] + a_G3])
inv_gu_theta = sp.expand(sum(gG3[i] * u[i] for i in range(3)) * theta)
inv_r1p = sp.expand(inv_gu_theta.subs(P_subs).subs(R1_subs))
print("G3_G_DOT_U_THETA_ON_R1_P", inv_r1p)
print("G3_ODD_ENERGY", split_parity(inv_r1p)[1])
print("G3_ODD_IF_A_TRANSFORMED_TO_MINUS", split_parity(inv_r1p.subs(a_G3, -a_G3))[1])
