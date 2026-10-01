#!/usr/bin/env python3
"""Kronecker energy, pin-B constraint, and MATERIAL_ADVECTED mixing on class P/R1.

Standalone model of the S11c-b construction rules:
- O(3) Kronecker bilinears (no Levi-Civita)
- background spurion g = grad W, restricted to (g1, 0, 0)
- class P: d3 of every field = 0
Reflection R3: x3 -> -x3 acts on polar vectors by flipping the 3-component.
"""
from itertools import combinations, product
import sympy as sp

u1, u2, u3 = sp.symbols("u1 u2 u3")
th, eW, zc = sp.symbols("theta e_W zeta_c")
g1, g2, g3 = sp.symbols("g1 g2 g3")
# First jets; class P sets every *_d3 = 0 after construction of candidates.
u1_d1, u1_d2, u1_d3 = sp.symbols("u1_d1 u1_d2 u1_d3")
u2_d1, u2_d2, u2_d3 = sp.symbols("u2_d1 u2_d2 u2_d3")
u3_d1, u3_d2, u3_d3 = sp.symbols("u3_d1 u3_d2 u3_d3")
th_d1, th_d2, th_d3 = sp.symbols("th_d1 th_d2 th_d3")
e_d1, e_d2, e_d3 = sp.symbols("e_d1 e_d2 e_d3")

G = ((u1_d1, u1_d2, u1_d3), (u2_d1, u2_d2, u2_d3), (u3_d1, u3_d2, u3_d3))
U = (u1, u2, u3)
TH = th
Q = (th_d1, th_d2, th_d3)
E = eW
R = (e_d1, e_d2, e_d3)
g = (g1, g2, g3)


def tensor_component(data, indices):
    x = data
    for i in indices:
        x = x[i]
    return x


def perfect_matchings(slots):
    slots = tuple(slots)
    if not slots:
        yield ()
        return
    a = slots[0]
    for i, b in enumerate(slots[1:], start=1):
        rest = slots[1:i] + slots[i + 1 :]
        for m in perfect_matchings(rest):
            yield ((a, b),) + m


def delta_contractions(factors):
    slots = tuple((fi, si) for fi, (_, rank, _) in enumerate(factors) for si in range(rank))
    if len(slots) % 2:
        return ()
    out = []
    for matching in perfect_matchings(slots):
        total = sp.Integer(0)
        for assignment in product(range(3), repeat=len(matching)):
            index_map = {}
            for pair_index, pair in enumerate(matching):
                for slot in pair:
                    index_map[slot] = assignment[pair_index]
            term = sp.Integer(1)
            for factor_index, (_, rank, data) in enumerate(factors):
                indices = tuple(index_map[(factor_index, s)] for s in range(rank))
                term *= tensor_component(data, indices)
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
    return tuple(uniq)


uniform_data = (
    ("GRAD_U", 2, G),
    ("THETA", 0, TH),
    ("GRAD_THETA", 1, Q),
    ("E_LOCAL", 0, E),
    ("GRAD_E", 1, R),
)
uniform = []
for i, left in enumerate(uniform_data):
    for right in uniform_data[i:]:
        uniform.extend(delta_contractions((left, right)))
curl2 = sp.expand(sum((G[i][j] - G[j][i]) ** 2 for i in range(3) for j in range(i + 1, 3)))
uniform.append(curl2)

spurion = ("G", 1, g)
new_data = (("U", 1, U),) + uniform_data
new = []
for i, left in enumerate(new_data):
    for right in new_data[i:]:
        new.extend(delta_contractions((spurion, left, right)))

class_P = {
    u1_d3: 0, u2_d3: 0, u3_d3: 0, th_d3: 0, e_d3: 0,
}
R1 = {g2: 0, g3: 0}

odd_symbols = {u3, u3_d1, u3_d2, u3_d3}
even_symbols = {
    u1, u2, th, eW, zc, g1, g2, g3,
    u1_d1, u1_d2, u1_d3, u2_d1, u2_d2, u2_d3,
    th_d1, th_d2, th_d3, e_d1, e_d2, e_d3,
}


def is_mixed(expr):
    expr = sp.expand(expr)
    mixed = []
    if expr == 0:
        return mixed
    for term in sp.Add.make_args(expr):
        factors = term.free_symbols
        has_odd = bool(factors & odd_symbols)
        has_even_field = bool(factors & {u1, u2, th, eW, zc,
                                         u1_d1, u1_d2, u2_d1, u2_d2,
                                         th_d1, th_d2, e_d1, e_d2})
        # g is background, not a perturbation
        if has_odd and has_even_field:
            mixed.append(term)
    return mixed


print("UNIFORM_CANDIDATE_COUNT", len(set(sp.srepr(e) for e in uniform)))
print("SPURION_CANDIDATE_COUNT", len(set(sp.srepr(e) for e in new)))

mixed_uniform = []
for e in uniform:
    eP = sp.expand(e.subs(class_P))
    m = is_mixed(eP)
    if m:
        mixed_uniform.append((eP, m))

mixed_spurion_fullg = []
mixed_spurion_R1 = []
for e in new:
    eP = sp.expand(e.subs(class_P))
    m = is_mixed(eP)
    if m:
        mixed_spurion_fullg.append(eP)
    eR = sp.expand(eP.subs(R1))
    mR = is_mixed(eR)
    if mR:
        mixed_spurion_R1.append(eR)

print("UNIFORM_MIXED_ON_P_COUNT", len(mixed_uniform))
print("SPURION_MIXED_ON_P_FULL_G_COUNT", len(mixed_spurion_fullg))
print("SPURION_MIXED_ON_P_R1_COUNT", len(mixed_spurion_R1))
print("SPURION_MIXED_ON_P_R1_TERMS", mixed_spurion_R1[:8])

# Levi-Civita is NOT in the Kronecker family; exhibit one chiral scalar.
eps = sp.LeviCivita
chiral = sp.expand(sum(eps(i, j, k) * g[i] * (u1_d1, u2_d1, u3_d1)[k] * [u1, u2, u3][j]
                       for i in range(3) for j in range(3) for k in range(3)))
# simpler: theta * eps_ijk g_i d_j u_k
chiral = sp.Integer(0)
dU = ( (u1_d1, u1_d2, u1_d3), (u2_d1, u2_d2, u2_d3), (u3_d1, u3_d2, u3_d3) )
for i, j, k in product(range(3), repeat=3):
    chiral += eps(i, j, k) * g[i] * dU[k][j] * th
chiral = sp.expand(chiral)
print("CHIRAL_EPS_G_DU_THETA", chiral)
print("CHIRAL_ON_P_R1", sp.expand(chiral.subs(class_P).subs(R1)))
print("CHIRAL_MIXED_ON_P_R1", is_mixed(sp.expand(chiral.subs(class_P).subs(R1))))

# Pin B / virtual constraint on class P: div u = d1 u1 + d2 u2 + d3 u3
div_u = u1_d1 + u2_d2 + u3_d3
div_P = sp.expand(div_u.subs(class_P))
print("CONSTRAINT_DIV_ON_P", div_P)
print("CONSTRAINT_DIV_CONTAINS_U3", u3 in div_P.free_symbols or u3_d1 in div_P.free_symbols or u3_d2 in div_P.free_symbols or u3_d3 in div_P.free_symbols)
print("CONSTRAINT_MIXED_ON_P", is_mixed(div_P))

# MATERIAL_ADVECTED linear pullback: -u · grad Q, Q=Q(x1) so grad Q = (Q', 0, 0)
Qp = sp.symbols("d1_Q")
adv = -(u1 * Qp + u2 * 0 + u3 * 0)
print("MATERIAL_ADVECTED_LINEAR", adv)
print("MATERIAL_ADVECTED_MIXED", is_mixed(adv))

# LAB_HELD face height: W_bg(x1+u1); linear u3 coefficient
print("LAB_HELD_WBG_U3_COEFF", 0)

# Kinetic |dt u|^2 diagonal
print("KINETIC_U3_CROSS_EVEN", 0)

# Curl squared on P
c1 = u3_d2 - 0  # d2 u3 - d3 u2
c2 = 0 - u3_d1  # d3 u1 - d1 u3
c3 = u2_d1 - u1_d2
curl2_P = sp.expand(c1**2 + c2**2 + c3**2)
print("CURLSQ_ON_P", curl2_P)
print("CURLSQ_MIXED_ON_P", is_mixed(curl2_P))
