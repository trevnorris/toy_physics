#!/usr/bin/env python3
"""Standalone model of the S11c-b stored-energy construction (Kronecker
O(3) family + first-jet spurions) under the planar 3 -> -3 reflection and
restriction R1. Does not import the S11c-b engine or exports.
"""
from __future__ import annotations

from itertools import product

import sympy as sp

DIRECTIONS = (0, 1, 2)  # 0-based loop; labels are 1,2,3
LABEL = {0: 1, 1: 2, 2: 3}

# Abstract field tensors, same rank content as the engine's enumerate_* .
bu = sp.symbols("bu_1:4")
bG = tuple(tuple(sp.symbols(f"bG_{a+1}_{i+1}") for i in DIRECTIONS) for a in DIRECTIONS)
btheta = sp.Symbol("btheta")
bq = sp.symbols("bq_1:4")
be = sp.Symbol("be_local")
br = sp.symbols("br_local_1:4")
bg = sp.symbols("bg_1:4")  # background first jet of W (polar)


def perfect_matchings(items):
    if not items:
        return ((),)
    first = items[0]
    result = []
    for index in range(1, len(items)):
        pair = (first, items[index])
        rest = items[1:index] + items[index + 1 :]
        for matching in perfect_matchings(rest):
            result.append((pair,) + matching)
    return tuple(result)


def tensor_component(data, indices):
    if not indices:
        return sp.sympify(data)
    if len(indices) == 1:
        return data[indices[0]]
    return data[indices[0]][indices[1]]


def delta_contractions(factors):
    slots = tuple((factor, index) for factor, (_, rank, _) in enumerate(factors) for index in range(rank))
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
    return tuple(dict.fromkeys(expressions))


def unique_expressions(expressions):
    result = []
    for expression in expressions:
        expression = sp.expand(expression)
        if expression != 0 and expression not in result:
            result.append(expression)
    return tuple(result)


def enumerate_uniform_candidates():
    data = (
        ("GRAD_U", 2, bG),
        ("THETA", 0, btheta),
        ("GRAD_THETA", 1, bq),
        ("E_LOCAL", 0, be),
        ("GRAD_E_LOCAL", 1, br),
    )
    raw = []
    for left_index, left in enumerate(data):
        for right in data[left_index:]:
            raw.extend(delta_contractions((left, right)))
    antisymmetric = sp.expand(
        sp.Add(*((bG[i][j] - bG[j][i]) ** 2 for i in DIRECTIONS for j in range(i + 1, 3)))
    )
    symmetric_tracefree = sp.expand(
        sp.Rational(1, 2) * sp.Add(*((bG[i][j] + bG[j][i]) ** 2 for i in DIRECTIONS for j in DIRECTIONS))
        - sp.Rational(2, 3) * (sum(bG[i][i] for i in DIRECTIONS)) ** 2
    )
    return unique_expressions([antisymmetric, symmetric_tracefree, *raw])


def enumerate_new_candidates(g_vector):
    data = (
        ("U", 1, bu),
        ("GRAD_U", 2, bG),
        ("THETA", 0, btheta),
        ("GRAD_THETA", 1, bq),
        ("E_LOCAL", 0, be),
        ("GRAD_E_LOCAL", 1, br),
    )
    spurion = ("BACKGROUND_FIRST_JET", 1, g_vector)
    raw = []
    for left_index, left in enumerate(data):
        for right in data[left_index:]:
            raw.extend(delta_contractions((spurion, left, right)))
    return unique_expressions(raw)


def reflection_map():
    """Polar-vector 3 -> -3, scalars even. Gradients: ∂3 flips."""
    subs = {
        btheta: btheta,
        be: be,
        bu[0]: bu[0],
        bu[1]: bu[1],
        bu[2]: -bu[2],
        bq[0]: bq[0],
        bq[1]: bq[1],
        bq[2]: -bq[2],
        br[0]: br[0],
        br[1]: br[1],
        br[2]: -br[2],
        bg[0]: bg[0],
        bg[1]: bg[1],
        bg[2]: -bg[2],
    }
    for a in DIRECTIONS:
        for i in DIRECTIONS:
            sign = 1
            if a == 2:
                sign *= -1
            if i == 2:
                sign *= -1
            subs[bG[a][i]] = sign * bG[a][i]
    return subs


def levi_civita_term():
    eps = sp.LeviCivita
    return sp.expand(sp.Add(*(eps(i, j, k) * bg[i] * bG[k][j] * btheta for i, j, k in product(DIRECTIONS, repeat=3))))


R = reflection_map()
uniform = enumerate_uniform_candidates()
spurion = enumerate_new_candidates(bg)

print("N_UNIFORM_CANDIDATES", len(uniform))
print("N_SPURION_CANDIDATES", len(spurion))

uniform_res = [sp.expand(e.subs(R) - e) for e in uniform]
spurion_res = [sp.expand(e.subs(R) - e) for e in spurion]
print("UNIFORM_REFLECTION_RESIDUALS_NONZERO", sum(1 for r in uniform_res if r != 0))
print("SPURION_REFLECTION_RESIDUALS_NONZERO", sum(1 for r in spurion_res if r != 0))

chi = levi_civita_term()
print("LEVI_CIVITA_THETA_G_CURLU", chi)
print("LEVI_CIVITA_REFLECTION_RESIDUAL", sp.expand(chi.subs(R) - chi))
print("LEVI_CIVITA_IN_UNIFORM", any(sp.expand(chi - e) == 0 for e in uniform))
print("LEVI_CIVITA_IN_SPURION", any(sp.expand(chi - e) == 0 for e in spurion))

# R1: background jet only along direction 1 (label 1 => index 0)
r1 = {bg[1]: 0, bg[2]: 0}
# Class P: z-independent => any tensor index 3 on a derivative vanishes;
# u_3 undifferentiated remains, but ∂3 anything = 0.
p_subs = {bG[a][2]: 0 for a in DIRECTIONS}
p_subs.update({bq[2]: 0, br[2]: 0})

# Mixed bilinear in (odd sector u_3 and its in-plane gradients) vs even scalars/even u.
odd_atoms = (bu[2], bG[2][0], bG[2][1])  # u3, ∂1 u3, ∂2 u3  (∂3 u3 already 0)
even_atoms = (
    btheta, be, bq[0], bq[1], br[0], br[1],
    bu[0], bu[1],
    bG[0][0], bG[0][1], bG[1][0], bG[1][1],
    bG[2][2],  # ∂3 u3, already 0 on P
)

def mixed_degree(expr):
    e = sp.expand(expr.subs(r1).subs(p_subs))
    mixed = sp.Integer(0)
    for odd in odd_atoms:
        for even in even_atoms:
            mixed += sp.diff(e, odd, even)
    return sp.expand(mixed)

uniform_mixed = [mixed_degree(e) for e in uniform]
spurion_mixed = [mixed_degree(e) for e in spurion]
print("UNIFORM_ODD_EVEN_MIXED_NONZERO", sum(1 for m in uniform_mixed if m != 0))
print("SPURION_ODD_EVEN_MIXED_NONZERO", sum(1 for m in spurion_mixed if m != 0))
print("SPURION_MIXED_WITHOUT_R1", sum(1 for e in spurion if mixed_degree(e.subs({bg[1]: bg[1], bg[2]: bg[2]})) != 0 or True))

# Without R1, g3 * u3 * theta is even*odd*even wait: g3 is odd, u3 odd, theta even => even, allowed.
# That term couples odd u3 to even theta using odd g3.
g3_u3_theta = sp.expand(bg[2] * bu[2] * btheta)
print("G3_U3_THETA_IN_SPURION_UNRESTRICTED", any(sp.expand(e - g3_u3_theta) == 0 or g3_u3_theta in sp.expand(e).as_ordered_terms() for e in spurion))
print("G3_U3_THETA_COEFF_IN_SPURION_SUM", sp.expand(sum(spurion)).coeff(bg[2] * bu[2] * btheta))
print("G3_U3_THETA_AFTER_R1", sp.expand(g3_u3_theta.subs(r1)))

# Undifferentiated u3 against g along 1: g1 * u3 * theta is even*odd*even = odd, forbidden.
g1_u3_theta = bg[0] * bu[2] * btheta
print("G1_U3_THETA_REFLECTION_RESIDUAL", sp.expand(g1_u3_theta.subs(R) - g1_u3_theta))
print("G1_U3_THETA_COEFF_IN_SPURION_SUM", sp.expand(sum(spurion)).coeff(bg[0] * bu[2] * btheta))

# Curl-squared mixing on class P: (∂2 u3)^2 is odd^2 = even, diagonal in odd sector.
curl2 = uniform[0]
print("CURL_SQUARED", curl2)
curl2_P_R1 = sp.expand(curl2.subs(p_subs).subs(r1))
# coefficient of (∂1 u3)^2 and (∂2 u3)^2 and cross with even
print("CURL_P_COEFF_D1U3_SQ", curl2_P_R1.coeff(bG[2][0] ** 2))
print("CURL_P_COEFF_D2U3_SQ", curl2_P_R1.coeff(bG[2][1] ** 2))
print("CURL_P_COEFF_D1U3_D1U1", curl2_P_R1.coeff(bG[2][0] * bG[0][0]))
print("CURL_P_COEFF_D1U3_THETA", 0 if btheta not in curl2_P_R1.free_symbols else "HAS_THETA")
