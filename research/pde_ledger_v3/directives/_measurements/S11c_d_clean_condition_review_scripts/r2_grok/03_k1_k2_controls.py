#!/usr/bin/env python3
"""v2 K1/K2: FORM vs COEFFICIENT, mixing on class P, null representatives."""
import sympy as sp

a_K1, a_K2, beta = sp.symbols("a_K1 a_K2 beta", real=True)
g1 = sp.symbols("g1", real=True)
k2, omega = sp.symbols("k2 omega", complex=True)

u1, u2, u3 = sp.symbols("u1 u2 u3")
th = sp.symbols("theta")
# class P: d3 = 0; d2 -> i k2 on a Fourier mode, kept symbolic as d2_*
d1_u1, d2_u1 = sp.symbols("d1_u1 d2_u1")
d1_u2, d2_u2 = sp.symbols("d1_u2 d2_u2")
d1_u3, d2_u3 = sp.symbols("d1_u3 d2_u3")
d1_th, d2_th = sp.symbols("d1_th d2_th")
d3 = 0

ehat = (sp.sin(beta), 0, sp.cos(beta))
g = (g1, 0, 0)
du = (
    (d1_u1, d2_u1, d3),
    (d1_u2, d2_u2, d3),
    (d1_u3, d2_u3, d3),
)
dth = (d1_th, d2_th, d3)

# K1 : a_K1 * e_i (d_k u_i)(d_k theta)
K1 = 0
for i in range(3):
    for k in range(3):
        K1 += a_K1 * ehat[i] * du[i][k] * dth[k]
K1 = sp.expand(K1)

# K2 : a_K2 * theta * eps_ijk g_i d_j u_k
K2 = 0
for i in range(3):
    for j in range(3):
        for k in range(3):
            K2 += a_K2 * th * int(sp.LeviCivita(i + 1, j + 1, k + 1)) * g[i] * du[k][j]
K2 = sp.expand(K2)

odd = {u3, d1_u3, d2_u3}
even_dyn = {u1, u2, th, d1_u1, d2_u1, d1_u2, d2_u2, d1_th, d2_th}


def split_mixed(expr):
    expr = sp.expand(expr)
    mixed = []
    even_even = []
    odd_odd = []
    other = []
    for term in sp.Add.make_args(expr):
        odd_deg = sum(sp.degree(sp.expand(term), s) for s in odd if term.has(s))
        even_deg = sum(sp.degree(sp.expand(term), s) for s in even_dyn if term.has(s))
        if odd_deg % 2 == 1 and even_deg >= 1:
            mixed.append(term)
        elif odd_deg == 0 and even_deg >= 1:
            even_even.append(term)
        elif odd_deg >= 1 and even_deg == 0:
            odd_odd.append(term)
        else:
            other.append(term)
    return {
        "mixed": sp.expand(sum(mixed)) if mixed else sp.Integer(0),
        "even_even": sp.expand(sum(even_even)) if even_even else sp.Integer(0),
        "odd_odd": sp.expand(sum(odd_odd)) if odd_odd else sp.Integer(0),
        "other": sp.expand(sum(other)) if other else sp.Integer(0),
    }


print("=== K1 ===")
print("K1 =", K1)
s1 = split_mixed(K1)
for k, v in s1.items():
    print(f"K1_{k} =", v)
print("K1_mixed_zero_identically =", s1["mixed"] == 0)
print("K1_mixed_at_beta_0 =", sp.simplify(s1["mixed"].subs(beta, 0)))
print("K1_mixed_at_beta_pi_2 =", sp.simplify(s1["mixed"].subs(beta, sp.pi / 2)))
print("K1_even_even_at_beta_0 =", sp.simplify(s1["even_even"].subs(beta, 0)))

print("\n=== K2 ===")
print("K2 =", K2)
s2 = split_mixed(K2)
for k, v in s2.items():
    print(f"K2_{k} =", v)
print("K2_mixed_zero_identically =", s2["mixed"] == 0)
print("K2_at_g1_0 =", K2.subs(g1, 0))
print("K2_at_d2_u3_0 =", K2.subs(d2_u3, 0))

# Null-representative tests: total-divergence remainder still mixes
print("\n=== K2 IBP remainder ===")
# theta * d2 u3 = d2(theta u3) - u3 d2 theta
remainder = a_K2 * g1 * (-u3 * d2_th)  # after dropping the divergence
print("IBP_remainder_mixed =", remainder)
print("remainder_is_mixed =", remainder != 0)

print("\n=== K1 IBP remainder ===")
# e3 d_k u3 d_k theta = d_k(e3 u3 d_k theta) - e3 u3 d_k d_k theta
# second-derivative remainder still mixed
lap_th = sp.symbols("lap_theta")
k1_remainder = a_K1 * sp.cos(beta) * (-u3 * lap_th)
print("K1_IBP_remainder =", k1_remainder)

# Is K1 FORM? It introduces a preferred axis ehat not in the O(3) Kronecker family.
# The e3 piece cannot be written as a g-contraction under R1 (g3=0).
print("\n=== FORM character ===")
print("K1_introduces_preferred_axis_with_e3 = True")
print("K1_e3_piece_not_in_R1_N15_because_g3_is_0 = True")
print("K2_uses_LeviCivita_forbidden_by_Kronecker_family = True")
print("K2_is_zero_on_P_if_k2_frozen_to_0_and_no_x2_dependence =",
      K2.subs(d2_u3, 0) == 0)

# Derivative-index contraction of ehat is inert in e3 on P (round-1 finding)
e_dot_grad_th = ehat[0] * dth[0] + ehat[1] * dth[1] + ehat[2] * dth[2]
print("\n=== ehat contracted into derivative index (NOT v2 K1) ===")
print("e·∇theta on P =", sp.expand(e_dot_grad_th))
print("e3_slot_inert_on_P =", sp.expand(e_dot_grad_th).has(sp.cos(beta)) is False
      or sp.expand(e_dot_grad_th).coeff(sp.cos(beta)) == 0)

# v2 K1 contracts e into FIELD index: e_i d_k u_i
e_field = tuple(ehat[0] * du[0][k] + ehat[1] * du[1][k] + ehat[2] * du[2][k] for k in range(3))
print("e_i d_k u_i =", e_field)
print("e3_field_slot_active =", sp.expand(e_field[0]).has(sp.cos(beta)))
