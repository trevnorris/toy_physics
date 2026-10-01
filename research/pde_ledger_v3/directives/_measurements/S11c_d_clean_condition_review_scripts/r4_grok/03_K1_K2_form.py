#!/usr/bin/env python3
"""K1/K2 FORM controls on class P after R1.

K1: a_K1 e_hat_i (d_k u_i)(d_k theta), e_hat=(sin beta, 0, cos beta)
K2: a_K2 theta epsilon_ijk g_i d_j u_k, g = grad W
"""
import sympy as sp

aK1, aK2, beta = sp.symbols("a_K1 a_K2 beta", real=True)
g1, g2, g3 = sp.symbols("g1 g2 g3")
u1d = sp.symbols("u1_d1 u1_d2 u1_d3")
u2d = sp.symbols("u2_d1 u2_d2 u2_d3")
u3d = sp.symbols("u3_d1 u3_d2 u3_d3")
thd = sp.symbols("th_d1 th_d2 th_d3")
th = sp.symbols("theta")
u_jets = (u1d, u2d, u3d)
ehat = (sp.sin(beta), 0, sp.cos(beta))
g = (g1, g2, g3)

K1 = sp.Integer(0)
for i in range(3):
    for k in range(3):
        K1 += aK1 * ehat[i] * u_jets[i][k] * thd[k]
K1 = sp.expand(K1)

eps = sp.LeviCivita
K2 = sp.Integer(0)
for i in range(3):
    for j in range(3):
        for k in range(3):
            K2 += aK2 * th * eps(i, j, k) * g[i] * u_jets[k][j]
K2 = sp.expand(K2)

class_P = {u1d[2]: 0, u2d[2]: 0, u3d[2]: 0, thd[2]: 0}
R1 = {g2: 0, g3: 0}

K1_P = sp.expand(K1.subs(class_P))
K2_P = sp.expand(K2.subs(class_P))
K2_P_R1 = sp.expand(K2_P.subs(R1))

odd = {u3d[0], u3d[1], u3d[2]}
even_jets = {thd[0], thd[1], thd[2], th, u1d[0], u1d[1], u2d[0], u2d[1]}


def mixed_terms(expr):
    out = []
    for term in sp.Add.make_args(sp.expand(expr)):
        fs = term.free_symbols
        if (fs & odd) and (fs & even_jets):
            out.append(term)
    return out


print("K1", K1)
print("K1_ON_P", K1_P)
print("K1_MIXED_ON_P", mixed_terms(K1_P))
K1_if_e3_zeroed = sp.expand(K1_P.subs({sp.cos(beta): 0}))
print("K1_MIXED_IF_E3_ZEROED", mixed_terms(K1_if_e3_zeroed))
print("K1_EVEN_EVEN_IF_E3_ZEROED", K1_if_e3_zeroed)

print("K2", K2)
print("K2_ON_P_R1", K2_P_R1)
print("K2_MIXED_ON_P_R1", mixed_terms(K2_P_R1))

# K1 is FORM: e_hat is a new preferred axis not parallel to g=(g1,0,0)
# unless cos(beta)=0. Coefficient rescaling of mu_R would not produce e_hat_3.
print("K1_AXIS_PARALLEL_TO_G_ONLY_IF", "cos(beta)=0")
print("K2_OUTSIDE_KRONECKER_FAMILY", True)

# Reflection R3 on K1 energy density: u3 jets odd, theta even, ehat fixed (control datum)
# The ê_3 piece is odd in the fields (mixes sectors). The ê_1 piece is even-even.
print("K1_E3_PIECE", sp.expand(aK1 * sp.cos(beta) * (u3d[0] * thd[0] + u3d[1] * thd[1] + u3d[2] * thd[2]).subs(class_P)))
print("K1_E1_PIECE", sp.expand(aK1 * sp.sin(beta) * (u1d[0] * thd[0] + u1d[1] * thd[1]).subs(class_P)))
