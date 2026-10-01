#!/usr/bin/env python3
"""K1/K2 as stored-energy additions: reflection parity, mixing on class P,
R1 interaction, and FORM vs COEFFICIENT.
"""
from __future__ import annotations

import sympy as sp

# Indices 0,1,2 <-> labels 1,2,3
u = sp.symbols("u1 u2 u3")
du = tuple(tuple(sp.symbols(f"d{a+1}u{i+1}") for a in range(3)) for i in range(3))
# du[i][a] = ∂_{a+1} u_{i+1}
theta = sp.Symbol("theta")
dtheta = sp.symbols("d1theta d2theta d3theta")
g = sp.symbols("g1 g2 g3")
a_K1, a_K2, beta = sp.symbols("a_K1 a_K2 beta")
ehat = (sp.sin(beta), 0, sp.cos(beta))

# Class P: ∂3 = 0
P = {du[i][2]: 0 for i in range(3)}
P[dtheta[2]] = 0
# R1: g2 = g3 = 0
R1 = {g[1]: 0, g[2]: 0}

K1 = a_K1 * sum(ehat[i] * sum(du[i][k] * dtheta[k] for k in range(3)) for i in range(3))
eps = sp.LeviCivita
K2 = a_K2 * theta * sum(
    eps(i, j, k) * g[i] * du[k][j] for i in range(3) for j in range(3) for k in range(3)
)
K1 = sp.expand(K1)
K2 = sp.expand(K2)
print("K1", K1)
print("K2", K2)

# Reflection 3 -> -3 of polar tensors
R = {
    u[0]: u[0], u[1]: u[1], u[2]: -u[2],
    theta: theta,
    dtheta[0]: dtheta[0], dtheta[1]: dtheta[1], dtheta[2]: -dtheta[2],
    g[0]: g[0], g[1]: g[1], g[2]: -g[2],
    # ê is a fixed axis written in components; if it is a polar background
    # vector it must flip e3. Two treatments:
}
for i in range(3):
    for a in range(3):
        sign = 1
        if i == 2:
            sign *= -1
        if a == 2:
            sign *= -1
        R[du[i][a]] = sign * du[i][a]

print("K1_RESIDUAL_E_FIXED", sp.expand(K1.subs(R) - K1))
print("K2_RESIDUAL_G_POLAR", sp.expand(K2.subs(R) - K2))

# If ê is treated as a polar vector, e3 -> -e3, then K1 is even.
R_polar_e = dict(R)
# ehat enters as numbers sinβ, 0, cosβ. Polar transformation of ê:
# replace cosβ -> -cosβ in the transformed energy, equivalent to evaluating
# K1 with e3 flipped then comparing to K1 with original e.
K1_e_flipped = sp.expand(K1.subs({sp.cos(beta): -sp.cos(beta)}))
print("K1_RESIDUAL_IF_E_POLAR_FLIPPED_ONLY", sp.expand(K1_e_flipped - K1))
K1_full_polar = sp.expand(K1.subs(R).subs({sp.cos(beta): -sp.cos(beta)}))
print("K1_RESIDUAL_E_AS_POLAR_VECTOR", sp.expand(K1_full_polar - K1))

K1_P = sp.expand(K1.subs(P))
K2_P = sp.expand(K2.subs(P))
K2_P_R1 = sp.expand(K2_P.subs(R1))
print("K1_ON_P", K1_P)
print("K2_ON_P_R1", K2_P_R1)

# Mixed odd-even: terms linear in (u3 jets) and linear in (theta jets or theta)
odd_jets = (du[2][0], du[2][1], u[2])
even_src = (theta, dtheta[0], dtheta[1], u[0], u[1], du[0][0], du[0][1], du[1][0], du[1][1])

def mixed_linear(expr):
    pieces = []
    for o in odd_jets:
        for e in even_src:
            c = sp.diff(expr, o, e)
            if c != 0:
                pieces.append((str(o), str(e), c))
    return pieces

print("K1_MIXED_ON_P", mixed_linear(K1_P))
print("K2_MIXED_ON_P_R1", mixed_linear(K2_P_R1))
print("K1_MIXED_IF_BETA_PI_2", mixed_linear(sp.expand(K1_P.subs(beta, sp.pi / 2))))
print("K1_MIXED_IF_R1_ZEROS_E3", mixed_linear(sp.expand(K1_P.subs(sp.cos(beta), 0))))

# Baseline Kronecker energy has no Levi-Civita: K2 is outside the family (FORM).
# Coefficient rescale of μ_R |curl u|^2 cannot produce K2.
curl2 = (du[2][1] - du[1][2]) ** 2 + (du[0][2] - du[2][0]) ** 2 + (du[1][0] - du[0][1]) ** 2
curl2_P = sp.expand(curl2.subs(P))
print("CURL2_ON_P", curl2_P)
print("CURL2_MIXED_ON_P", mixed_linear(curl2_P))

# K1 with ê not parallel to g is a new preferred-axis structure (FORM).
# If ê || e1 (β=π/2), K1 lives in the even sector only.
print("K1_EVEN_SECTOR_AT_BETA_PI2", sp.expand(K1_P.subs(beta, sp.pi / 2)))

# μ_θ = δU/δθ holding u jets fixed: K2 contributes a_K2 * ε_ijk g_i ∂j u_k
mu_theta_K2 = sp.diff(K2_P_R1, theta)
print("MU_THETA_FROM_K2", mu_theta_K2)
print("MU_THETA_K2_HAS_D2U3", mu_theta_K2.has(du[2][1]))
# K1: U contains a_K1 ê_i (∂k ui)(∂k θ), so μ_θ involves -∂k (a ê_i ∂k ui) after IBP;
# the held-fixed variational derivative wrt θ (not θ-jets) of the density as written
# (before IBP) is 0 because θ appears only through dtheta. After IBP it hits u3.
# Report both.
print("K1_DENSITY_D_DTHETA", sp.diff(K1_P, theta))
print("K1_DENSITY_D_D1THETA", sp.diff(K1_P, dtheta[0]))
print("K1_DENSITY_D_D1THETA_HAS_D1U3", sp.diff(K1_P, dtheta[0]).has(du[2][0]))
