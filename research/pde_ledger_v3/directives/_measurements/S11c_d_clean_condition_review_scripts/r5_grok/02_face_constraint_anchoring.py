#!/usr/bin/env python3
"""Face normal, material advection, constraint, and V_s / J_s parity on class P/R1.

Uses the S11c-a face-map equations, not the engine. Class P: d3 of perturbations
= 0. R1: W_bg = W_bg(x1) only.
"""
import sympy as sp

x1, x2, x3 = sp.symbols("x1 x2 x3", real=True)
s, W0, eta = sp.symbols("s W0 eta", real=True)
w1 = sp.Function("w1")
W_bg = W0 * (1 + eta * w1(x1))

u1, u2, u3 = sp.Function("u1"), sp.Function("u2"), sp.Function("u3")
zeta_c = sp.Function("zeta_c")
delta_W = sp.Function("delta_W")
theta = sp.Function("theta")

# Class P: no x3 dependence of perturbations.
args = (x1, x2)
u = (u1(*args), u2(*args), u3(*args))
zc = zeta_c(*args)
dW = delta_W(*args)
th = theta(*args)
zeta_s = zc + s * dW / 2

# LAB_HELD graph height (linear, Eulerian): h = s W_bg(x)/2 + zeta_s
h_L = s * W_bg / 2 + zeta_s
# MATERIAL_ADVECTED: h = s W_bg(chi)/2 + zeta_s(chi); linear: chi = x - u
# W_bg(x - u) = W_bg(x) - u·grad W  + O(u^2); grad W = (d1 W, 0, 0)
d1W = sp.diff(W_bg, x1)
h_M = s * (W_bg - u[0] * d1W) / 2 + zeta_s

eps = sp.symbols("eps")


def unit_normal(h):
    gh = (sp.diff(h, x1), sp.diff(h, x2), sp.diff(h, x3))
    denom = sp.sqrt(1 + sum(c**2 for c in gh))
    # S11c-a engine: n_i = -s * d_i h / denom, n_w = s / denom
    n = tuple(-s * c / denom for c in gh) + (s / denom,)
    return n


def linear_n3(h):
    n = unit_normal(h)
    # first order in wave amplitude: treat u, zeta, deltaW as O(eps)
    wave = {
        u1(*args): eps * u1(*args),
        u2(*args): eps * u2(*args),
        u3(*args): eps * u3(*args),
        zeta_c(*args): eps * zc,
        delta_W(*args): eps * dW,
    }
    n3 = n[2].subs(wave)
    return sp.simplify(sp.series(n3, eps, 0, 2).removeO())


print("N3_LAB_HELD_LINEAR", linear_n3(h_L))
print("N3_MATERIAL_LINEAR", linear_n3(h_M))

# Material pullback of a scalar Q(x1): Q(chi) = Q(x1 - u1) at linear order.
# Coefficient of u3:
Q = sp.Function("Q")
# Q(x - u) ~ Q(x) - u·grad Q, grad Q = (Q', 0, 0)
pull = -u[0] * sp.diff(Q(x1), x1)  # linear correction, no u3
print("MATERIAL_PULLBACK_LINEAR", pull)
print("MATERIAL_PULLBACK_COEFF_U3", sp.Integer(0))

# Uniform virtual constraint: d_v theta + d_v e_W + div d_v u = 0
# On class P, div u = d1 u1 + d2 u2; no u3.
div_u = sp.diff(u[0], x1) + sp.diff(u[1], x2) + sp.diff(u[2], x3)
print("DIV_U_ON_P", div_u)
print("DIV_U_HAS_U3", sp.diff(div_u, u3(*args)))

# V_s = n · v_face, v_face = (dt u, dt zeta_s)
# At background (n = (0,0,0,s) up to tilt along x1 from d1 W)
# Linear even scalar: n_bg · (dt u, dt zeta) = s * dt zeta_s plus tilt n1 * dt u1.
n_bg_L = unit_normal(s * W_bg / 2)
# evaluate at zero wave
n_bg = tuple(sp.simplify(c.subs({zeta_c(*args): 0, delta_W(*args): 0, u1(*args): 0})) for c in n_bg_L)
print("N_BG_LAB", tuple(sp.simplify(c) for c in n_bg))

# Linear V from background normal dotted into face velocity, plus n^(1) dotted into 0.
dt_u = sp.symbols("dt_u1 dt_u2 dt_u3")
dt_zeta = sp.symbols("dt_zeta_s")
v_face = dt_u + (dt_zeta,)
V_bg = sum(sp.simplify(n_bg[i]) * v_face[i] for i in range(4))
print("V_FROM_BACKGROUND_NORMAL", sp.simplify(V_bg))
print("V_COEFF_DT_U3", sp.diff(V_bg, dt_u[2]))

# J uses the same n. Relative (v_bulk - v_face)·n.
# On P, v_bulk_3 = d3 phi = 0 from a scalar potential independent of x3.
v_b = sp.symbols("vb1 vb2 vb3 vbw")
J_bg = sum(sp.simplify(n_bg[i]) * (v_b[i] - v_face[i]) for i in range(4))
print("J_FROM_BACKGROUND_NORMAL", sp.simplify(J_bg))
print("J_COEFF_VB3", sp.diff(J_bg, v_b[2]))
print("J_COEFF_DT_U3", sp.diff(J_bg, dt_u[2]))
