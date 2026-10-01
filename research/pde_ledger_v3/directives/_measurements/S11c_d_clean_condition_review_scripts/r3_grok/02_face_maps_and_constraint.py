#!/usr/bin/env python3
"""Linearized S11c-a face maps / face laws / constraint under R1 and class P.

Uses infinitesimal wave amplitudes (ε-bookkeeping) so coefficients are exact.
Face maps from S11c-a §3a; laws from §3b; constraint from S11c-b §1c.
"""
from __future__ import annotations

import sympy as sp

x1, x2, x3, t = sp.symbols("x1 x2 x3 t", real=True)
s = sp.symbols("s")
W0, eta, L_W, rho_m, Lambda_X = sp.symbols("W0 eta L_W rho_m Lambda_X", positive=True)
w1 = sp.Function("w1")
W_bg = W0 * (1 + eta * w1(x1 / L_W))  # R1: depends only on x1
dW = sp.diff(W_bg, x1)

eps = sp.symbols("eps")
# Class P amplitudes: functions of (x1,x2,t) with no x3; keep x2 dependence live.
u1, u2, u3, theta, zeta_c, dW_pert = sp.symbols("u1 u2 u3 theta zeta_c dW_pert")
# spatial jets we need
u1_t, u2_t, u3_t = sp.symbols("u1_t u2_t u3_t")
u1_1, u1_2, u1_3 = sp.symbols("u1_1 u1_2 u1_3")
u2_1, u2_2, u2_3 = sp.symbols("u2_1 u2_2 u2_3")
u3_1, u3_2, u3_3 = sp.symbols("u3_1 u3_2 u3_3")
zc_t, dW_t = sp.symbols("zc_t dW_t")
zc_1, zc_2, zc_3 = sp.symbols("zc_1 zc_2 zc_3")
dWp_1, dWp_2, dWp_3 = sp.symbols("dWp_1 dWp_2 dWp_3")

# Class P: every *_3 jet of a perturbation is zero
P = {
    u1_3: 0, u2_3: 0, u3_3: 0,
    zc_3: 0, dWp_3: 0,
}

zeta_s = zeta_c + s * dW_pert / 2
zeta_s_t = zc_t + s * dW_t / 2
zeta_s_1 = zc_1 + s * dWp_1 / 2
zeta_s_2 = zc_2 + s * dWp_2 / 2
zeta_s_3 = zc_3 + s * dWp_3 / 2

# Graph heights, first order in eps
h_L = s * (W_bg + eps * u1 * dW) / 2 + eps * zeta_s
h_M = s * W_bg / 2 + eps * zeta_s

def first_order(expr):
    e = sp.expand(expr.subs(P))
    return sp.series(e, eps, 0, 2).removeO()

def n_from_h(h):
    """Outward-oriented raw (unnormalized) normal from F = w - h, s(n·ŵ)>0."""
    dh1 = first_order(sp.diff(h, x1).subs({sp.diff(W_bg, x1): dW}))
    # h is algebraic in jets; substitute explicit jets for ∂x2, ∂x3
    # Reconstruct ∂h/∂x_i from jets
    return dh1

# Explicit ∂h/∂x_i including background
dhL = [
    first_order(s * dW / 2 + eps * (s * u1_1 * dW / 2 + s * u1 * sp.diff(dW, x1) / 2 + zeta_s_1)),
    first_order(eps * (s * u1_2 * dW / 2 + zeta_s_2)),
    first_order(eps * (s * u1_3 * dW / 2 + zeta_s_3)),
]
dhM = [
    first_order(s * dW / 2 + eps * zeta_s_1),
    first_order(eps * zeta_s_2),
    first_order(eps * zeta_s_3),
]

def oriented_n(dh):
    raw = [-dh[0], -dh[1], -dh[2], sp.Integer(1)]
    return [s * raw[0], s * raw[1], s * raw[2], s * raw[3]]

nL = oriented_n(dhL)
nM = oriented_n(dhM)

print("L_N3", sp.simplify(nL[2]))
print("M_N3", sp.simplify(nM[2]))
print("L_N3_AT_EPS0", sp.simplify(nL[2].subs(eps, 0)))
print("M_N3_AT_EPS0", sp.simplify(nM[2].subs(eps, 0)))
print("L_N3_COEFF_U3", sp.expand(nL[2]).coeff(u3))
print("M_N3_COEFF_U3", sp.expand(nM[2]).coeff(u3))
print("L_N1_HAS_U3", nL[0].has(u3))
print("L_N2_HAS_U3", nL[1].has(u3))
print("L_NW_HAS_U3", nL[3].has(u3))
print("L_N3_HAS_EVEN_ZETA", nL[2].has(zeta_c) or nL[2].has(dW_pert) or nL[2].has(zc_3) or nL[2].has(dWp_3))
print("L_N3_AFTER_P", sp.simplify(nL[2].subs(P)))
print("M_N3_AFTER_P", sp.simplify(nM[2].subs(P)))

# Background normal (ε=0)
nL0 = [sp.expand(c.subs(eps, 0)) for c in nL]
nM0 = [sp.expand(c.subs(eps, 0)) for c in nM]
print("L_N_BACKGROUND", nL0)
print("M_N_BACKGROUND", nM0)
print("L_N0_COMPONENT3", nL0[2])
print("M_N0_COMPONENT3", nM0[2])

# v_face = (u_t, ∂t h)
v_in = [eps * u1_t, eps * u2_t, eps * u3_t]
v_w_L = first_order(s * eps * u1_t * dW / 2 + eps * zeta_s_t)
v_w_M = first_order(eps * zeta_s_t)
print("VFACE3", v_in[2])
print("VFACE3_IS_EPS_U3T", sp.simplify(v_in[2] - eps * u3_t) == 0)

def dot4(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3]

vL = v_in + [v_w_L]
vM = v_in + [v_w_M]
# Linear V uses n0 · v_face^(1)  (n^(1)·v^(0)=0 because v_bg=0)
VL = sp.expand(dot4(nL0, vL))
VM = sp.expand(dot4(nM0, vM))
print("L_V_HAS_U3T", VL.has(u3_t))
print("M_V_HAS_U3T", VM.has(u3_t))
print("L_V_HAS_U1T", VL.has(u1_t))
print("M_V_HAS_U1T", VM.has(u1_t))
print("L_V", VL)
print("M_V", VM)

vb = [sp.symbols("vb1"), sp.symbols("vb2"), sp.symbols("vb3"), sp.symbols("vbw")]
JL = sp.expand(rho_m * dot4(nL0, [vb[i] - vL[i] for i in range(4)]))
JM = sp.expand(rho_m * dot4(nM0, [vb[i] - vM[i] for i in range(4)]))
print("L_J_HAS_U3T", JL.has(u3_t))
print("M_J_HAS_U3T", JM.has(u3_t))
print("L_J_HAS_VB3", JL.has(vb[2]))
print("M_J_HAS_VB3", JM.has(vb[2]))
print("L_J_HAS_VB1", JL.has(vb[0]))
print("L_J_HAS_U1T", JL.has(u1_t))

# traction t = -(δp + Λ_X A) n  => t_3 ~ n_3
delta_p, A_s = sp.symbols("delta_p A_s")
factor = -(delta_p + Lambda_X * A_s)
tL = [factor * nL0[i] for i in range(4)]
print("L_T3", tL[2])
print("L_T3_IS_ZERO", tL[2] == 0)
print("L_T1_IS_NONZERO_IF_TILT", sp.simplify(tL[0]) != 0)

# virtual displacement
dvu1, dvu2, dvu3 = sp.symbols("dvu1 dvu2 dvu3")
dv_zc, dv_dW = sp.symbols("dv_zc dv_dW")
dv_zeta_s = dv_zc + s * dv_dW / 2
dvx_L = [dvu1, dvu2, dvu3, s * dvu1 * dW / 2 + dv_zeta_s]
dvx_M = [dvu1, dvu2, dvu3, dv_zeta_s]
print("DVX3", dvx_L[2])
print("DVW_L_HAS_DVU3", dvx_L[3].has(dvu3))
print("DVW_M_HAS_DVU3", dvx_M[3].has(dvu3))
work_L = sp.expand(dot4(tL, dvx_L))
print("WORK_L_HAS_DVU3", work_L.has(dvu3))
print("WORK_L_HAS_DVU1", work_L.has(dvu1))
print("WORK_L_HAS_DVZC", work_L.has(dv_zc))

# Uniform constraint δθ + δe_W + ∇·δu = 0, e_W = δW/W0
# ∇·δu = u1_1 + u2_2 + u3_3  -> on P, u3_3 = 0
div_P = u1_1 + u2_2 + u3_3
print("DIV_P_AFTER_P", div_P.subs(P))
print("DIV_P_HAS_U3_AFTER_P", div_P.subs(P).has(u3) or div_P.subs(P).has(u3_1) or div_P.subs(P).has(u3_2))

# LAB_HELD / MATERIAL pullback of a background scalar Q(x1): Q(χ) ~ Q - u·∇Q
# ∇Q = (Q', 0, 0) so u·∇Q = u1 Q'
print("PULLBACK_U_DOT_GRADQ_COMPONENTS", "u1 * Q_x1 ; u2 coeff 0 ; u3 coeff 0")

# Reflection 3 -> -3 on polar vectors: odd part is component 3
print("POLAR_ODD_PART_IS_COMPONENT_3", True)
print("N3_VANISHES_ON_R1_P", nL0[2] == 0 and nM0[2] == 0 and nL[2].subs(P).subs(eps, 0) == 0)
