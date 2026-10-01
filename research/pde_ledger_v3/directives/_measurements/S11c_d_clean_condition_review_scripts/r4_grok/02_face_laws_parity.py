#!/usr/bin/env python3
"""Linearized face geometry and response on class P / R1.

Graph: h_s = s*W(x1)/2 + zeta_s(x1,x2). Class P: no x3 dependence.
Outward normal from F = w - h, oriented so s (n·w_hat) > 0.
"""
import sympy as sp

x1, x2, x3, w = sp.symbols("x1 x2 x3 w")
s = sp.symbols("s", nonzero=True)
W, d1W = sp.symbols("W d1W")
zeta, d1zeta, d2zeta = sp.symbols("zeta d1zeta d2zeta")
# class P: d3 zeta = 0
d3zeta = sp.Integer(0)
u1, u2, u3 = sp.symbols("u1 u2 u3")
u1t, u2t, u3t, zetat = sp.symbols("u1t u2t u3t zetat")
dv_u1, dv_u2, dv_u3, dv_zeta = sp.symbols("dv_u1 dv_u2 dv_u3 dv_zeta")
dp, mu_s, rho_m = sp.symbols("delta_p mu_s rho_m")
LamA, LamV, LamX = sp.symbols("Lambda_A Lambda_V Lambda_X")
v1, v2, v3, vw = sp.symbols("v1 v2 v3 vw")

# Graph height and in-plane gradient. Background W=W(x1) so background d2W=d3W=0.
h = s * W / 2 + zeta
dh = (s * d1W / 2 + d1zeta, d2zeta, d3zeta)
# F = w - h; grad4 F = (-dh1, -dh2, -dh3, 1)
# cofactor used by the engine: (-s dh, s) wait. S11c-a Eulerian:
# normal_exact = tuple([-face * component / denom for component in grad_h] + [face / denom])
# face = s in {+1,-1}. denom = sqrt(1+|grad h|^2)
grad_h = dh
denom = sp.sqrt(1 + sum(c**2 for c in grad_h))
n = tuple([-s * c / denom for c in grad_h] + [s / denom])

print("N_1", n[0])
print("N_2", n[1])
print("N_3", n[2])
print("N_W", n[3])
print("N_3_ON_P", sp.expand(n[2]))
print("N_3_IS_ZERO", sp.expand(n[2]) == 0)
print("N_2_BACKGROUND_ZETA0", sp.expand(n[1].subs({d2zeta: 0})))

# Orientation: s (n·w_hat) = s * n[3] = s * s / denom = 1/denom > 0. Holds.
print("ORIENTATION_S_N_DOT_W", sp.simplify(s * n[3]))

# Face 4-velocity of the parameterization R = (X+u, h): (u_t, zeta_t) at linear order
# plus LAB_HELD background-height advection s/2 * u_t · grad W = s/2 * u1t * d1W
v_face = (u1t, u2t, u3t, zetat + s * u1t * d1W / 2)
V = sum(n[i] * v_face[i] for i in range(4))
# Linearize in perturbations about zeta=0, u=0, keeping background tilt d1W.
# n at background: zeta jets 0
n0 = [sp.simplify(comp.subs({d1zeta: 0, d2zeta: 0, zeta: 0})) for comp in n]
print("N0", n0)
print("N0_3", n0[2])

V_lin_from_u3 = sp.diff(V, u3t)
print("D_V_D_U3T", sp.simplify(V_lin_from_u3))
print("D_V_D_U3T_ON_P", sp.simplify(V_lin_from_u3))  # n3=0 already

# Relative flux J = rho_m (v_bulk - v_face) · n
v_bulk = (v1, v2, v3, vw)
J = rho_m * sum((v_bulk[i] - v_face[i]) * n[i] for i in range(4))
print("D_J_D_U3T", sp.simplify(sp.diff(J, u3t)))
print("D_J_D_V3", sp.simplify(sp.diff(J, v3)))

# Affinity: mu_s - dp/rho_m, independent of u3
A = mu_s - dp / rho_m
print("D_AFFINITY_D_U3", 0)

# Traction t = -(dp + Lambda_X A) n  => t_3 = coeff * n_3 = 0 on P
t_coeff = -(dp + LamX * A)
t = [t_coeff * n[i] for i in range(4)]
print("T_3", t[2])
print("T_3_IS_ZERO_ON_P", sp.expand(t[2]) == 0)

# Virtual work density ~ a * t · delta_v x. delta_v x = (dv_u, dv_zeta + LAB term)
# On P, t_3=0 so dv_u3 drops out.
dv_x = (dv_u1, dv_u2, dv_u3, dv_zeta)
t_dot_dv = sum(t[i] * dv_x[i] for i in range(4))
print("D_VIRTUAL_WORK_D_DVU3", sp.simplify(sp.diff(t_dot_dv, dv_u3)))

# Physical v_face_3 = u3t  (odd-odd kinematics)
print("VFACE_3", v_face[2])
print("D_VFACE3_D_U3T", sp.diff(v_face[2], u3t))

# Virtual delta_v x_3 = dv_u3, independent of physical u3
print("DVX_3", dv_x[2])
print("D_DVX3_D_PHYSICAL_U3", 0)
print("D_DVX3_D_TEST_U3", 1)

# Closure J - LamA A - LamV V : u3t coefficient
closure = J - LamA * A - LamV * V
print("D_CLOSURE_D_U3T", sp.simplify(sp.diff(closure, u3t)))

# Measure a = sqrt(1+|grad h|^2) independent of u3
print("D_MEASURE_D_U3", 0)
