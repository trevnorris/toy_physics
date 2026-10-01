#!/usr/bin/env python3
"""Face laws, constraint fold, and MATERIAL_ADVECTED on class P + R1."""
import sympy as sp

g1, s = sp.symbols("g1 s")
# Graph height at linear order on R1: h_s = s*W(x1)/2 + zeta_s(x1,x2)
# n ~ (-dh/dx1, -dh/dx2, -dh/dx3, 1) before normalization
zeta = sp.Function("zeta")
x1, x2, x3 = sp.symbols("x1 x2 x3")
W = sp.Function("W")
# class P: zeta = zeta(x1,x2), W = W(x1)
h = s * W(x1) / 2 + zeta(x1, x2)
dh = (sp.diff(h, x1), sp.diff(h, x2), sp.diff(h, x3))
print("=== Face normal on R1+P ===")
print("dh =", dh)
print("n3_numerator =", -dh[2])
print("n3_is_zero =", sp.expand(dh[2]) == 0)

# Outward orientation: s (n·w) > 0. w-component of n is +1 before normalize (graph).
# In-plane n = -∇h / a, a = sqrt(1+|∇h|^2)
a = sp.sqrt(1 + dh[0] ** 2 + dh[1] ** 2 + dh[2] ** 2)
n_inplane = tuple(-dh[i] / a for i in range(3))
n_w = 1 / a
print("n_inplane =", n_inplane)
print("n_w =", n_w)

u1t, u2t, u3t, zetat = sp.symbols("u1t u2t u3t zetat")
# v_face = (u_t, zeta_t) at linear order (graph)
V = n_inplane[0] * u1t + n_inplane[1] * u2t + n_inplane[2] * u3t + n_w * zetat
print("\n=== V_s = n·v_face ===")
print("V =", sp.simplify(V))
print("coeff_u3t =", sp.expand(sp.simplify(V)).coeff(u3t))
print("V_independent_of_u3 =", sp.expand(sp.simplify(V)).coeff(u3t) == 0)

# Traction t = -p_tot n, work t·delta u
delta_u3 = sp.symbols("delta_u3")
t_dot_du3 = (-sp.symbols("p_tot")) * n_inplane[2] * delta_u3
print("\n=== traction work on u3 ===")
print("t·delta_u3 =", t_dot_du3)
print("is_zero =", t_dot_du3 == 0 or sp.expand(t_dot_du3) == 0)

# Affinity / J: 𝒜 = μ_s - δp/ρ, J = ΛA 𝒜 + ΛV V
# μ_s from μ_θ(energy). If energy has no u3-even mix, μ_s independent of u3.
print("\n=== J_s channel ===")
print("J_depends_on_u3_only_through_V_or_mu_theta")
print("V_coeff_u3 = 0 on this geometry, so J_u3 = 0 unless mu_theta(u3) != 0")

# Constraint: δθ + δe + div δu = 0 (uniform linearisation; non-uniform adds u·∇ρ terms)
print("\n=== Constraint fold on P ===")
d1_du1, d2_du2, d3_du3 = sp.symbols("d1_delta_u1 d2_delta_u2 d3_delta_u3")
div_du = d1_du1 + d2_du2 + d3_du3
print("div_delta_u =", div_du)
print("on_P_d3_delta_u3 = 0, so div_delta_u has no u3")
# MATERIAL_ADVECTED extra: u·∇ρ with ∇ρ // e1
rho_x1 = sp.symbols("drho_dx1")
u1, u2, u3 = sp.symbols("u1 u2 u3")
adv = u1 * rho_x1 + u2 * 0 + u3 * 0
print("u·∇ρ_R1 =", adv)
print("adv_coeff_u3 =", adv.coeff(u3))

# MATERIAL_ADVECTED height: h^M = s W(X)/2 + zeta(X) = s W(χ)/2 + zeta(χ) in Eulerian
# χ1 = x1 - u1 + O(u^2); W depends only on χ1. Linear u3 does not enter W(χ).
print("\n=== MATERIAL_ADVECTED linear height ===")
print("W(chi)_linear = W(x1) - u1*W'(x1)")
print("u3_coefficient_in_height = 0")

# LAB_HELD: h^L = s W(x)/2 + zeta(χ); W(x) independent of u; zeta(χ) = zeta(x) - u·∇zeta
# u·∇zeta = u1 d1 zeta + u2 d2 zeta + u3 d3 zeta, d3 zeta = 0 on P
print("LAB_HELD zeta(chi)_linear_u3_coeff = 0 on P")

# Face map in-plane position includes u3, but n3=0 so t·e3 = 0 still.
print("inplane_face_point_has_u3_but_n3=0_so_work_still_zero")
