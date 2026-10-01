#!/usr/bin/env python3
"""K1/K2 on class P: which added terms actually break z-reflection and
open odd-even blocks, and which are inert or vanishing.
"""
import sympy as sp

sinb, cosb = sp.symbols("sin_beta cos_beta", real=True)
ehat = (0, sinb, cosb)  # (x,y,z) as in directive B
gy = sp.symbols("g_y", real=True)
g = (0, gy, 0)

# z-independent jets
ux, uy, uz = sp.symbols("u_x u_y u_z")
th, ew = sp.symbols("theta e_W")
# derivatives: no z derivatives
dx_ux, dy_ux = sp.symbols("dx_ux dy_ux")
dx_uy, dy_uy = sp.symbols("dx_uy dy_uy")
dx_uz, dy_uz = sp.symbols("dx_uz dy_uz")
dx_th, dy_th = sp.symbols("dx_th dy_th")
dx_ew, dy_ew = sp.symbols("dx_ew dy_ew")
dz = lambda _f: 0  # class P

du = (
    (dx_ux, dy_ux, 0),
    (dx_uy, dy_uy, 0),
    (dx_uz, dy_uz, 0),
)
dth = (dx_th, dy_th, 0)
dew = (dx_ew, dy_ew, 0)
u = (ux, uy, uz)

odd = {uz, dx_uz, dy_uz}
even = {ux, uy, dx_ux, dy_ux, dx_uy, dy_uy, th, dx_th, dy_th, ew, dx_ew, dy_ew, gy, sinb, cosb}


def odd_even_mixed(expr):
    expr = sp.expand(expr)
    mixed = []
    for term in sp.Add.make_args(expr):
        odd_deg = 0
        even_field_deg = 0
        for s in odd:
            if term.has(s):
                odd_deg += sp.Poly(sp.expand(term), s).degree()
        even_dyn = {ux, uy, dx_ux, dy_ux, dx_uy, dy_uy, th, dx_th, dy_th, ew, dx_ew, dy_ew}
        for s in even_dyn:
            if term.has(s):
                even_field_deg += sp.Poly(sp.expand(term), s).degree()
        if odd_deg % 2 == 1 and even_field_deg >= 1:
            mixed.append(term)
    return sp.expand(sum(mixed)) if mixed else sp.Integer(0)


def report(name, expr):
    expr = sp.expand(expr)
    mixed = odd_even_mixed(expr)
    print(f"\nTERM[{name}]")
    print(f"  expr  = {expr}")
    print(f"  mixed_odd_even = {mixed}")
    print(f"  is_zero = {expr == 0}")
    print(f"  mixed_is_zero = {mixed == 0}")


print("=== K1 candidate terms, ê=(0, sinβ, cosβ), class P (d_z=0) ===")

# A. ê contracted into derivative index only: (ê · ∇φ)^2
e_dot_grad_th = ehat[0] * dth[0] + ehat[1] * dth[1] + ehat[2] * dth[2]
report("K1_e_dot_grad_theta_squared", e_dot_grad_th**2)

e_dot_grad_u_sq = sum(
    (ehat[0] * du[a][0] + ehat[1] * du[a][1] + ehat[2] * du[a][2]) ** 2 for a in range(3)
)
report("K1_e_dot_grad_each_u_squared", e_dot_grad_u_sq)

# B. ê contracted into field index of a gradient: ê_i ê_j d_k u_i d_k u_j
k1_field_index = 0
for k in range(3):
    e_dot_dk_u = ehat[0] * du[0][k] + ehat[1] * du[1][k] + ehat[2] * du[2][k]
    k1_field_index += e_dot_dk_u**2
report("K1_e_i_e_j_dk_ui_dk_uj", k1_field_index)

# C. ê · (grad u) · grad theta : ê_i d_k u_i d_k theta
k1_u_theta = 0
for k in range(3):
    e_dot_dk_u = ehat[0] * du[0][k] + ehat[1] * du[1][k] + ehat[2] * du[2][k]
    k1_u_theta += e_dot_dk_u * (dth[k])
report("K1_e_i_dk_ui_dk_theta", k1_u_theta)

# D. undifferentiated (ê · u) theta  -- may not qualify as 'field gradients'
report("K1_e_dot_u_times_theta", (ehat[0] * u[0] + ehat[1] * u[1] + ehat[2] * u[2]) * th)

# E. ê · ∇W already along y: (ê · g)^2 is background-only
report("K1_e_dot_g_squared_background_only", (ehat[0] * g[0] + ehat[1] * g[1] + ehat[2] * g[2]) ** 2)

print("\n=== K2 Levi-Civita candidates ===")
eps = lambda i, j, k: int(sp.LeviCivita(i, j, k))

# vanishing: eps_ijk g_i u_j u_k
k2_guu = sum(eps(i, j, k) * g[i] * u[j] * u[k] for i in range(3) for j in range(3) for k in range(3))
report("K2_eps_g_u_u_VANISHES", k2_guu)

# vanishing: eps_ijk g_i g_j anything
k2_gg = sum(eps(i, j, k) * g[i] * g[j] * u[k] for i in range(3) for j in range(3) for k in range(3))
report("K2_eps_g_g_u_VANISHES", k2_gg)

# g · curl u * theta = eps_ijk g_i d_j u_k * theta
curl_comp = [
    sum(eps(i, j, k) * du[k][j] for j in range(3) for k in range(3))
    for i in range(3)
]
# curl_i = eps_ijk d_j u_k
curl = [
    sum(eps(i, j, k) * du[k][j] for j in range(3) for k in range(3))
    for i in range(3)
]
g_dot_curl = g[0] * curl[0] + g[1] * curl[1] + g[2] * curl[2]
report("K2_g_dot_curl_u_times_theta", g_dot_curl * th)

# eps_ijk d_i u_j d_k theta  (no g)
k2_u_th = sum(
    eps(i, j, k) * du[j][i] * dth[k] for i in range(3) for j in range(3) for k in range(3)
)
report("K2_eps_du_dtheta_no_g", k2_u_th)

# helicity eps_ijk u_i d_j u_k (u-u, not u-scalar)
helicity = sum(
    eps(i, j, k) * u[i] * du[k][j] for i in range(3) for j in range(3) for k in range(3)
)
report("K2_helicity_u_curl_u", helicity)

print("\n=== ê_z SLOT ON CLASS P ===")
print("ehat · nabla = sinb * d_y + cosb * d_z = sinb * d_y   on P")
print("Any ê contraction into a derivative index has inert ê_z on class P.")
print("Mixing through ê_z requires ê contracted into a field (vector) index.")
