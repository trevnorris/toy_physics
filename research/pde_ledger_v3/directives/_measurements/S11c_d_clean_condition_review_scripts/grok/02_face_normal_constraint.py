#!/usr/bin/env python3
"""Linear face normal, V_s, traction pairing, and virtual-constraint content
for W_bg = W(y) and z-independent perturbations.

Reconstructed from S11c-a §§3a–3b and T-g, not from engine output.
"""
import sympy as sp

x, y, z, t = sp.symbols("x y z t", real=True)
W = sp.Function("W_bg")
ux, uy, uz = sp.symbols("u_x u_y u_z", cls=sp.Function)
th = sp.Function("theta")
zp, zm = sp.symbols("zeta_plus zeta_minus", cls=sp.Function)
# z-independent: fields(x,y,t) only
u = (ux(x, y, t), uy(x, y, t), uz(x, y, t))
zeta_p = zp(x, y, t)
zeta_m = zm(x, y, t)

# Graph heights, LAB_HELD background W(y), linear faces
h_p = sp.Rational(1, 2) * W(y) + zeta_p
h_m = -sp.Rational(1, 2) * W(y) + zeta_m

def graph_normal(h, s):
    """Outward unit normal from F = w - h, with s (n·w_hat) > 0.
    4-vector (n_x, n_y, n_z, n_w)."""
    dh_x = sp.diff(h, x)
    dh_y = sp.diff(h, y)
    dh_z = sp.diff(h, z)
    # ∇_4 F = (-dh_x, -dh_y, -dh_z, 1)
    raw = sp.Matrix([-dh_x, -dh_y, -dh_z, 1])
    # For s=+1 need n_w > 0: use +raw. For s=-1 need n_w < 0: use -raw.
    vec = raw if s == 1 else -raw
    a = sp.sqrt(sp.simplify(vec.dot(vec)))
    n = sp.simplify(vec / a)
    return n, a

n_p, a_p = graph_normal(h_p, +1)
n_m, a_m = graph_normal(h_m, -1)

print("=== BACKGROUND + LINEAR n_hat (exact graph formula) ===")
print(f"n_plus = {n_p}")
print(f"n_minus = {n_m}")

# Linearise in face displacements about zeta=0, keep W' live
eps = sp.symbols("eps_face")
subs0 = {
    zeta_p: 0,
    zeta_m: 0,
    sp.diff(zeta_p, x): 0,
    sp.diff(zeta_p, y): 0,
    sp.diff(zeta_m, x): 0,
    sp.diff(zeta_m, y): 0,
}

def lin_n(n):
    # n is already linear in zeta jets + W'; set zeta jets to 0 for background n
    return sp.simplify(n.subs(subs0))

n_p0 = lin_n(n_p)
n_m0 = lin_n(n_m)
print("\n=== BACKGROUND n_hat (zeta=0) ===")
print(f"n_plus0  = {n_p0.T}")
print(f"n_minus0 = {n_m0.T}")
print(f"n_plus0_z  = {n_p0[2]}")
print(f"n_minus0_z = {n_m0[2]}")

# n · u_inplane
n_p0_inplane = sp.Matrix([n_p0[0], n_p0[1], n_p0[2]])
n_m0_inplane = sp.Matrix([n_m0[0], n_m0[1], n_m0[2]])
u_col = sp.Matrix(list(u))
n_dot_u_p = sp.simplify(n_p0_inplane.dot(u_col))
n_dot_u_m = sp.simplify(n_m0_inplane.dot(u_col))
print("\n=== BACKGROUND n_inplane · u ===")
print(f"n_plus0 · u  = {n_dot_u_p}")
print(f"n_minus0 · u = {n_dot_u_m}")
print(f"COEFF_u_z_in_n_plus_dot_u  = {sp.expand(n_dot_u_p).coeff(uz(x, y, t))}")
print(f"COEFF_u_z_in_n_minus_dot_u = {sp.expand(n_dot_u_m).coeff(uz(x, y, t))}")

# Face velocity: v_face = (d_t u, d_t zeta); V = n · v_face
# At linear order use background n
dt_up = sp.Matrix([sp.diff(u[0], t), sp.diff(u[1], t), sp.diff(u[2], t), sp.diff(zeta_p, t)])
dt_um = sp.Matrix([sp.diff(u[0], t), sp.diff(u[1], t), sp.diff(u[2], t), sp.diff(zeta_m, t)])
V_p = sp.simplify(n_p0.dot(dt_up))
V_m = sp.simplify(n_m0.dot(dt_um))
print("\n=== LINEAR V_s = n0 · v_face ===")
print(f"V_plus  = {V_p}")
print(f"V_minus = {V_m}")
print(f"V_plus  coeff of d_t u_z = {sp.expand(V_p).coeff(sp.diff(uz(x, y, t), t))}")
print(f"V_minus coeff of d_t u_z = {sp.expand(V_m).coeff(sp.diff(uz(x, y, t), t))}")

# Traction t = -(dp + Lambda_X * A) n  => t · delta u = scalar * n·delta u
print("\n=== TRACTION PAIRING t · delta u uses n·delta u, same n_z=0 ===")
print(f"t_z background ~ n_z = {n_p0[2]} (plus), {n_m0[2]} (minus)")

# First-order n_z from zeta: dh/dz = 0 for z-independent zeta, so delta n_z = 0
print("\n=== FIRST-ORDER n_z FROM FACE PERTURBATION ===")
print(f"d h_plus / dz = {sp.diff(h_p, z)}")
print(f"d h_minus / dz = {sp.diff(h_m, z)}")
print(f"n_plus[2] exact = {n_p[2]}")
print(f"n_minus[2] exact = {n_m[2]}")

# --- virtual constraint linearisation ---
# Σ_mat = Σ_E(x(X,t),t) * J , J = det(I + ∇u) ≈ 1 + div u
# LAB_HELD: Σ_E(x) = rho(y) * W(y) * (1+theta) * (thickness factor)
# At linear virtual order, delta_v Σ / Σ0 = delta_v theta + delta_v e_W + div delta_v u
#   + delta_v u · ∇ log Σ0
# Σ0 = Σ0(y) => ∇ log Σ0 = (0, d_y log Σ0, 0)
rho0, W0s = sp.symbols("rho0 W0", positive=True)
# representative-independent: Σ0 = Σ0(y)
lnS = sp.Function("logSigma0")
dv_ux, dv_uy, dv_uz = sp.symbols("dv_ux dv_uy dv_uz", cls=sp.Function)
dv_th, dv_e = sp.symbols("dv_theta dv_eW")
div_dv = sp.diff(dv_ux(x, y), x) + sp.diff(dv_uy(x, y), y) + sp.diff(dv_uz(x, y), z)
# z-independent virtual fields: dv_uz = dv_uz(x,y), d_z = 0
div_dv_P = sp.diff(dv_ux(x, y), x) + sp.diff(dv_uy(x, y), y)
adv_P = dv_uy(x, y) * sp.diff(lnS(y), y)  # only y component
constraint_P = dv_th + dv_e + div_dv_P + adv_P
print("\n=== VIRTUAL CONSTRAINT on class P (LAB_HELD, Σ0=Σ0(y)) ===")
print(f"constraint = {constraint_P}")
print(f"has dv_uz = {constraint_P.has(dv_uz)}")
print(f"div_P has d_z dv_uz = {div_dv_P.has(dv_uz)}")

# MATERIAL_ADVECTED: Q_bg(chi) = Q(y - u_y) to linear order if Q=Q(y)
print("\n=== MATERIAL_ADVECTED pullback Q(y - u_y) ===")
print("linearised Q(chi) - Q(y) = -u_y * d_y Q")
print("COEFF_u_z in material pullback of Q(y) = 0")
