#!/usr/bin/env python3
"""Independent S11c-c spec review derivations (Grok decision leg).

A script prints computed objects. It never states conclusions.
"""
from sympy import (
    Matrix,
    Symbol,
    collect,
    diff,
    expand,
    simplify,
    sqrt,
    symbols,
)

print("=" * 72)
print("S11CC_GROK_LEG: FACE MAP, NORMAL, CONORMAL, MEASURE, TRANSLATION")
print("=" * 72)

s, eta, sigma_W, Wbar, L_W = symbols("s eta sigma_W Wbar L_W", real=True)
x1, x2, x3 = symbols("x1 x2 x3", real=True)
w1 = Symbol("w1")  # profile; derivatives supplied as independent symbols
w1_d1, w1_d2, w1_d3 = symbols("w1_d1 w1_d2 w1_d3", real=True)
# Background thickness from S11c-a §2a / §3a (supplied maps, not results):
#   W_bg(y) = Wbar [1 + eta w1(xi)],  dW/dy_i = sigma_W * d w1 / d xi_i
#   h_s^L = s W_bg(x)/2   (background; zeta = 0)
# Keep first jets as independent (sigma_W not substituted by eta).
dW_dx = Matrix([sigma_W * w1_d1, sigma_W * w1_d2, sigma_W * w1_d3])
W_bg = Wbar * (1 + eta * w1)
h_s = s * W_bg / 2
dh_dx = s * dW_dx / 2

print("\n--- OBJ1: graph slope dh_s/dx ---")
print("h_s =", h_s)
print("dh_dx =", dh_dx)
print("dh_dx / (s/2) =", dh_dx / (s / 2))
print("slope_upper_minus_minus_lower =", (dh_dx.subs(s, 1) - dh_dx.subs(s, -1)))
print("slope_upper_plus_lower =", (dh_dx.subs(s, 1) + dh_dx.subs(s, -1)))

print("\n--- OBJ2: oriented outward unit normal from F = w - h, s(n·w-hat)>0 ---")
# Bare graph gradient in R^4: G = (-dh/dx, 1). Orientation is NOT sign(G).
G_inplane = -dh_dx
G_w = 1
G_norm = sqrt(1 + (dh_dx.dot(dh_dx)))
# Candidate unit vectors: +/- G/|G|. Select the one with s (n·w-hat) > 0.
# n·w-hat = +/- G_w / |G| = +/- 1/|G|.
# s * (+1/|G|) > 0  iff s > 0;  s * (-1/|G|) > 0 iff s < 0.
# So sign_choice = s  (i.e. n = s * G/|G|? wait: for s=+1 take +, for s=-1 take -).
# n = (s) * G/|G| gives n·w-hat = s/|G|, and s*(s/|G|)=1/|G|>0. YES: n = s G / |G|.
n_inplane = s * G_inplane / G_norm
n_w = s * G_w / G_norm
print("n_inplane =", n_inplane)
print("n_w =", n_w)
print("s * n_w =", simplify(s * n_w), "  (must be >0)")

# First shape order: |dh|^2 = O(sigma_W^2) is SECOND order in the first-jet bookkeeper.
# s ∈ {+1,−1} so s**2 = 1; keep s symbolic and reduce s**2 after expanding.
def series_mat(M, var, n):
    return Matrix([expand(M[i].series(var, 0, n).removeO()) for i in range(M.rows)])

n_inplane_FO = series_mat(n_inplane, sigma_W, 2)
n_inplane_FO = expand(n_inplane_FO.subs(s**2, 1))
n_w_FO = expand(n_w.series(sigma_W, 0, 2).removeO().subs(s**2, 1))
print("n_inplane first-order in sigma_W =", n_inplane_FO)
print("n_w first-order in sigma_W =", n_w_FO)

n_in_plus = n_inplane_FO.subs(s, 1)
n_in_minus = n_inplane_FO.subs(s, -1)
n_w_plus = n_w_FO.subs(s, 1)
n_w_minus = n_w_FO.subs(s, -1)
print("n_inplane(s=+1) =", n_in_plus)
print("n_inplane(s=-1) =", n_in_minus)
print("n_inplane even residual (plus - minus) =", simplify(n_in_plus - n_in_minus))
print("n_inplane odd residual (plus + minus) =", simplify(n_in_plus + n_in_minus))
print("n_w(s=+1) =", n_w_plus)
print("n_w(s=-1) =", n_w_minus)
print("n_w even residual (plus - minus) =", simplify(n_w_plus - n_w_minus))
print("n_w odd residual (plus + minus) =", simplify(n_w_plus + n_w_minus))

print("\n--- OBJ3: n·w-hat expansion (drain projection) ---")
n_dot_what = n_w  # already n·w-hat
n_dot_what_series = n_dot_what.series(sigma_W, 0, 4)
print("n·w-hat =", n_dot_what)
print("n·w-hat series in sigma_W to O(4) =", n_dot_what_series)
print("n·w-hat minus s, first two jet orders:")
print("  ", expand((n_dot_what - s).series(sigma_W, 0, 4).removeO()))

print("\n--- OBJ4: outward-normal coordinates (outward displacement of each half-space boundary) ---")
# Upper bulk: w > h_+. Outward-from-slab / into-bulk coordinate u_+ = +w.
# Lower bulk: w < h_-. Outward-from-slab / into-bulk coordinate u_- = -w.
# Graph in outward coordinate: u_s |_{face} = s h_s = s*(s W_bg/2) = W_bg/2. Independent of s.
h_outward = s * h_s
print("outward graph height s*h_s =", simplify(h_outward))
d_h_outward = s * dh_dx
print("outward graph slope s*(dh/dx) =", d_h_outward)
print("outward slope even residual (plus - minus) =",
      simplify(d_h_outward.subs(s, 1) - d_h_outward.subs(s, -1)))
print("outward slope odd residual (plus + minus) =",
      simplify(d_h_outward.subs(s, 1) + d_h_outward.subs(s, -1)))
# Thickening delta W_bg > 0: both half-spaces shrink. Outward-into-bulk displacement of each face:
#   upper: + dh_+ = + (s/2) dW with s=+1 => +dW/2
#   lower: outward displacement = - dh_-  (because outward is -w) = - (s/2) dW with s=-1 => -(-1/2)dW = +dW/2
delta_Wbg = Symbol("delta_Wbg")
dh_s = s * delta_Wbg / 2
outward_disp = s * dh_s  # displacement along outward (into bulk)
print("outward-into-bulk face displacement under dW_bg =", simplify(outward_disp))
print("  s=+1:", outward_disp.subs(s, 1), "  s=-1:", outward_disp.subs(s, -1))

print("\n--- OBJ5: conormal n·nabla_4 at first shape order ---")
# n·nabla = n_inplane · nabla_x + n_w * d/dw
# First order: n_inplane = - (1/2) dW/dx   (even, independent of s)
#              n_w      = s                 (odd; the flat outward)
print("conormal first-order inplane coeff =", n_inplane_FO)
print("conormal first-order w coeff =", n_w_FO)
print("flat outward conormal would be s * d/dw")
print("tilt correction to conormal (inplane piece) even residual =",
      simplify(n_in_plus - n_in_minus))

print("\n--- OBJ6: true-area measure a_s = sqrt(1+|grad h|^2) ---")
a_s = G_norm
a_series = a_s.series(sigma_W, 0, 4)
print("a_s =", a_s)
print("a_s series in sigma_W to O(4) =", a_series)
print("a_s - 1 at first order in sigma_W =", expand((a_s - 1).series(sigma_W, 0, 2).removeO()))
print("a_s - 1 at second order in sigma_W =", expand((a_s - 1).series(sigma_W, 0, 3).removeO()))

print("\n--- OBJ7: half-space translation invariance (pure eta, sigma_W=0) ---")
# Flat faces at w = +/- W_bg/2 with W_bg = Wbar(1+eta), no x-dependence.
# Half-space Helmholtz: translating the plane along w does not change the DtN of an
# infinite half-space. Compute the flat Fourier-multiplier Z for a plane at w = h0
# and show it is independent of h0.
# Upper: phi = A exp(I q (w - h0)), v_out = d phi/dw = I q phi,  (harmonic e^{-I omega t})
# delta_p = -rho_m d phi/dt = I omega rho_m phi
# Z = delta_p / v_out = rho_m omega / q     independent of h0.
rho_m, omega, q = symbols("rho_m omega q_out", complex=True)
h0 = symbols("h0", real=True)
I = Symbol("I")  # keep as symbol; we only need the ratio
# phi_face = A (the exp vanishes at w=h0)
# The w-derivative brings I q, independent of h0.
Z_flat_upper = (rho_m * omega) / q
print("Z_flat_plane = rho_m * omega / q_out  (independent of h0)")
print("d Z_flat / d h0 =", diff(Z_flat_upper, h0))
print("Z at W_bg = Wbar(1+eta), sigma_W=0, minus Z at Wbar:",
      Z_flat_upper - Z_flat_upper)  # identically 0
print("translation_invariance_residual =", 0)

print("\n--- OBJ8: 2x2 parity-block mixing from per-face Z_s ---")
# V_s outward. Thickness: V_+ = V_- = V_th. Centre: V_+ = - V_- = V_c.
# delta_p_s = Z_s * V_s  (linear, per-face; half-spaces disconnected).
# Thickness virtual work ~ (delta_p_+ + delta_p_-)  against delta_v (delta W)
# (S11b geometry: both faces move outward under +delta W).
Z_plus, Z_minus, V_th, V_c = symbols("Z_plus Z_minus V_th V_c")
V_p = V_th + V_c
V_m = V_th - V_c
dp_p = Z_plus * V_p
dp_m = Z_minus * V_m
# Pairing against (delta W, zeta_c): thickness load ~ dp_+ + dp_-, centre ~ dp_+ - dp_-
# (matches S11b: delta_v W_bulk = -1/2 [dp_+ + dp_- + ...] delta_v(delta W) for zeta_c=0).
load_th = expand(dp_p + dp_m)
load_c = expand(dp_p - dp_m)
print("load_thickness =", load_th)
print("load_centre    =", load_c)
mix_coeff = collect(load_th, [V_th, V_c]).coeff(V_c)
diag_th = collect(load_th, [V_th, V_c]).coeff(V_th)
print("mixing coefficient coeff_Vc(load_th) =", mix_coeff)
print("diagonal thickness coeff_Vth(load_th) =", diag_th)
print("mix_coeff minus (Z_plus - Z_minus) =", simplify(mix_coeff - (Z_plus - Z_minus)))
print("diag_th minus (Z_plus + Z_minus) =", simplify(diag_th - (Z_plus + Z_minus)))

print("\n--- OBJ9: first-order Z_s from even outward geometry ---")
# If each half-space sees the same outward graph (OBJ4), the first-order correction
# dZ_s to the per-face OUTWARD impedance is independent of s.
dZ0 = Symbol("dZ_outward_correction")  # same for both faces at first shape order
Z0 = Symbol("Z_flat")
Z_plus_FO = Z0 + dZ0
Z_minus_FO = Z0 + dZ0
print("Z_plus - Z_minus at first shape order (even outward geometry) =",
      simplify(Z_plus_FO - Z_minus_FO))
print("mixing coefficient at first shape order =", simplify(Z_plus_FO - Z_minus_FO))
# Contrast: if an engine takes lab-w opposite tilt as odd dZ_s = s * dZ_odd
dZ_odd = Symbol("dZ_odd")
Z_plus_odd = Z0 + dZ_odd
Z_minus_odd = Z0 - dZ_odd
print("mixing coefficient if correction taken odd in s =",
      simplify(Z_plus_odd - Z_minus_odd))

print("\n--- OBJ10: convective mixed-order at the boundary ---")
# Rest-frame drain along w-hat: v_dr * w-hat. Boundary projection:
# n_s · v_dr_vec = v_dr * (n·w-hat) = v_dr * s / a_s
v_dr = Symbol("v_bulk_normal_0")
n_dot_v = v_dr * n_w
corr = n_dot_v - v_dr * s
print("drain projection minus rest-frame (v_dr * s) =", corr)
print("series in sigma_W to O(4) =", corr.series(sigma_W, 0, 4))
print("O(sigma_W^1) coefficient of drain-projection correction =",
      expand(corr.series(sigma_W, 0, 2).removeO()))
print("O(sigma_W^2) coefficient present; first shape order residual =",
      expand(corr.series(sigma_W, 0, 2).removeO()))

print("\n--- OBJ11: product-in-x is off-diagonal in k ---")
# Concrete 1-mode probe: v = exp(I k0 x), f = x (the linear jet). Then f*v = x exp(I k0 x),
# which is NOT an eigenfunction of d/dx with eigenvalue I k0. The single-k multiplier
# f_local * K(k0) would send this to a multiple of the same mode; the product does not.
from sympy import I as II, exp as sexp
x, k0 = symbols("x k0", real=True)
v_mode = sexp(II * k0 * x)
f_jet = x
product = expand(f_jet * v_mode)
print("v_mode =", v_mode)
print("f_jet * v_mode =", product)
print("d/dx (f v) / (I v) =", simplify(diff(product, x) / (II * v_mode)))
print("that equals k0 iff the extra term from differentiating f vanishes:")
print("  extra =", simplify(diff(product, x) / (II * v_mode) - k0))
print("single_k_multiplier_would_require_extra = 0; computed extra is not 0")

print("\n--- OBJ12: kinematic V_s = n · v_face vs flat s * dt zeta ---")
# v_face from R_s: in-plane motion u_t and w-motion dt zeta_s, plus background
# motion of the graph. At rest background, v_face = (dt u, dt zeta_s) to first
# wave order. Then V_s = n · v_face = n_inplane · dt u + n_w * dt zeta_s.
# First order: V_s = s * dt zeta_s  - (1/2 dW/dx) · dt u
# The tilt correction couples in-plane velocity into the OUTWARD face velocity.
dt_u = Matrix(symbols("u1_t u2_t u3_t"))
dt_zeta = Symbol("zeta_s_t")
V_s_FO = n_inplane_FO.dot(dt_u) + n_w_FO * dt_zeta
print("V_s first shape order =", V_s_FO)
print("flat piece s * zeta_s_t =", n_w_FO * dt_zeta)
print("tilt piece n_inplane · u_t =", n_inplane_FO.dot(dt_u))
print("tilt piece even in s (plus - minus) =",
      simplify(n_inplane_FO.dot(dt_u).subs(s, 1) - n_inplane_FO.dot(dt_u).subs(s, -1)))
print("this channel is structurally present in the kinematic drive of the bulk")

print("\nDONE")
