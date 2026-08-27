#!/usr/bin/env python3
"""Independent §3c grounded density-trace derivation from S11c_a_SHARED_PHYSICS.

Sources (spec only — not engine operands):
  §2b  RHO4-CONSTANT / RHOBR-CONSTANT representatives for ρ_4D,bg⁰(y)
  §2d  𝔅⁰ contains ρ_4D,bg⁰ (per anchoring and representative)
  §3c  δ[f(x,h_s)] = δf(x,h_s⁰) + δh_s · ∂_w f⁰(x,h_s⁰)
       background face value / normal derivative from differentiating 𝔅⁰;
       supplied density background depends on the in-plane anchor, not on w.

Emits both representatives' background face values, the differentiated
normal derivatives (must *emerge*), and the §3c shift term δh·∂_w ρ⁰.
"""
from __future__ import annotations

import sympy as sp

x1, x2, x3, w, t = sp.symbols("x1 x2 x3 w t", real=True)
W0, rho_br, eta_bg = sp.symbols("W_0 rho_br eta_bg", positive=True)
w1 = sp.Function("w1_profile")(x1, x2, x3)
e_W = sp.Function("e_W")(x1, x2, x3, t)
zeta_c = sp.Function("zeta_c")(x1, x2, x3, t)

# §2a width background — in-plane only.
W_bg = W0 * (1 + eta_bg * w1)

# §2b members of 𝔅⁰ — functions of the in-plane anchor; neither depends on w.
rho4_ref = rho_br / W0
reps = {
    "RHO4_CONSTANT": rho4_ref,
    "RHOBR_CONSTANT": rho_br / W_bg,
}

print("=== §2b ρ_4D,bg⁰ representatives ===")
for name, rho0 in reps.items():
    print(f"{name}: {rho0}")

print("\n=== background normal derivative by differentiation (must emerge) ===")
for name, rho0 in reps.items():
    dw = sp.diff(rho0, w)
    print(f"{name}: D[rho4_bg, w] = {dw}")
    print(f"{name}: D[rho4_bg, w] == 0 ? {dw == 0}")

print("\n=== background face values at h_s⁰ = s W_bg/2 (LAB_HELD, s=-1) ===")
s = -1
h0 = s * W_bg / 2
for name, rho0 in reps.items():
    face = sp.simplify(rho0.subs(w, h0))
    print(f"{name}: rho4_bg(h0) = {face}")

print("\n=== §3c shift term δh · ∂_w ρ⁰ |_{h0}  (DOF=DELTA_W, ζ_c=0) ===")
# ζ_s = ζ_c + s δW/2; δW = W0 e_W; for DELTA_W take ζ_c=0.
dh = s * W0 * e_W / 2
for name, rho0 in reps.items():
    shift = sp.simplify(dh * sp.diff(rho0, w).subs(w, h0))
    print(f"{name}: delta_h * D_w(rho0)|_h0 = {shift}")

print("\n=== representative-independence of face value (not a hardcoded 0) ===")
face_r4 = reps["RHO4_CONSTANT"]
face_rb = reps["RHOBR_CONSTANT"]
print(f"RHO4 face - RHOBR face = {sp.simplify(face_r4 - face_rb)}")
print(f"faces_equal? {sp.simplify(face_r4 - face_rb) == 0}")

print("\n=== composed first-order δ[ρ] at LAB_HELD/MINUS/DELTA_W ===")
# f = ρ⁰ + ε δρ; evaluate at w = h0 + ε dh; take O(ε).
# With ∂_w ρ⁰ = 0 the shift vanishes; perturbation retains η through h0.
eps = sp.symbols("epsilon_shape", real=True)
delta_rho_face = sp.symbols("delta_rho_4D_face_minus", real=True)
dw_delta_rho = sp.symbols("delta_rho_4D_face_minus_dw", real=True)
w_flat = s * W0 / 2
delta_rho_bulk = delta_rho_face + (w - w_flat) * dw_delta_rho
for name, rho0 in reps.items():
    bulk = rho0 + eps * delta_rho_bulk
    exact = bulk.subs(w, h0 + eps * dh)
    first = sp.series(exact, eps, 0, 2).removeO().coeff(eps)
    first_lin_eta = sp.series(first, eta_bg, 0, 2).removeO()
    print(f"{name}: first-order delta[rho] = {sp.simplify(first_lin_eta)}")
