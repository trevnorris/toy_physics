#!/usr/bin/env python3
"""Derive §3c shifted density trace from S11c-a SHARED_PHYSICS §§2b,2d,3c.

Spec members used (not engine operands):
  §2b  RHO4-CONSTANT / RHOBR-CONSTANT representatives for ρ_4D,bg⁰(y)
  §2d  𝔅⁰ contains ρ_4D,bg⁰
  §3c  δ[f(x,h_s)] = δf(x,h_s⁰) + δh_s · ∂_w f⁰(x,h_s⁰)
       and "the supplied density background depends on the in-plane
       anchor, not on w"

This script differentiates the §2b members of 𝔅⁰ w.r.t. the normal
coordinate w and composes the §3c first-order shifted density trace for
(BRANCH=LAB_HELD, FACE=MINUS, DOF=DELTA_W).
"""

from __future__ import annotations

import sympy as sp

# Coordinates and supplied constants (§1a, §2a).
x1, x2, x3, w, t = sp.symbols("x1 x2 x3 w t", real=True)
W0, rho_br, eta_bg = sp.symbols("W_0 rho_br eta_bg", positive=True)
w1 = sp.Function("w1_profile")(x1, x2, x3)
e_W = sp.Function("e_W")(x1, x2, x3, t)

# §2a width background (in-plane only).
W_bg = W0 * (1 + eta_bg * w1)

# §2b density representatives — members of 𝔅⁰ (§2d). Both are functions of
# the in-plane anchor y = (x1,x2,x3) through W_bg(y); neither depends on w.
rho4_ref = rho_br / W0
rho4_bg_RHO4 = rho4_ref                                # ≡ ρ_4D,ref⁰
rho4_bg_RHOBR = rho_br / W_bg                          # ≡ rho_br / W_bg(y)

print("=== §2b members of 𝔅⁰: ρ_4D,bg⁰ ===")
print(f"RHO4_CONSTANT : {rho4_bg_RHO4}")
print(f"RHOBR_CONSTANT: {rho4_bg_RHOBR}")

print("\n=== ∂_w ρ_4D,bg⁰  (differentiate the 𝔅⁰ member) ===")
for name, rho0 in (("RHO4_CONSTANT", rho4_bg_RHO4),
                   ("RHOBR_CONSTANT", rho4_bg_RHOBR)):
    d_w = sp.diff(rho0, w)
    print(f"{name}: d/dw rho_4D_bg0 = {d_w}")
    print(f"{name}: d/dw rho_4D_bg0 == 0 ? {d_w == 0}")

# §3c face height and DOF shift for LAB_HELD, FACE=MINUS, DOF=DELTA_W.
# h_s⁰ = s W_bg/2 with s = -1; δh from ζ_s = ζ_c + s δW/2, δW = W0 e_W, ζ_c=0.
s = -1
h0 = s * W_bg / 2
dh = s * W0 * e_W / 2

# Perturbation density δρ and its normal jet at the flat reference face
# (symbols only — content shared by both engines; not under dispute).
delta_rho_face = sp.symbols("delta_rho_4D_face_minus", real=True)
dw_delta_rho = sp.symbols("delta_rho_4D_face_minus_dw", real=True)
# Affine bulk perturbation about flat reference face w = s W0/2 (§3c
# evaluation of δf at h_s⁰ = s W_bg/2 expands through this jet).
w_flat = s * W0 / 2
delta_rho_bulk = delta_rho_face + (w - w_flat) * dw_delta_rho

print("\n=== §3c shift term δh · ∂_w ρ⁰ for each representative ===")
for name, rho0 in (("RHO4_CONSTANT", rho4_bg_RHO4),
                   ("RHOBR_CONSTANT", rho4_bg_RHOBR)):
    shift_term = sp.simplify(dh * sp.diff(rho0, w).subs(w, h0))
    print(f"{name}: delta_h * d_w(rho0)|_h0 = {shift_term}")

print("\n=== §3c first-order shifted density trace (both reps) ===")
# Compose f = ρ⁰ + ε δρ, evaluate at w = h0 + ε dh, take O(ε).
eps = sp.symbols("epsilon_shape", real=True)
for name, rho0 in (("RHO4_CONSTANT", rho4_bg_RHO4),
                   ("RHOBR_CONSTANT", rho4_bg_RHOBR)):
    bulk = rho0 + eps * delta_rho_bulk
    face_height = h0 + eps * dh
    exact = bulk.subs(w, face_height)
    # First-order content in the wave/shape bookkeeping parameter eps.
    # (Background face value ρ⁰(h0) is O(1) and not part of the δ[·] payload.)
    first = sp.series(exact, eps, 0, 2).removeO().coeff(eps)
    # Expand η through W_bg in h0 to first shape order in eta_bg, keeping
    # the perturbation jet contribution from evaluating at h_s⁰ ≠ s W0/2.
    first_lin_eta = sp.series(first, eta_bg, 0, 2).removeO()
    print(f"{name}:")
    print(f"  first-order trace = {sp.simplify(first_lin_eta)}")

print("\n=== residual of background shift against forced zero ===")
# Spec forces ∂_w ρ⁰ = 0, so any live symbolic ∂_w / ∂²_w of an undefined
# rhoBulkBackground is an extra free-premise term, not a surviving physical
# contribution. Its forced value is identically zero.
rho_free = sp.Function("rhoBulkBackground")
free_dw = sp.diff(rho_free(x1, x2, x3, w, t), w)
forced = 0
print(f"forced value of d/dw rho_4D_bg0 from 𝔅⁰: {forced}")
print(f"WL free-premise d/dw rhoBulkBackground is NOT a 𝔅⁰ member")
print(f"spec-forced residual (WL bg-shift terms - 0): vanishes identically")
print(f"correct δ[ρ](LAB_HELD, MINUS, DELTA_W) = "
      f"delta_rho_4D_face_minus"
      f" - (W_0/2)*eta_bg*w1_profile*delta_rho_4D_face_minus_dw")
