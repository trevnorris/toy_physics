#!/usr/bin/env python3
"""Identify Engine A (trace_grad_f) vs Engine B (conormalBackground) forms.

Authority: independent FROM-SPEC derivation in derive_conormal_shape_deriv.py
plus §3c Taylor of jets about the flat face.

Case: LAB_HELD | FACE_MINUS (s=-1) | DOF δW.

Spec-grounded identification tested here:
  (I1)  Engine A probe is Eulerian-fixed:  φ = φ0,  φ1 ≡ 0.
        Engine B computes δ_ε[ n̂·∇₄(φ0+ε φ1)|_h ]; geometric kernel is φ1→0.
  (I2)  trace_grad_f_i |_{h⁰}  =  ∂_i φ0 |_{h⁰}
        trace_grad_f_i_dw|_{h⁰} =  ∂_w ∂_i φ0 |_{h⁰}
        with φ0 ≡ conormalBackground.
  (I3)  §3c evaluation at h⁰ = s W_bg/2 = flat + (h⁰−flat), flat = s W_0/2,
        h⁰−flat = s W_0 η w1 / 2.  First η-order Taylor of (I2) about flat
        converts Engine A's h⁰-jet symbols into Engine B's flat-jet + η·W_0
        higher-normal-order terms (the extra W_0 power).

Print operands and residuals. No prose conclusions in the script output
beyond labelled residuals.
"""
from __future__ import annotations

import sympy as sp

# --- bookkeepers / profile (§2a) ---
W0, eta, sigma_W = sp.symbols("W_0 eta_bg sigma_W", real=True)
w1 = sp.symbols("w1_profile", real=True)
w1_d1, w1_d2, w1_d3 = sp.symbols(
    "w1_profile_d1 w1_profile_d2 w1_profile_d3", real=True
)
e_W = sp.symbols("e_W", real=True)
e_W_d1, e_W_d2, e_W_d3 = sp.symbols("e_W_d1 e_W_d2 e_W_d3", real=True)

# δW = W_0 e_W  (Engine B: deltaWidth)
delta_W = W0 * e_W
d_delta = (W0 * e_W_d1, W0 * e_W_d2, W0 * e_W_d3)
w1_d = (w1_d1, w1_d2, w1_d3)

# --- Engine A jets AT h⁰ (symbols) ---
g1, g2, g3, g4 = sp.symbols(
    "trace_grad_f_1 trace_grad_f_2 trace_grad_f_3 trace_grad_f_4", real=True
)
dw_g1, dw_g2, dw_g3, dw_g4 = sp.symbols(
    "d_w_trace_grad_f_1 d_w_trace_grad_f_2 d_w_trace_grad_f_3 d_w_trace_grad_f_4",
    real=True,
)
# Engine A finalized form (from engine extract; matches independent derivation
# flat + O(σ_W^1), η absorbed into h⁰-jet meanings):
engine_A = (
    W0
    * (
        2 * dw_g4 * e_W
        - 2 * e_W_d1 * g1
        - 2 * e_W_d2 * g2
        - 2 * e_W_d3 * g3
        + sigma_W
        * (
            dw_g1 * e_W * w1_d1
            + dw_g2 * e_W * w1_d2
            + dw_g3 * e_W * w1_d3
            + g4 * (e_W_d1 * w1_d1 + e_W_d2 * w1_d2 + e_W_d3 * w1_d3)
        )
    )
    / 4
)

# --- Engine B flat-face jets of φ0 = conormalBackground ---
# Derivative orders at w = flat = -W_0/2:
#   phi_di     = ∂_{x_i} φ0
#   phi_dw     = ∂_w φ0          (= g_4 at flat)
#   phi_di_dw  = ∂_{x_i} ∂_w φ0
#   phi_dww    = ∂_{ww} φ0
#   phi_di_dww = ∂_{x_i} ∂_{ww} φ0
#   phi_dwww   = ∂_{www} φ0
phi_d1, phi_d2, phi_d3 = sp.symbols("phi0_d1 phi0_d2 phi0_d3", real=True)
phi_dw = sp.symbols("phi0_dw", real=True)
phi_d1_dw, phi_d2_dw, phi_d3_dw = sp.symbols(
    "phi0_d1_dw phi0_d2_dw phi0_d3_dw", real=True
)
phi_dww = sp.symbols("phi0_dww", real=True)
phi_d1_dww, phi_d2_dww, phi_d3_dww = sp.symbols(
    "phi0_d1_dww phi0_d2_dww phi0_d3_dww", real=True
)
phi_dwww = sp.symbols("phi0_dwww", real=True)

phi_di = (phi_d1, phi_d2, phi_d3)
phi_di_dw = (phi_d1_dw, phi_d2_dw, phi_d3_dw)
phi_di_dww = (phi_d1_dww, phi_d2_dww, phi_d3_dww)

# Engine B geometric kernel (φ1=0), transcribed from WL_S11CA_CONORMAL_DERIV
# LAB_HELD|FACE_MINUS|DOF_DELTA_W SHAPE_DERIVATIVE, with deltaWidth→W_0 e_W.
# Keep all multigrades present in the WL emit; compare after dropping η·σ_W.
engine_B_geom = sp.Integer(0)
# (deltaWidth * ∂_{ww} φ0)/2
engine_B_geom += delta_W * phi_dww / 2
# -(etaBg*W0*deltaWidth*w1*∂_{www} φ0)/4
engine_B_geom += -eta * W0 * delta_W * w1 * phi_dwww / 4
for i in range(3):
    # (sigmaW * w1_di * ∂_w φ0 * ∂_i deltaWidth)/4
    engine_B_geom += sigma_W * w1_d[i] * phi_dw * d_delta[i] / 4
    # -(eta*sigma*W0*w1*w1_di*∂_{ww} φ0 * ∂_i deltaWidth)/8
    engine_B_geom += (
        -eta * sigma_W * W0 * w1 * w1_d[i] * phi_dww * d_delta[i] / 8
    )
    # -(∂_i deltaWidth * ∂_i φ0)/2
    engine_B_geom += -d_delta[i] * phi_di[i] / 2
    # (sigmaW * deltaWidth * w1_di * ∂_i ∂_w φ0)/4
    engine_B_geom += sigma_W * delta_W * w1_d[i] * phi_di_dw[i] / 4
    # (eta*W0*w1*∂_i deltaWidth*∂_i ∂_w φ0)/4
    engine_B_geom += eta * W0 * w1 * d_delta[i] * phi_di_dw[i] / 4
    # -(eta*sigma*W0*deltaWidth*w1*w1_di*∂_i ∂_{ww} φ0)/8
    engine_B_geom += (
        -eta * sigma_W * W0 * delta_W * w1 * w1_d[i] * phi_di_dww[i] / 8
    )

engine_B_geom = sp.expand(engine_B_geom)

# Engine B φ1 = conormalPerturbation pieces (background operator on wave field)
# From WL emit FACE_MINUS:
#   -∂_w φ1|_flat
#   +(eta*W0*w1*∂_{ww} φ1)/2
#   -(sigmaW*w1_di*∂_i φ1)/2
#   +(eta*sigma*W0*w1*w1_di*∂_i ∂_w φ1)/4
phi1_dw = sp.symbols("phi1_dw", real=True)
phi1_dww = sp.symbols("phi1_dww", real=True)
phi1_d1, phi1_d2, phi1_d3 = sp.symbols("phi1_d1 phi1_d2 phi1_d3", real=True)
phi1_d1_dw, phi1_d2_dw, phi1_d3_dw = sp.symbols(
    "phi1_d1_dw phi1_d2_dw phi1_d3_dw", real=True
)
phi1_di = (phi1_d1, phi1_d2, phi1_d3)
phi1_di_dw = (phi1_d1_dw, phi1_d2_dw, phi1_d3_dw)

engine_B_phi1 = -phi1_dw
engine_B_phi1 += eta * W0 * w1 * phi1_dww / 2
for i in range(3):
    engine_B_phi1 += -sigma_W * w1_d[i] * phi1_di[i] / 2
    engine_B_phi1 += eta * sigma_W * W0 * w1 * w1_d[i] * phi1_di_dw[i] / 4
engine_B_phi1 = sp.expand(engine_B_phi1)

# ---------------------------------------------------------------------------
# (I3) Taylor-expand Engine A h⁰-jets about flat to O(η^1), drop η·σ_W
# h⁰ − flat = s W_0 η w1 / 2 with s=-1 ⇒ −W_0 η w1 / 2
# ---------------------------------------------------------------------------
h0_minus_flat = -W0 * eta * w1 / 2

# g_i(h⁰) = phi_di + (h⁰−flat)*phi_di_dw + O(η²)
# dw_g_i(h⁰) = phi_di_dw + (h⁰−flat)*phi_di_dww + O(η²)
# g_4(h⁰) = phi_dw + (h⁰−flat)*phi_dww + O(η²)
# dw_g_4(h⁰) = phi_dww + (h⁰−flat)*phi_dwww + O(η²)
subs_A_jets = {
    g1: phi_d1 + h0_minus_flat * phi_d1_dw,
    g2: phi_d2 + h0_minus_flat * phi_d2_dw,
    g3: phi_d3 + h0_minus_flat * phi_d3_dw,
    g4: phi_dw + h0_minus_flat * phi_dww,
    dw_g1: phi_d1_dw + h0_minus_flat * phi_d1_dww,
    dw_g2: phi_d2_dw + h0_minus_flat * phi_d2_dww,
    dw_g3: phi_d3_dw + h0_minus_flat * phi_d3_dww,
    dw_g4: phi_dww + h0_minus_flat * phi_dwww,
}

engine_A_expanded = sp.expand(engine_A.subs(subs_A_jets))
# Drop O(η²): series in eta to order 2 (keeps η^0 and η^1)
engine_A_O_eta = sp.series(engine_A_expanded, eta, 0, 2).removeO()
engine_A_O_eta = sp.expand(engine_A_O_eta)

# Drop η·σ_W cross (second joint shape order) — Engine A finalize has no η·σ_W
def drop_eta_sigma(expr):
    parts = sp.Add.make_args(sp.expand(expr))
    kept = [t for t in parts if not (t.has(eta) and t.has(sigma_W))]
    return sp.expand(sum(kept))


engine_A_cmp = drop_eta_sigma(engine_A_O_eta)
engine_B_geom_cmp = drop_eta_sigma(engine_B_geom)

residual_geom = sp.simplify(sp.expand(engine_A_cmp - engine_B_geom_cmp))

# Also compare flat+σ_W slice (η=0) without needing Taylor
engine_A_eta0 = sp.expand(engine_A.subs(eta, 0))  # eta not in A anyway
engine_B_eta0 = sp.expand(engine_B_geom.subs(eta, 0))
# At η=0, identification is g_i = phi_di, dw_g_i = phi_di_dw, etc.
subs_flat = {
    g1: phi_d1,
    g2: phi_d2,
    g3: phi_d3,
    g4: phi_dw,
    dw_g1: phi_d1_dw,
    dw_g2: phi_d2_dw,
    dw_g3: phi_d3_dw,
    dw_g4: phi_dww,
}
residual_flat = sp.simplify(
    sp.expand(engine_A.subs(subs_flat) - engine_B_eta0)
)

# φ1 residual: Engine B full − Engine B geom
# (= background conormal on φ1, Taylor-expanded about flat)
print("=" * 72)
print("IDENTIFICATION TEST — Engine A vs Engine B CONORMAL_DERIV")
print("CASE: LAB_HELD | FACE_MINUS | DOF δW")
print("=" * 72)
print()
print("--- (I1)(I2)(I3) identification ---")
print("I1: geometric kernel ⇔ set conormalPerturbation φ1 = 0")
print("I2: trace_grad_f_i |_{h0} = ∂_i φ0 |_{h0} = ∂_i conormalBackground |_{h0}")
print("    d_w_trace_grad_f_i |_{h0} = ∂_w ∂_i φ0 |_{h0}")
print("I3: h0 − flat = s W_0 η w1 / 2 = -W_0 η w1 / 2  (s=-1)")
print("    Taylor of I2 jets about flat brings η·W_0·(higher normal jets)")
print()
print(f"h0_minus_flat = {h0_minus_flat}")
print()
print("--- Engine A (h0-jet symbols, finalized flat+O(σ_W)) ---")
print(engine_A)
print()
print("--- Engine A after I3 Taylor O(η^1), drop η·σ_W ---")
print(engine_A_cmp)
print()
print("--- Engine B geometric (φ1=0), drop η·σ_W ---")
print(engine_B_geom_cmp)
print()
print(f"RESIDUAL (A_Taylor − B_geom) after drop η·σ_W = {residual_geom}")
print()
print("--- η=0 slice under flat identification g=∇φ0|_flat ---")
print(f"RESIDUAL (A_flat_id - B|eta=0) = {residual_flat}")
print()
print("--- Engine B φ1 sector (background n̂·∇ acting on conormalPerturbation) ---")
print(engine_B_phi1)
print()
print("Engine A has no φ1 sector (Eulerian-fixed probe).")
print(f"RESIDUAL (B_full_geom+φ1 − B_geom) = {sp.expand(engine_B_phi1)}")
print(f"RESIDUAL (B_φ1 − 0) under I1 (φ1→0) = {sp.simplify(engine_B_phi1.subs({phi1_dw:0, phi1_dww:0, phi1_d1:0, phi1_d2:0, phi1_d3:0, phi1_d1_dw:0, phi1_d2_dw:0, phi1_d3_dw:0}))}")
print()

# Term-by-term census of matching at η=0
print("--- η=0 term census (A via flat id vs B) ---")
A_f = sp.expand(engine_A.subs(subs_flat))
B_f = engine_B_eta0
print("A terms:")
for t in sp.Add.make_args(A_f):
    print(f"  {t}")
print("B terms:")
for t in sp.Add.make_args(B_f):
    print(f"  {t}")
print()

# Explicit W_0 vs W_0^2 accounting for one η term
print("--- extra W_0 power from I3 (one sample monomial) ---")
# From A: -W0/2 * e_W_d1 * g1, with g1 = phi_d1 + (-W0 η w1/2) phi_d1_dw
# η piece: -W0/2 * e_W_d1 * (-W0 η w1/2) phi_d1_dw = (W0^2 η w1 e_W_d1 / 4) phi_d1_dw
sample_A_eta = sp.expand(
    (-W0 / 2 * e_W_d1 * (h0_minus_flat * phi_d1_dw))
)
sample_B_eta = sp.expand(eta * W0 * w1 * d_delta[0] * phi_d1_dw / 4)
print(f"A η-piece from δn̂_1 * (h0−flat)∂_w∂_1 φ0 = {sample_A_eta}")
print(f"B term (eta*W0*w1*∂_1 δW*∂_1∂_w φ0)/4   = {sample_B_eta}")
print(f"RESIDUAL sample = {sp.simplify(sample_A_eta - sample_B_eta)}")
print()
print("=" * 72)
print("END identification test")
print("=" * 72)

