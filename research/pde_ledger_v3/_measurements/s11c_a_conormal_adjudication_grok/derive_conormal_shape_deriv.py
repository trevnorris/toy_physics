#!/usr/bin/env python3
"""Independent FROM-SPEC derivation of S11c-a T-a' CONORMAL_DERIV.

Object: face operator n̂_s · ∇₄, INCLUDING evaluation on the supplied graph,
first ε-order shape derivative for ONE case:
  LAB_HELD | FACE_MINUS (s=-1) | DOF δW

Authority (read, not engines):
  research/pde_ledger_v3/directives/S11c_a_SHARED_PHYSICS.md
  §2a  W_bg, η, σ_W, w1
  §3a  R_s, h_s, F_s, n̂_s, a_s
  §3c  δ[f(x,h_s)] = δf(x,h_s⁰) + δh_s ∂_w f⁰(x,h_s⁰); evaluate at h_s⁰ = s W_bg/2
  §4   T-a′ · Conormal ⇒ S11CA_CONORMAL_DERIV

Method:
  Build n̂ from F=w−h; form C[f] = n̂·(∇₄f)|_{w=h}; shape-differentiate with
  bulk field held Eulerian-fixed (probe). Print operands and residuals.
"""
from __future__ import annotations

import sympy as sp

# ---------------------------------------------------------------------------
# Coordinates / bookkeepers
# ---------------------------------------------------------------------------
eps = sp.symbols("epsilon", real=True, positive=True)
eta, sigma_W = sp.symbols("eta_bg sigma_W", real=True)
W0 = sp.symbols("W_0", real=True, positive=True)

x1, x2, x3, w = sp.symbols("x1 x2 x3 w", real=True)
# Profile and its in-plane first jets (§2a). Treat as independent symbols
# at the evaluation point (local Taylor data).
w1 = sp.symbols("w1_profile", real=True)
w1_d1, w1_d2, w1_d3 = sp.symbols(
    "w1_profile_d1 w1_profile_d2 w1_profile_d3", real=True
)

# Wave DOF for this case: δW only (ζ_c held at 0). e_W ≡ δW/W_0 (§1a).
e_W = sp.symbols("e_W", real=True)
e_W_d1, e_W_d2, e_W_d3 = sp.symbols("e_W_d1 e_W_d2 e_W_d3", real=True)

# Face sign for FACE_MINUS
s = sp.Integer(-1)

# ---------------------------------------------------------------------------
# §2a background thickness and its in-plane gradient
#   W_bg = W_0 [1 + η w1(ξ)]
#   ∂_{y_i} W_bg = σ_W ∂_{ξ_i} w1
# ---------------------------------------------------------------------------
W_bg = W0 * (1 + eta * w1)
dW_bg = (sigma_W * w1_d1, sigma_W * w1_d2, sigma_W * w1_d3)

# ---------------------------------------------------------------------------
# §3a LAB_HELD Eulerian graph height
#   h_s^L = s W_bg(x)/2 + ζ_s ,   ζ_s = ζ_c + s δW/2
# DOF δW only ⇒ ζ_c = 0, δW = W_0 e_W (amplitude tracked by eps):
#   h = s W_bg/2 + eps * (s W_0 e_W / 2)
# ---------------------------------------------------------------------------
dh_DOF = s * W0 * e_W / 2  # first ε variation of height from δW
d_grad_h_DOF = (
    s * W0 * e_W_d1 / 2,
    s * W0 * e_W_d2 / 2,
    s * W0 * e_W_d3 / 2,
)

h = s * W_bg / 2 + eps * dh_DOF
grad_h = tuple(
    s * dW_bg[i] / 2 + eps * d_grad_h_DOF[i] for i in range(3)
)

# Graph measure a = √(1+|∇_x h|²)
a2 = 1 + sum(g**2 for g in grad_h)
a = sp.sqrt(a2)

# ---------------------------------------------------------------------------
# §3a outward unit normal to F=w−h=0 with orientation s(n̂·ŵ)>0
#   ∇₄ F = (−∇_x h, 1), |∇₄ F| = a
#   n̂ = s * ∇₄F / |∇₄F|   satisfies s(n̂·ŵ) = s*(s/a) = 1/a > 0
# ---------------------------------------------------------------------------
n_vec = (
    s * (-grad_h[0]) / a,
    s * (-grad_h[1]) / a,
    s * (-grad_h[2]) / a,
    s * (sp.Integer(1)) / a,
)

# ---------------------------------------------------------------------------
# Generic bulk probe: traced gradient g = (∇₄ f)|_face and its normal jet
# Symbols stand for values AT THE BACKGROUND FACE h⁰ = s W_bg/2 (§3c).
# ---------------------------------------------------------------------------
g1, g2, g3, g4 = sp.symbols(
    "trace_grad_f_1 trace_grad_f_2 trace_grad_f_3 trace_grad_f_4", real=True
)
dw_g1, dw_g2, dw_g3, dw_g4 = sp.symbols(
    "trace_grad_f_1_dw trace_grad_f_2_dw trace_grad_f_3_dw trace_grad_f_4_dw",
    real=True,
)
g = (g1, g2, g3, g4)
dw_g = (dw_g1, dw_g2, dw_g3, dw_g4)

# ---------------------------------------------------------------------------
# Face operator evaluated on the graph (§4 T-a′, §3c):
#   C = n̂(h) · (∇₄ f)(x, h)
# With Eulerian-fixed probe, (∇₄ f)(x,h) = g_at_h0 + (h−h0) dw_g + O((h−h0)²)
# and h−h0 = eps*dh_DOF, so to O(eps):
#   (∇₄ f)|_h = g + eps*dh_DOF*dw_g
# ---------------------------------------------------------------------------
n_dot_g_at_h = sum(
    n_vec[i] * (g[i] + eps * dh_DOF * dw_g[i]) for i in range(4)
)

# Background value (eps→0)
C0 = sp.simplify(n_dot_g_at_h.subs(eps, 0))

# First ε-order shape derivative: ∂_eps of C at eps=0
# (exact quotient (C−C0)/eps, then eps→0)
C_shape = sp.simplify(
    sp.series(n_dot_g_at_h, eps, 0, 2).removeO().coeff(eps)
)

# ---------------------------------------------------------------------------
# Multigrade in (η, σ_W): keep first shape order in each background bookkeeper.
# Spec §2a: η and σ_W are independent grades; requested truncation is first
# shape order in each. Expand jointly to total degree ≤1 in {η, σ_W} for the
# LEADING flat-background identification, AND also emit the untruncated
# exact-in-background form (still O(eps^1)).
# ---------------------------------------------------------------------------
# Exact-in-background (retain full W_bg / σ_W dependence in n̂⁰, a⁰):
delta_C_exact_bg = sp.simplify(sp.expand(C_shape))

# First joint shape order in (η, σ_W): series in eta then sigma_W to O(1) each,
# dropping products η·σ_W and higher.
delta_C_O1 = sp.series(delta_C_exact_bg, eta, 0, 2).removeO()
delta_C_O1 = sp.series(delta_C_O1, sigma_W, 0, 2).removeO()
# Drop the η·σ_W cross term (second shape order)
delta_C_O1 = sp.expand(delta_C_O1)
cross = [
    t
    for t in sp.Add.make_args(delta_C_O1)
    if t.has(eta) and t.has(sigma_W)
]
delta_C_first_shape = sp.expand(delta_C_O1 - sum(cross))

# Flat-background (η=0, σ_W=0) slice — leading W_0 powers:
delta_C_flat = sp.simplify(
    delta_C_exact_bg.subs({eta: 0, sigma_W: 0})
)

# ---------------------------------------------------------------------------
# Analytic reconstruction from the product / chain rule (operands):
#   δC = δn̂ · g + n̂⁰ · (δh ∂_w g)
# with δn̂ = d/dε|_{0} n̂, n̂⁰ = n̂|_{ε=0}, δh = dh_DOF
# ---------------------------------------------------------------------------
n0 = tuple(sp.simplify(ni.subs(eps, 0)) for ni in n_vec)
dn = tuple(
    sp.simplify(sp.series(ni, eps, 0, 2).removeO().coeff(eps)) for ni in n_vec
)
delta_C_reconstructed = sp.simplify(
    sum(dn[i] * g[i] for i in range(4))
    + dh_DOF * sum(n0[i] * dw_g[i] for i in range(4))
)
recon_residual = sp.simplify(
    sp.expand(delta_C_exact_bg - delta_C_reconstructed)
)

# ---------------------------------------------------------------------------
# Print — operands, not conclusions
# ---------------------------------------------------------------------------
def dump(label: str, expr) -> None:
    print(f"{label} =")
    print(f"  {sp.factor(sp.together(expr))}")
    print(f"  expanded: {sp.expand(expr)}")
    print()


print("=" * 72)
print("S11c-a T-a′ CONORMAL_DERIV — independent FROM-SPEC derivation")
print("CASE: LAB_HELD | FACE_MINUS (s=-1) | DOF δW")
print("=" * 72)
print()
print("--- supplied inputs (§2a, §3a) ---")
print(f"s            = {s}")
print(f"W_bg         = {W_bg}")
print(f"dW_bg        = {dW_bg}")
print(f"h            = {h}")
print(f"grad_h       = {grad_h}")
print(f"dh_DOF (δh)  = {dh_DOF}")
print(f"d_grad_h_DOF = {d_grad_h_DOF}")
print()
print("--- outward unit normal n̂ = s ∇₄F / |∇₄F| ---")
print(f"a            = {a}")
for i, ni in enumerate(n_vec, 1):
    print(f"n̂_{i}         = {sp.simplify(ni)}")
print(f"n̂·ŵ check    = {sp.simplify(n_vec[3])} ;  s(n̂·ŵ) = {sp.simplify(s*n_vec[3])}")
print()
print("--- background face values (eps→0) ---")
print(f"h⁰           = {sp.simplify(h.subs(eps,0))}")
print(f"a⁰           = {sp.simplify(a.subs(eps,0))}")
for i, ni in enumerate(n0, 1):
    print(f"n̂⁰_{i}        = {ni}")
print(f"C⁰ = n̂⁰·g    = {C0}")
print()
print("--- first ε-order δn̂ components ---")
for i, dni in enumerate(dn, 1):
    print(f"δn̂_{i}        = {dni}")
print()
print("--- OBJECT: δ(n̂·∇₄ f)|_graph  to first ε ---")
dump("δC exact-in-background (full η,σ_W in n̂⁰)", delta_C_exact_bg)
dump("δC first joint shape order in (η,σ_W) [drop η·σ_W]", delta_C_first_shape)
dump("δC flat-background slice (η=σ_W=0)", delta_C_flat)
dump("δC reconstructed = δn̂·g + δh (n̂⁰·∂_w g)", delta_C_reconstructed)
print(f"RESIDUAL (exact_bg − reconstructed) = {recon_residual}")
print()

# Flat-slice term census (natural variables)
print("--- flat-background term census (natural variables) ---")
flat_exp = sp.expand(delta_C_flat)
for term in sp.Add.make_args(flat_exp):
    print(f"  {term}")
print()

# Orientation / measure factors at flat background
print("--- flat-background closed form ---")
print("At η=σ_W=0: a⁰=1, n̂⁰ = (0,0,0,s) = (0,0,0,-1) for FACE_MINUS")
print("δh = s W_0 e_W / 2 = -W_0 e_W / 2")
print("δ(∇h) = s W_0 ∇e_W / 2")
print("δa = 0 at flat (∇h⁰=0)")
print("δn̂_i (i=1,2,3) = −s δ(∂_i h) = −s*(s W_0 e_W_di/2) = −W_0 e_W_di/2")
print("δn̂_w = 0 at flat")
print("δC_flat = Σ_{i=1..3} δn̂_i g_i + δh n̂⁰_w ∂_w g_w")
print("         = Σ_{i=1..3} (−W_0 e_W_di/2) g_i + (−W_0 e_W/2)*(−1)*∂_w g_4")
print("         = −(W_0/2) Σ_i e_W_di g_i + (W_0/2) e_W ∂_w g_4")
pred_flat = (
    -(W0 / 2) * (e_W_d1 * g1 + e_W_d2 * g2 + e_W_d3 * g3)
    + (W0 / 2) * e_W * dw_g4
)
print(f"predicted flat = {pred_flat}")
print(f"RESIDUAL (δC_flat − predicted) = {sp.simplify(sp.expand(delta_C_flat - pred_flat))}")
print()

# First-shape σ_W terms at η=0 (tilted background from σ_W only)
delta_C_sigma = sp.series(
    delta_C_exact_bg.subs(eta, 0), sigma_W, 0, 2
).removeO()
delta_C_sigma = sp.expand(delta_C_sigma)
print("--- O(σ_W^1) slice at η=0 (exact − flat) ---")
sigma_only = sp.expand(delta_C_sigma - delta_C_flat.subs(sigma_W, 0))
# delta_C_flat already has σ_W=0; delta_C_sigma includes flat + O(σ_W)
sigma_only = sp.expand(delta_C_sigma - pred_flat)
dump("δC |_{η=0} to O(σ_W^1) minus flat piece", sigma_only)

# O(η^1) at σ_W=0: W_bg = W_0(1+η w1) enters h⁰ and δh? 
# δh = s W_0 e_W/2 is independent of η (δW = W_0 e_W by definition).
# At σ_W=0, ∇W_bg=0 so a⁰=1, n̂⁰=(0,0,0,s) still — η alone (zero jet) does not
# tilt the face. Expected: δC|_{σ_W=0} = δC_flat with no η.
delta_C_eta = sp.simplify(delta_C_exact_bg.subs(sigma_W, 0))
print("--- σ_W=0 slice (η retained): expect identical to flat ---")
dump("δC|_{σ_W=0}", delta_C_eta)
print(
    f"RESIDUAL (δC|σ_W=0 − flat) = {sp.simplify(sp.expand(delta_C_eta - delta_C_flat))}"
)
print()

print("=" * 72)
print("END independent derivation (engines not yet consulted)")
print("=" * 72)

