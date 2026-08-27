#!/usr/bin/env python3
"""
Q1 — Does the perturbation current's w-dependence change the first-shape-order
projection of ∂_t δρ + ∇₄·δj against the dynamic window (IBP in w per §1b)?

Choice ②: δj_a = δj_a(x,w,t) full bulk field inside the w-integral.
Choice ①: δj_a = δj_a^⊙(x,t) := value at background face, constant in w,
          with in-plane jets carried as w-independent fields.

Prints operands and the difference (② − ①). No conclusions in the script.
"""
from __future__ import annotations

import sympy as sp

x1, x2, x3, w, t = sp.symbols("x1 x2 x3 w t", real=True)
spatial = (x1, x2, x3)
bounds = (w, -sp.oo, sp.oo)

# Background faces and first-shape face rates (already evaluated at ε=0).
h_plus0 = sp.Function("h_plus0")(*spatial, t)
h_minus0 = sp.Function("h_minus0")(*spatial, t)
dh_plus = sp.Function("dh_plus")(*spatial, t)
dh_minus = sp.Function("dh_minus")(*spatial, t)

# Two-argument window on background args; shape variation via chain rule.
O = sp.Function("O")
g1, g2 = sp.symbols("g1 g2", real=True)
G_plus0 = w - h_plus0
G_minus0 = h_minus0 - w
Omega0 = O(G_plus0, G_minus0)

# Partial jets of O at background arguments (canonical Subs form).
O_1 = sp.Subs(sp.diff(O(g1, g2), g1), (g1, g2), (G_plus0, G_minus0))
O_2 = sp.Subs(sp.diff(O(g1, g2), g2), (g1, g2), (G_plus0, G_minus0))
O_11 = sp.Subs(sp.diff(O(g1, g2), g1, 2), (g1, g2), (G_plus0, G_minus0))
O_22 = sp.Subs(sp.diff(O(g1, g2), g2, 2), (g1, g2), (G_plus0, G_minus0))
O_12 = sp.Subs(sp.diff(O(g1, g2), g1, g2), (g1, g2), (G_plus0, G_minus0))

# ∂ε Ω |_0 and ∂ε(∂w Ω)|_0 by chain rule.
# G+ = w - (h+0 + ε dh+), G- = (h-0 + ε dh-) - w
# ⇒ ∂ε G+|_0 = -dh+, ∂ε G-|_0 = dh-
# ∂w G+ = 1, ∂w G- = -1
# ∂w Ω = O_1 - O_2
# ∂ε Ω = O_1*(-dh+) + O_2*(dh-)
# ∂ε(∂w Ω) = ∂ε(O_1 - O_2) = (O_11*(-dh+) + O_12*dh-) - (O_12*(-dh+) + O_22*dh-)
dOm0 = O_1 * (-dh_plus) + O_2 * (dh_minus)
d_eps_dw_Om0 = (O_11 * (-dh_plus) + O_12 * dh_minus) - (
    O_12 * (-dh_plus) + O_22 * dh_minus
)
# = -O_11 dh+ + O_12 dh- + O_12 dh+ - O_22 dh-
d_eps_dw_Om0 = sp.expand(d_eps_dw_Om0)

# Bulk perturbation fields (choice ②) and face representatives (choice ①).
drho = sp.Function("delta_rho")(*spatial, w, t)
dj = [sp.Function(f"delta_j_{i}")(*spatial, w, t) for i in (1, 2, 3)]
dj_w = sp.Function("delta_j_w")(*spatial, w, t)
dj_face = [sp.Function(f"delta_j_{i}_face")(*spatial, t) for i in (1, 2, 3)]
dj_w_face = sp.Function("delta_j_w_face")(*spatial, t)

div_x_bulk = sum(sp.diff(dj[i], spatial[i]) for i in range(3))
div_face = sum(sp.diff(dj_face[i], spatial[i]) for i in range(3))

# ---------------------------------------------------------------------------
# IBP in w of ∫ Ω (∂t δρ + ∇₄·δj) dw at fixed window:
#   ∫ Ω ∂w j_w = [Ω j_w]_{±∞} − ∫ (∂w Ω) j_w
# PROPERTY: Ω → 0 as w → ±∞  ⇒  boundary vanishes.
# ---------------------------------------------------------------------------
integrand_before_2 = Omega0 * (sp.diff(drho, t) + div_x_bulk + sp.diff(dj_w, w))
integrand_ibp_2 = Omega0 * (sp.diff(drho, t) + div_x_bulk) - sp.diff(Omega0, w) * dj_w
ibp_check_2 = sp.expand(integrand_before_2 - integrand_ibp_2 - sp.diff(Omega0 * dj_w, w))

integrand_before_1 = Omega0 * (sp.diff(drho, t) + div_face + sp.diff(dj_w_face, w))
integrand_ibp_1 = Omega0 * (sp.diff(drho, t) + div_face) - sp.diff(Omega0, w) * dj_w_face
ibp_check_1 = sp.expand(
    integrand_before_1 - integrand_ibp_1 - sp.diff(Omega0 * dj_w_face, w)
)

print("=== Q1 STAGE 0: IBP algebra identity (local densities) ===")
print("PROPERTY_USED: O(G+,G-) -> 0 as w -> +/- oo (compact/decaying window support)")
print("before_ibp_2 - ibp_2 - d/dw(Omega*delta_j_w) =", ibp_check_2)
print("before_ibp_1 - ibp_1 - d/dw(Omega*delta_j_w_face) =", ibp_check_1)
print("ibp_identity_holds_choice2:", ibp_check_2 == 0)
print("ibp_identity_holds_choice1:", ibp_check_1 == 0)
print()

# First-shape-order integrands: only the window is shape-differentiated.
# d/dε|_{0} of IBP form = (∂εΩ)(∂tδρ + div_x j) - (∂ε ∂w Ω) j_w
shape_int_2 = sp.expand(dOm0 * (sp.diff(drho, t) + div_x_bulk) - d_eps_dw_Om0 * dj_w)
shape_int_1 = sp.expand(dOm0 * (sp.diff(drho, t) + div_face) - d_eps_dw_Om0 * dj_w_face)
diff_int = sp.expand(shape_int_2 - shape_int_1)

print("=== Q1 STAGE 1: first-shape-order integrands (after IBP, boundary dropped) ===")
print("dOmega/deps|_0 =")
print(dOm0)
print()
print("d/deps(dOmega/dw)|_0 =")
print(d_eps_dw_Om0)
print()
print("shape_integrand_choice2 =")
print(shape_int_2)
print()
print("shape_integrand_choice1 =")
print(shape_int_1)
print()
print("difference_integrand (choice2 - choice1) =")
print(diff_int)
print()

diff_expected = sp.expand(
    dOm0 * (div_x_bulk - div_face) - d_eps_dw_Om0 * (dj_w - dj_w_face)
)
print("=== Q1 STAGE 2: structured difference ===")
print("expected_diff = (dOm0)*(div_x delta_j - div_face) - (d_eps_dw_Om0)*(delta_j_w - delta_j_w_face)")
print(diff_expected)
print()
print("diff_int - expected_diff =", sp.expand(diff_int - diff_expected))
print("matches_structured_form:", sp.expand(diff_int - diff_expected) == 0)
print("diff_has_delta_rho:", diff_int.has(drho))
print()

# Second IBP on the difference, using dOm0 → 0 at ±∞:
# −∫ (∂w δΩ) Δj_w = ∫ δΩ ∂w Δj_w − [δΩ Δj_w]_{±∞}
diff_after_second_ibp = sp.expand(dOm0 * (div_x_bulk - div_face + sp.diff(dj_w, w)))
second_ibp_local = sp.expand(
    diff_expected - diff_after_second_ibp + sp.diff(dOm0 * (dj_w - dj_w_face), w)
)
# Subs spellings of the same O-jet need .doit() before == 0 is meaningful.
second_ibp_local_doit = sp.expand(second_ibp_local.doit())

print("=== Q1 STAGE 3: optional second IBP on the difference (same decay of dOm0) ===")
print("diff_after_second_ibp = dOm0 * (div_4(delta_j) - div_face)")
print(diff_after_second_ibp)
print()
print("(expected_diff - diff_after_second_ibp) + d/dw[dOm0*(delta_j_w - delta_j_w_face)] =")
print(second_ibp_local)
print("second_ibp_local.doit() =", second_ibp_local_doit)
print("second_ibp_identity_holds (after doit):", second_ibp_local_doit == 0)
print()

# Probe: no w-dependence ⇒ difference vanishes.
probe = diff_int
for i in range(3):
    probe = probe.subs(dj[i], dj_face[i])
probe = sp.expand(probe.subs(dj_w, dj_w_face))

print("=== Q1 STAGE 4: probe — replace bulk current by face values ===")
print("diff |_{delta_j -> delta_j_face} =", probe)
print("vanishes_when_no_w_dependence:", probe == 0)
print()

# Concrete check: O = exp(-G+^2 - G-^2), δj_w = w, face value = 0 (at h0=0 slab center
# is not the face; take face value of δj_w at w=h_plus0 as symbol, and a linear-in-w field).
print("=== Q1 STAGE 5: concrete nonzero probe (decaying O, linear-in-w current) ===")
W0 = sp.symbols("W0", positive=True)
# Flat faces h±0 = ±W0/2, dh+ = A, dh- = B constants; δj_i=0, δj_w = c*w, face = c*h_plus0
A, B, c = sp.symbols("A B c", real=True)
Oc = sp.exp(-(G_plus0**2) - (G_minus0**2))
# Recompute dOm0, d_eps_dw for concrete O with h±0 = ±W0/2
Gp = w - W0 / 2
Gm = -W0 / 2 - w
Oc_flat = sp.exp(-(Gp**2) - (Gm**2))
# ∂ε Ω with h+ = W0/2+εA, h- = -W0/2+εB:
Oc_dyn = sp.exp(-((w - (W0 / 2 + sp.Symbol("eps") * A)) ** 2) - (((-W0 / 2 + sp.Symbol("eps") * B) - w) ** 2))
eps_s = sp.Symbol("eps")
dOc0 = sp.diff(Oc_dyn, eps_s).subs(eps_s, 0)
d_eps_dw_Oc0 = sp.diff(sp.diff(Oc_dyn, w), eps_s).subs(eps_s, 0)
# Choice ②: δj_w = c*w; choice ①: δj_w_face = c*(W0/2); div_x = 0 both
diff_conc = sp.expand(dOc0 * 0 - d_eps_dw_Oc0 * (c * w - c * W0 / 2))
# Integrate the difference
I_diff = sp.simplify(sp.integrate(diff_conc, (w, -sp.oo, sp.oo)))
print("concrete_diff_integrand =", diff_conc)
print("Integral(concrete_diff) =", I_diff)
print("concrete_diff_Integral_is_zero:", I_diff == 0)
print()

print("=== Q1 LITERAL RESIDUAL REPORT ===")
print("DIFF_INTEGRAND_choice2_minus_choice1 =")
print(diff_int)
print()
print("DIFF_IDENTICALLY_ZERO:", diff_int == 0)
print("DIFF_VANISHES_IFF_CURRENT_W_INDEPENDENT_EQUAL_FACE:", probe == 0)
print("CONCRETE_LINEAR_IN_W_CURRENT_INTEGRAL_DIFF_ZERO:", I_diff == 0)
print("SURVIVING_STRUCTURED_FORM:")
print("  Integral[ (dOmega/deps|_0)*(div_x delta_j(w) - div_x delta_j_face)")
print("           - (d/deps dOmega/dw|_0)*(delta_j_w(w) - delta_j_w_face) , {w,-oo,oo} ]")
print("EQUIVALENT_AFTER_SECOND_IBP_AND_DECAY (second_ibp_identity_holds=%s):" % (second_ibp_local_doit == 0))
print("  Integral[ (dOmega/deps|_0)*(div_4 delta_j(w) - div_x delta_j_face) , {w,-oo,oo} ]")
print("DEPENDS_ON: window shape variation (dh+/-, through dOm0), and w-variation of delta_j")
print("            relative to its background-face value; density sector cancels in the diff.")
print("BOUNDARY_PROPERTY_USED: Omega -> 0 and (dOmega/deps)|_0 -> 0 as w -> +/- oo")
print("  (drops [Omega*j_w] and [dOm0*(j_w - j_w_face)]).")
print("BOUNDARY_CHECK: Stage0 IBP residual =", ibp_check_2, "; Stage3 second-IBP residual after doit =", second_ibp_local_doit)
