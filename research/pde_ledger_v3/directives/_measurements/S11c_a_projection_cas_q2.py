#!/usr/bin/env python3
"""
Q2 — Does the analytic chain-rule window shape-derivative equal the
differentiate-under-the-integral form?

Object: I(ε) = ∫_{-∞}^{∞} O(G+(ε), G-(ε)) dw
  G+ = w - h+(x,t;ε)
  G- = h-(x,t;ε) - w

(a) Analytic chain rule through G+, G−, then integrate the resulting integrand.
(b) Differentiate under the integral sign, keeping ∫ formal.

Also expands the in-plane profile dependence (η-jet of W_bg) to expose the
second-derivative × profile terms that appear in an analytic expansion about
the flat background faces.

Prints both operands and the difference. No conclusions in the script.
"""
from __future__ import annotations

import sympy as sp

x1, x2, x3, w, t = sp.symbols("x1 x2 x3 w t", real=True)
eps, eta = sp.symbols("epsilon eta", real=True)
W0 = sp.symbols("W0", positive=True)
spatial = (x1, x2, x3)
bounds = (w, -sp.oo, sp.oo)

zeta_c = sp.Function("zeta_c")(*spatial, t)
delta_W = sp.Function("delta_W")(*spatial, t)
w1 = sp.Function("w1")(*spatial)

# §3a: h_s = s W_bg/2 + ζ_s, ζ_s = ζ_c + s δW/2, W_bg = W0 (1 + η w1)
W_bg = W0 * (1 + eta * w1)
zeta_plus = eps * (zeta_c + delta_W / 2)
zeta_minus = eps * (zeta_c - delta_W / 2)
h_plus = W_bg / 2 + zeta_plus
h_minus = -W_bg / 2 + zeta_minus

O = sp.Function("O")
G_plus = w - h_plus
G_minus = h_minus - w
Omega = O(G_plus, G_minus)

# ---------------------------------------------------------------------------
# (b) Differentiate under the integral: dI/dε = ∫ ∂ε Ω dw  (formal).
# (a) Analytic chain rule through G+, G−.
# ---------------------------------------------------------------------------
dOmega_deps = sp.diff(Omega, eps)

g1, g2 = sp.symbols("g_plus g_minus", real=True)
O_1 = sp.diff(O(g1, g2), g1).subs({g1: G_plus, g2: G_minus})
O_2 = sp.diff(O(g1, g2), g2).subs({g1: G_plus, g2: G_minus})
dOmega_chain = O_1 * sp.diff(G_plus, eps) + O_2 * sp.diff(G_minus, eps)

diff_ab = sp.expand(dOmega_deps - dOmega_chain)
diff_ab_at0 = sp.expand(dOmega_deps.subs(eps, 0) - dOmega_chain.subs(eps, 0))

I_form_b = sp.Integral(dOmega_deps, bounds)
I_form_a = sp.Integral(dOmega_chain, bounds)
I_form_b_at0 = sp.Integral(sp.expand(dOmega_deps.subs(eps, 0)), bounds)
I_form_a_at0 = sp.Integral(sp.expand(dOmega_chain.subs(eps, 0)), bounds)

print("=== Q2 STAGE 0: chain rule vs raw ∂ε of Omega (integrands) ===")
print("dOmega/deps (raw) =")
print(sp.expand(dOmega_deps))
print()
print("dOmega/deps (chain rule through G+,G-) =")
print(sp.expand(dOmega_chain))
print()
print("raw - chain =", diff_ab)
print("integrands_identical:", diff_ab == 0)
print("integrands_identical_at_eps0:", diff_ab_at0 == 0)
print()

print("=== Q2 STAGE 1: formal integrals (a) and (b) ===")
print("form_(b)_diff_under_integral =")
print(I_form_b)
print("form_(a)_chain_then_Integral =")
print(I_form_a)
print("Integral(raw - chain) =", sp.Integral(diff_ab, bounds))
print("Integral_difference_identically_zero:", diff_ab == 0)
print()

print("=== Q2 STAGE 2: evaluate at eps=0 (first shape order) ===")
print("integrand_(b)|_eps0 =")
print(sp.expand(dOmega_deps.subs(eps, 0)))
print()
print("integrand_(a)|_eps0 =")
print(sp.expand(dOmega_chain.subs(eps, 0)))
print()
print("form_(b)|_eps0 =", I_form_b_at0)
print("form_(a)|_eps0 =", I_form_a_at0)
print("(a)-(b) integrand at eps0 =", diff_ab_at0)
print()

# ---------------------------------------------------------------------------
# η-expansion about flat faces: exposes O_ij × w1 profile terms.
# ---------------------------------------------------------------------------
G_plus_flat = w - W0 / 2
G_minus_flat = -w - W0 / 2

dG_plus_eta = sp.diff(G_plus, eta).subs({eps: 0, eta: 0})
dG_minus_eta = sp.diff(G_minus, eta).subs({eps: 0, eta: 0})
dG_plus_eps0 = sp.diff(G_plus, eps).subs({eps: 0, eta: 0})
dG_minus_eps0 = sp.diff(G_minus, eps).subs({eps: 0, eta: 0})

dOm_eps_exact_eta = sp.diff(Omega, eps).subs(eps, 0)
dOm_series = sp.expand(sp.series(dOm_eps_exact_eta, eta, 0, 2).removeO())

xi1, xi2 = sp.Dummy("xi_1"), sp.Dummy("xi_2")
O1_f = sp.Subs(sp.diff(O(xi1, xi2), xi1), (xi1, xi2), (G_plus_flat, G_minus_flat))
O2_f = sp.Subs(sp.diff(O(xi1, xi2), xi2), (xi1, xi2), (G_plus_flat, G_minus_flat))
O11_f = sp.Subs(sp.diff(O(xi1, xi2), xi1, 2), (xi1, xi2), (G_plus_flat, G_minus_flat))
O22_f = sp.Subs(sp.diff(O(xi1, xi2), xi2, 2), (xi1, xi2), (G_plus_flat, G_minus_flat))
O12_f = sp.Subs(sp.diff(O(xi1, xi2), xi1, xi2), (xi1, xi2), (G_plus_flat, G_minus_flat))

analytic_O_eta = sp.expand(
    (O1_f + eta * (O11_f * dG_plus_eta + O12_f * dG_minus_eta)) * dG_plus_eps0
    + (O2_f + eta * (O12_f * dG_plus_eta + O22_f * dG_minus_eta)) * dG_minus_eps0
)

print("=== Q2 STAGE 3: analytic η-expansion of ∂εΩ|_ε0 (exposes O_ij × profile) ===")
print("dG+/deta|_0 =", dG_plus_eta)
print("dG-/deta|_0 =", dG_minus_eta)
print("dG+/deps|_0 =", dG_plus_eps0)
print("dG-/deps|_0 =", dG_minus_eps0)
print()
print("series_eta(dOmega/deps|_eps0) to O(eta) =")
print(dOm_series)
print()
print("analytic_chain_Taylor_O(eta) =")
print(analytic_O_eta)
print()
# Formal Subs/Derivative spellings from series vs manual Taylor need not match
# byte-for-byte. Verify on a concrete smooth O that both agree as functions.
O_conc = sp.exp(-(g1**2) - (g2**2))
Omega_conc = O_conc.subs({g1: G_plus, g2: G_minus})
dOm_conc = sp.diff(Omega_conc, eps).subs(eps, 0)
series_conc = sp.expand(sp.series(dOm_conc, eta, 0, 2).removeO())
# Manual Taylor with concrete partials at flat args:
def conc_partial(*vs):
    e = O_conc
    for v in vs:
        e = sp.diff(e, v)
    return e.subs({g1: G_plus_flat, g2: G_minus_flat})


analytic_conc = sp.expand(
    (
        conc_partial(g1)
        + eta
        * (conc_partial(g1, g1) * dG_plus_eta + conc_partial(g1, g2) * dG_minus_eta)
    )
    * dG_plus_eps0
    + (
        conc_partial(g2)
        + eta
        * (conc_partial(g1, g2) * dG_plus_eta + conc_partial(g2, g2) * dG_minus_eta)
    )
    * dG_minus_eps0
)
taylor_diff_conc = sp.simplify(series_conc - analytic_conc)
print("concrete_O=exp(-(G+)^2-(G-)^2): series - analytic_Taylor =", taylor_diff_conc)
print("eta_expansion_matches_chain_Taylor_on_concrete_O:", taylor_diff_conc == 0)
print()

dOm_eta1 = sp.expand(analytic_O_eta.coeff(eta))
d_dx1_eta1 = sp.expand(sp.diff(dOm_eta1.doit(), x1))
has_second = any(
    isinstance(a, sp.Derivative)
    and getattr(a.expr, "func", None) == O
    and sum(n for _, n in a.variable_count) >= 2
    for a in sp.preorder_traversal(d_dx1_eta1)
)
has_w1_grad = d_dx1_eta1.has(sp.diff(w1, x1))

print("=== Q2 STAGE 4: in-plane derivative of the η¹ shape integrand ===")
print("eta^1 coefficient of dOmega/deps|_eps0 =")
print(dOm_eta1)
print()
print("d/dx1 of that coefficient (after doit) =")
print(d_dx1_eta1)
print()
print("contains_second_derivative_of_O:", has_second)
print("contains_in_plane_profile_gradient_dw1/dx1:", has_w1_grad)
print()

# Concrete decaying window: commute check.
print("=== Q2 STAGE 5: concrete decaying window — commute d/deps and Integral ===")
A, B = sp.symbols("A B", real=True)
h_plus_c = W0 / 2 + eps * A
h_minus_c = -W0 / 2 + eps * B
Gp_c = w - h_plus_c
Gm_c = h_minus_c - w
Oc = sp.exp(-(Gp_c**2) - (Gm_c**2))
dOc = sp.diff(Oc, eps).subs(eps, 0)
dOc_chain = (
    sp.diff(sp.exp(-(g1**2) - (g2**2)), g1).subs({g1: Gp_c, g2: Gm_c}) * sp.diff(Gp_c, eps)
    + sp.diff(sp.exp(-(g1**2) - (g2**2)), g2).subs({g1: Gp_c, g2: Gm_c}) * sp.diff(Gm_c, eps)
)
dOc_chain0 = sp.expand(dOc_chain.subs(eps, 0))
print("concrete raw-chain integrand difference at eps0 =", sp.simplify(dOc - dOc_chain0))
print("concrete_integrands_identical:", sp.simplify(dOc - dOc_chain0) == 0)

I_closed_eps = sp.integrate(Oc, (w, -sp.oo, sp.oo))
dI_closed = sp.simplify(sp.diff(I_closed_eps, eps).subs(eps, 0))
I_from_integrand = sp.simplify(sp.integrate(dOc, (w, -sp.oo, sp.oo)))
print("I(eps) closed form =", I_closed_eps)
print("dI/deps|_0 from closed form =", dI_closed)
print("Integral(dOc/deps|_0) from under-integral =", I_from_integrand)
print("closed_minus_under_integral =", sp.simplify(dI_closed - I_from_integrand))
print("commute_on_concrete_decaying_window:", sp.simplify(dI_closed - I_from_integrand) == 0)
print()

print("=== Q2 STAGE 6: boundary / decay property check (concrete) ===")
print("limit w->+oo Oc|_eps0 =", sp.limit(Oc.subs(eps, 0), w, sp.oo))
print("limit w->-oo Oc|_eps0 =", sp.limit(Oc.subs(eps, 0), w, -sp.oo))
print("limit w->+oo dOc/deps|_0 =", sp.limit(dOc, w, sp.oo))
print("limit w->-oo dOc/deps|_0 =", sp.limit(dOc, w, -sp.oo))
print()

print("=== Q2 LITERAL RESIDUAL REPORT ===")
print("DIFF_integrand_(a)_minus_(b) =", diff_ab)
print("DIFF_IDENTICALLY_ZERO:", diff_ab == 0)
print("DIFF_Integral_(a)_minus_(b) =", sp.Integral(diff_ab, bounds))
print("CHAIN_RULE_EXACT_FOR_O_of_Gplus_Gminus:", diff_ab == 0)
print("d/deps_AND_Integral_COMMUTE_FORMALLY:", diff_ab == 0)
print(
    "d/deps_AND_Integral_COMMUTE_ON_CONCRETE_DECAYING_O:",
    sp.simplify(dI_closed - I_from_integrand) == 0,
)
print("PROPERTY_USED: O -> 0 as |arguments| -> +oo (equivalently |w|->oo at fixed faces);")
print("  used to justify dominated/compact-support passage of d/deps under ∫dw")
print("  and verified on O=exp(-(G+)^2-(G-)^2) by matching closed-form dI/deps.")
print("ETA_PROFILE_SECOND_DERIV_TERMS_ARE_PART_OF_THE_SAME_CHAIN_RULE:")
print("  concrete series-analytic residual =", taylor_diff_conc)
print("  eta_expansion_matches_chain_Taylor_on_concrete_O:", taylor_diff_conc == 0)
