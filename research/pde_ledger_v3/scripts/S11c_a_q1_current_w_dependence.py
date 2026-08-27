#!/usr/bin/env python3
"""CAS comparison for S11c-a Q1: bulk-w current versus a face-frozen current.

O_p and O_m in the displayed normal form denote partial derivatives of O
with respect to its first and second arguments.  D1, D2, D3 denote the
in-plane partial derivatives d_{x_i} J_i, with all spectator (x,t)
arguments suppressed.  D_i(h0) is the carried fixed-w partial derivative
evaluated at the background face, not a total derivative of a trace.
"""

import sympy as sp


eps, w = sp.symbols("epsilon w", real=True)
h_plus_0, h_minus_0, h0 = sp.symbols(
    "h_plus_0 h_minus_0 h0", real=True
)
v_plus, v_minus = sp.symbols("v_plus v_minus", real=True)

O = sp.Function("O")
O_p = sp.Function("O_p")
O_m = sp.Function("O_m")
Jw = sp.Function("Jw")
D1, D2, D3 = (sp.Function(name) for name in ("D1", "D2", "D3"))

# G_s = s(w-h_s), h_s(epsilon)=h_s^0+epsilon*v_s.
G_plus = w - h_plus_0 - eps * v_plus
G_minus = -w + h_minus_0 + eps * v_minus
Omega_eps = O(G_plus, G_minus)
Omega_0 = Omega_eps.subs(eps, 0)

# A generic first-order perturbation residual is epsilon*R.  This confirms
# that retaining the dynamic window before differentiation gives Omega_0*R
# in the first perturbation coefficient because the perturbation has no
# zeroth-order current.
R = sp.Symbol("R", commutative=True)
first_order_check = sp.simplify(
    sp.diff(Omega_eps * eps * R, eps).subs(eps, 0) - Omega_0 * R
)

# Direct CAS check of the local integration-by-parts identity in w.
ibp_local_check = sp.simplify(
    Omega_0 * sp.diff(Jw(w), w)
    - (-sp.diff(Omega_0, w) * Jw(w))
    - sp.diff(Omega_0 * Jw(w), w)
)

# The named chain-rule form of d_w Omega_0.  SymPy separately verifies its
# coefficients d_w G_+=1 and d_w G_-=-1.
g_plus_0 = G_plus.subs(eps, 0)
g_minus_0 = G_minus.subs(eps, 0)
dG_check = (
    sp.simplify(sp.diff(g_plus_0, w) - 1),
    sp.simplify(sp.diff(g_minus_0, w) + 1),
)
Omega_w_named = O_p(g_plus_0, g_minus_0) - O_m(g_plus_0, g_minus_0)

div_full = D1(w) + D2(w) + D3(w)
div_face = D1(h0) + D2(h0) + D3(h0)

# These are the current-dependent pieces after integration by parts in w
# and after the endpoint term has been set to zero.  The density-dependent
# part is identical in both choices and cancels from their difference.
integrand_full = Omega_0 * div_full - Omega_w_named * Jw(w)
integrand_face = Omega_0 * div_face - Omega_w_named * Jw(h0)
projection_full = sp.Integral(integrand_full, (w, -sp.oo, sp.oo))
projection_face = sp.Integral(integrand_face, (w, -sp.oo, sp.oo))
difference = sp.Integral(
    sp.collect(sp.expand(integrand_full - integrand_face), (Omega_0, Omega_w_named)),
    (w, -sp.oo, sp.oo),
)
difference_construction_check = sp.simplify(
    (integrand_full - integrand_face) - difference.function
)

# A concrete admissible witness proves the defect is not the zero functional.
# Use the standard smooth slab window
#   O(g_+,g_-)=1/((1+exp(g_+))*(1+exp(g_-))),
# take h_+^0=log(100), h_-^0=-log(100), h0=0, Jw(w)=w, and
# make the tangential divergence defect zero.  The window is about 0.9803 at
# the centre and decays at both ends.  Integration by parts makes the normal
# defect integral equal to integral(Omega,w).  Under y=exp(w), that integral
# is the rational integral below, evaluated by a CAS-checked antiderivative.
y = sp.symbols("y", positive=True)
a_witness = sp.Rational(1, 100)
witness_rational_integrand = 1 / (
    (1 + a_witness * y) * (y + a_witness)
)
witness_antiderivative = (
    sp.log(y + a_witness) - sp.log(1 + a_witness * y)
) / (1 - a_witness**2)
witness_antiderivative_check = sp.simplify(
    sp.diff(witness_antiderivative, y) - witness_rational_integrand
)
witness_normal_only = sp.simplify(
    sp.limit(witness_antiderivative, y, sp.oo)
    - sp.limit(witness_antiderivative, y, 0, dir="+")
)

print("Q1_S11c_a_CURRENT_W_DEPENDENCE")
print("FIRST_PERTURBATION_ORDER_CHECK =", first_order_check)
print("IBP_LOCAL_IDENTITY_CHECK =", ibp_local_check)
print("dG_PLUS_dw_MINUS_1, dG_MINUS_dw_PLUS_1 =", dG_check)
print("NOTATION: O_p=dO/dg_plus; O_m=dO/dg_minus; D_i(w)=partial_xi J_i(w)")
print("FULL_CURRENT_PROJECTION =", projection_full)
print("FACE_FROZEN_PROJECTION =", projection_face)
print("DIFFERENCE_2_MINUS_1 =", difference)
print("DIFFERENCE_CONSTRUCTION_CHECK =", difference_construction_check)
print("BOUNDARY_TERM_BEFORE_DECAY = [O(G_plus_0,G_minus_0)*(Jw(w)-Jw(h0))]_(w=-oo)^(w=oo)")
print("BOUNDARY_PROPERTY_USED = both endpoint products are 0; compact support suffices, or decay must dominate Jw")
print("FROZEN_NORMAL_TERM_AFTER_BOUNDARY = 0")
print("WITNESS_WINDOW = 1/((1+exp(g_plus))*(1+exp(g_minus))); h_plus_0=log(100), h_minus_0=-log(100), Jw(w)=w")
print("WITNESS_RATIONAL_ANTIDERIVATIVE_CHECK =", witness_antiderivative_check)
print("NONZERO_NORMAL_DEFECT_WITNESS =", witness_normal_only)
print("IS_DIFFERENCE_IDENTICALLY_ZERO =", sp.simplify(witness_normal_only) == 0)
