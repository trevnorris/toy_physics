#!/usr/bin/env python3
"""CAS comparison for S11c-a Q2: two forms of the window shape derivative."""

import sympy as sp


eps, w = sp.symbols("epsilon w", real=True)
h_plus_0, h_minus_0 = sp.symbols("h_plus_0 h_minus_0", real=True)
v_plus, v_minus = sp.symbols("v_plus v_minus", real=True)
a, b = sp.symbols("a b", real=True)
O = sp.Function("O")

# h_s(epsilon)=h_s^0+epsilon*v_s and G_s=s(w-h_s).
G_plus = w - h_plus_0 - eps * v_plus
G_minus = -w + h_minus_0 + eps * v_minus
Omega_eps = O(G_plus, G_minus)
bounds = (w, -sp.oo, sp.oo)

# Route A: expand analytically by the two-argument chain rule, evaluate at
# epsilon=0, and then wrap the result in the (formal) w integral.
G_plus_0 = G_plus.subs(eps, 0)
G_minus_0 = G_minus.subs(eps, 0)
chain_integrand = (
    -v_plus * sp.diff(O(a, b), a) + v_minus * sp.diff(O(a, b), b)
).subs({a: G_plus_0, b: G_minus_0})
route_a = sp.Integral(chain_integrand, bounds)

# Route B: ask SymPy to differentiate the unevaluated Integral.  Its fixed
# limits cause SymPy to differentiate under the integral while leaving the
# w integral unevaluated.
formal_integral = sp.Integral(Omega_eps, bounds)
route_b = sp.diff(formal_integral, eps).subs(eps, 0)

# Compare both the integrands and the two formal Integral objects.
under_integral = sp.diff(Omega_eps, eps).subs(eps, 0)
pointwise_difference = sp.simplify(chain_integrand - under_integral)
formal_difference = sp.simplify(route_a - route_b)

print("Q2_S11c_a_WINDOW_SHAPE_DERIVATIVE")
print("G_PLUS =", G_plus)
print("G_MINUS =", G_minus)
print("ROUTE_A_CHAIN_INTEGRAND =", chain_integrand)
print("ROUTE_A_FORMAL_INTEGRAL =", route_a)
print("ROUTE_B_DIFFERENTIATE_UNDER_INTEGRAL =", route_b)
print("POINTWISE_INTEGRAND_DIFFERENCE_A_MINUS_B =", pointwise_difference)
print("FORMAL_INTEGRAL_DIFFERENCE_A_MINUS_B =", formal_difference)
print("CHAIN_RULE_EXPANSION_EXACT =", pointwise_difference == 0)
print("DERIVATIVE_AND_INTEGRAL_COMMUTE =", formal_difference == 0)
print("COMMUTATION_CONDITION = fixed limits and an integrable epsilon-uniform dominator for dO/depsilon (uniform compact support suffices)")
