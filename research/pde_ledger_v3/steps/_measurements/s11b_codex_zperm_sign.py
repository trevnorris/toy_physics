#!/usr/bin/env python3
"""Independent S11b affinity-to-raw-pressure coefficient derivation."""

import sympy as sp


J, delta_p, mu_s, V = sp.symbols("J delta_p mu_s V")
Lambda_A0, Lambda_V0, Lambda_p0, rho_m = sp.symbols(
    "Lambda_A0 Lambda_V0 Lambda_p0 rho_m", nonzero=True
)

affinity = mu_s - delta_p / rho_m
canonical_rhs = Lambda_A0 * affinity + Lambda_V0 * V
slice_rhs = sp.expand(canonical_rhs.subs(mu_s, 0))
closure_coefficient = sp.expand(slice_rhs).coeff(delta_p)
zero_form = sp.expand(J - slice_rhs)
residual_coefficient = zero_form.coeff(delta_p)
mapping = sp.solve(sp.Eq(Lambda_p0, closure_coefficient), Lambda_p0)[0]

assert sp.simplify(closure_coefficient + Lambda_A0 / rho_m) == 0
assert sp.simplify(residual_coefficient - Lambda_A0 / rho_m) == 0

print(f"affinity = {affinity}")
print(f"mu_s=0 canonical RHS = {slice_rhs}")
print(f"coefficient in J=RHS = {closure_coefficient}")
print(f"zero-form residual J-RHS = {zero_form}")
print(f"coefficient of delta_p in residual = {residual_coefficient}")
print(f"forced map Lambda_p0 = {mapping}")
