# Independent orchestrator derivation (rule 13) of S11B_ZPERM_SLICE_MAP.
# Derives the mu_s=0 reduction coefficient Lambda_p from the SUPPLIED spec equations ONLY.
# Spec S11b_SHARED_PHYSICS.md §4:  A_pm = mu_s - dp_pm/rho_m ;  J_pm = Lambda_A*A_pm + Lambda_V*V_pm
# The raw-pressure closure (what Lambda_p coefficients) is:  J_pm = Lambda_p*dp_pm + Lambda_V*V_pm
import sympy as sp

Lambda_A, Lambda_V, rho_m, mu_s, dp, V = sp.symbols(
    'Lambda_A Lambda_V rho_m mu_s dp V', real=True)

# Supplied affinity and flux (verbatim from §4):
A   = mu_s - dp/rho_m
J   = Lambda_A*A + Lambda_V*V

# mu_s = 0 slice (slab-side contribution off):
J_slice = J.subs(mu_s, 0)
print("J on mu_s=0 slice            :", sp.expand(J_slice))

# The raw-pressure closure the reduction maps onto: J = Lambda_p*dp + Lambda_V*V.
# Lambda_p is, by definition, the coefficient of dp in J on the slice:
Lambda_p = sp.expand(J_slice).coeff(dp)
print("Lambda_p = coeff of dp in J  :", Lambda_p)

# Cross-check the velocity channel is untouched (spec: 'Lambda_V V unchanged'):
print("coeff of V in J on slice     :", sp.expand(J_slice).coeff(V), " (must equal Lambda_V)")

# The claimed relation to test:
print("Lambda_p + Lambda_A/rho_m    :", sp.simplify(Lambda_p + Lambda_A/rho_m), " (0 => Lambda_p = -Lambda_A/rho_m)")
print("Lambda_p - Lambda_A/rho_m    :", sp.simplify(Lambda_p - Lambda_A/rho_m), " (0 => Lambda_p = +Lambda_A/rho_m)")
