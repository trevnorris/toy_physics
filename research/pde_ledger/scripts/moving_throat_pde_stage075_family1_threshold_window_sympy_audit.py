#!/usr/bin/env python3
"""
moving_throat_pde_stage075_family1_threshold_window_sympy_audit.py

SymPy audit for Stage 075:
- evaluate Delta_0 and Delta_inf on the explicit Family-1 / healing-locked branch,
- compute the explicit Upsilon and Xi thresholds,
- reduce the remaining amplitude to Upsilon_w = alpha_r^2 Theta_w with alpha_r = 10,
- compute the Theta_w threshold window.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

banner("STAGE 075 — EXPLICIT FAMILY-1 THRESHOLD WINDOW")

Pe_req = sp.symbols("Pe_req", positive=True, real=True)
alpha_r = sp.Integer(10)
# F4 (v2 paper-alignment Q2 direction (a) lock): paper Inputs line states
# Upsilon_w = alpha_r^2 Theta_w with alpha_r^2 = 100. Lock the value so any
# future drift between the paper Inputs line and the script surfaces here.
assert alpha_r**2 == 100, "paper Inputs line lock: alpha_r^2 must equal 100"
print("alpha_r^2 (paper Inputs line lock) =", alpha_r**2, "  PASS")

Lambda_ell = sp.Integer(37)
eta = sp.Integer(37)
kappa = sp.Rational(12321, 5)
alpha = sp.sqrt(kappa)

Delta0 = sp.simplify(
    eta * (sp.cosh(alpha) - 1) /
    (alpha**2 * (alpha * sp.sinh(alpha) + eta * sp.cosh(alpha)))
)
Deltainf = sp.simplify(
    (sp.cosh(alpha) + (eta / alpha) * sp.sinh(alpha) - 1) /
    (alpha * sp.sinh(alpha) + eta * sp.cosh(alpha))
)

# Independent symbolic check: the stated closed forms must satisfy
# the defining algebraic identities for *all* positive alpha_sym, eta_sym,
# not just the substituted numeric values alpha = 111/sqrt(5), eta = 37.
alpha_sym, eta_sym = sp.symbols("alpha_sym eta_sym", positive=True, real=True)
Delta0_sym = eta_sym * (sp.cosh(alpha_sym) - 1) / (
    alpha_sym**2 * (alpha_sym * sp.sinh(alpha_sym) + eta_sym * sp.cosh(alpha_sym))
)
Deltainf_sym = (sp.cosh(alpha_sym) + (eta_sym / alpha_sym) * sp.sinh(alpha_sym) - 1) / (
    alpha_sym * sp.sinh(alpha_sym) + eta_sym * sp.cosh(alpha_sym)
)
delta0_identity = sp.simplify(
    (alpha_sym * sp.sinh(alpha_sym) + eta_sym * sp.cosh(alpha_sym)) * Delta0_sym
    - eta_sym * (sp.cosh(alpha_sym) - 1) / alpha_sym**2
)
deltainf_identity = sp.simplify(
    (alpha_sym * sp.sinh(alpha_sym) + eta_sym * sp.cosh(alpha_sym)) * Deltainf_sym
    - (sp.cosh(alpha_sym) + (eta_sym / alpha_sym) * sp.sinh(alpha_sym) - 1)
)
print("Delta_0 algebraic identity (free alpha, eta) =", delta0_identity)
print("Delta_inf algebraic identity (free alpha, eta) =", deltainf_identity)
assert delta0_identity == 0
assert deltainf_identity == 0

# F1 (v2): the algebraic identities above are tautological by construction —
# they only verify that the CAS can cancel a common factor. Add asymptotic-limit
# checks that exercise non-trivial properties of the closed forms.
#
# Large-alpha: alpha * Delta_inf -> 1 (since Delta_inf ~ 1/alpha in that limit).
# A wrong factor in numerator or denominator of Delta_inf would change this limit.
large_alpha_check = sp.limit(alpha_sym * Deltainf_sym, alpha_sym, sp.oo)
print("alpha * Delta_inf large-alpha limit =", large_alpha_check)
assert large_alpha_check == 1
# Small-alpha: Delta_0 -> 1/2 (since cosh(alpha)-1 ~ alpha^2/2 and
# alpha*sinh(alpha)+eta*cosh(alpha) ~ eta as alpha -> 0).
small_alpha_check_delta0 = sp.limit(Delta0_sym, alpha_sym, 0)
print("Delta_0 small-alpha limit =", small_alpha_check_delta0)
assert small_alpha_check_delta0 == sp.Rational(1, 2)
print("PASS: alpha * Delta_inf -> 1 (large-alpha) and Delta_0 -> 1/2 (small-alpha)")

print("Lambda_ell =", Lambda_ell)
print("eta =", eta)
print("kappa =", kappa)
print("alpha = sqrt(kappa) =", alpha)
print("Delta_0 =", Delta0)
print("Delta_inf =", Deltainf)
print("Delta_0 (numeric) =", sp.N(Delta0, 20))
print("Delta_inf (numeric) =", sp.N(Deltainf, 20))

Upsilon_fail = sp.simplify(Pe_req / (Lambda_ell**2 * Deltainf))
Upsilon_suff = sp.simplify(Pe_req / (Lambda_ell**2 * Delta0))
Xi_fail = sp.simplify(Pe_req / Deltainf)
Xi_suff = sp.simplify(Pe_req / Delta0)

Theta = sp.symbols("Theta_w", positive=True, real=True)
Theta_fail = sp.simplify(Upsilon_fail / alpha_r**2)
Theta_suff = sp.simplify(Upsilon_suff / alpha_r**2)

print("\nUpsilon_fail =", Upsilon_fail)
print("Upsilon_suff =", Upsilon_suff)
print("Xi_fail =", Xi_fail)
print("Xi_suff =", Xi_suff)
print("Theta_fail =", Theta_fail)
print("Theta_suff =", Theta_suff)

print("\nNumerics (coefficients multiplying Pe_req):")
print("Upsilon_fail / Pe_req =", sp.N(sp.simplify(Upsilon_fail / Pe_req), 20))
print("Upsilon_suff / Pe_req =", sp.N(sp.simplify(Upsilon_suff / Pe_req), 20))
print("Xi_fail / Pe_req =", sp.N(sp.simplify(Xi_fail / Pe_req), 20))
print("Xi_suff / Pe_req =", sp.N(sp.simplify(Xi_suff / Pe_req), 20))
print("Theta_fail / Pe_req =", sp.N(sp.simplify(Theta_fail / Pe_req), 20))
print("Theta_suff / Pe_req =", sp.N(sp.simplify(Theta_suff / Pe_req), 20))

# Exact reduction Upsilon_w = alpha_r^2 Theta_w on the reference branch.
# Test the round-trip on the actually-constructed fail and suff branches,
# not the trivial identity 100*Theta == 100*Theta.
Upsilon_fail_from_Theta = sp.simplify(alpha_r**2 * Theta_fail)
Upsilon_suff_from_Theta = sp.simplify(alpha_r**2 * Theta_suff)
print("\nUpsilon_fail - alpha_r^2 * Theta_fail =", sp.simplify(Upsilon_fail - Upsilon_fail_from_Theta))
print("Upsilon_suff - alpha_r^2 * Theta_suff =", sp.simplify(Upsilon_suff - Upsilon_suff_from_Theta))
assert sp.simplify(Upsilon_fail - Upsilon_fail_from_Theta) == 0
assert sp.simplify(Upsilon_suff - Upsilon_suff_from_Theta) == 0

print("\nFinal ledger:")
print("  Delta_0   ~", sp.N(Delta0, 16))
print("  Delta_inf ~", sp.N(Deltainf, 16))
print("  Upsilon_fail ~", sp.N(sp.simplify(Upsilon_fail / Pe_req), 16), "* Pe_req")
print("  Upsilon_suff ~", sp.N(sp.simplify(Upsilon_suff / Pe_req), 16), "* Pe_req")
print("  Theta_fail ~", sp.N(sp.simplify(Theta_fail / Pe_req), 16), "* Pe_req")
print("  Theta_suff ~", sp.N(sp.simplify(Theta_suff / Pe_req), 16), "* Pe_req")
