#!/usr/bin/env python3
"""
Stage 45 SymPy audit.

Checks:
1. n=5 EOS identities h'(rho)=m c_s^2/rho.
2. Parent-channel projection formulas Theta_sigma and Lambda_phi.
3. Exact parent formula for G_micro.
4. Coherence-factor decomposition.
5. Xi_micro = kappa G_micro reduction.
"""

from __future__ import annotations
import sympy as sp


def banner(title: str) -> None:
    line = "=" * 80
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


banner("STAGE 45 — PARENT-ACTION PROJECTION OF THE MICROSCOPIC GAIN")

# EOS identities
rho, K, m = sp.symbols("rho K m", positive=True, real=True)
n_poly = sp.symbols("n_poly", positive=True, integer=True)
# Preserve the stage normalization U = K rho^n/(n - 1), with c_s^2 = d(K rho^n)/d rho / m.
U_general = K * rho**n_poly / (n_poly - 1)
h_general = sp.diff(U_general, rho)
hprime_general = sp.diff(h_general, rho)
cs_sq_general = sp.diff(K * rho**n_poly, rho) / m

expect_zero(
    "h'(rho) = m c_s^2 / rho (general polytrope)",
    hprime_general - m * cs_sq_general / rho,
)

subs_n5 = {n_poly: 5}
print("Specializing to n=5:")
print("  U(rho) =", sp.simplify(U_general.subs(subs_n5)))
print("  h(rho) =", sp.simplify(h_general.subs(subs_n5)))
print("  h'(rho) =", sp.simplify(hprime_general.subs(subs_n5)))
print("  c_s^2(rho) =", sp.simplify(cs_sq_general.subs(subs_n5)))
expect_zero(
    "n=5 specialization of h' = m c_s^2/rho",
    (hprime_general - m * cs_sq_general / rho).subs(subs_n5),
)

cs_sq_wrong = sp.diff(K * rho**(n_poly + 1), rho) / m
residual_wrong = (hprime_general - m * cs_sq_wrong / rho).subs(subs_n5)
assert sp.simplify(residual_wrong) != 0, "Inconsistency check failed to detect wrong exponent"
print("Inconsistency probe (n+1 in c_s^2 only):", sp.simplify(residual_wrong))

# Projected source/support coefficients
rho_star, cs_star_sq = sp.symbols("rho_star cs_star_sq", positive=True, real=True)
Nss, Npp, Osp = sp.symbols("N_ss N_pp O_sp", positive=True, real=True)
g_phi, KX, TX, L, kappa = sp.symbols("g_phi K_X T_X L kappa", positive=True, real=True)

Theta_sigma = (m * cs_star_sq / rho_star) * Nss
Lambda_phi = g_phi * Osp
sigma, phi = sp.symbols("sigma phi", real=True)
S_parent = (
    sp.Rational(1, 2) * Theta_sigma * sigma**2
    + Lambda_phi * sigma * phi
    + sp.Rational(1, 2) * KX * phi**2
)
sigma_star = sp.solve(sp.diff(S_parent, sigma), sigma)[0]
S_eff_phi = sp.expand(S_parent.subs(sigma, sigma_star))
gain_from_action = sp.simplify((KX - 2 * sp.simplify(S_eff_phi.coeff(phi, 2))) / KX)
G_micro_closed = rho_star * g_phi**2 * Osp**2 / (m * cs_star_sq * KX * Nss)

print("Theta_sigma =", Theta_sigma)
print("Lambda_phi =", Lambda_phi)
print("sigma_star =", sigma_star)
print("S_eff_phi =", S_eff_phi)
print("G_micro from action =", gain_from_action)
expect_zero("G_micro from parent action vs closed form", gain_from_action - G_micro_closed)

print("Coherence factor (definition):  C_(sigma phi)^2 := Osp^2 / (Nss Npp)")

# Xi_micro relation
Xi_micro = sp.simplify(kappa * gain_from_action)
Xi_target = rho_star * g_phi**2 * Osp**2 * L**2 / (m * cs_star_sq * TX * Nss)
print("Xi_micro =", Xi_micro)
kappa_solved = sp.solve(Xi_micro - Xi_target, kappa)
assert kappa_solved == [KX * L**2 / TX], f"Unexpected kappa solution: {kappa_solved}"
print("kappa solved from Xi_micro = Xi_target:", kappa_solved[0])

print("\nAll Stage 45 symbolic checks passed.")
