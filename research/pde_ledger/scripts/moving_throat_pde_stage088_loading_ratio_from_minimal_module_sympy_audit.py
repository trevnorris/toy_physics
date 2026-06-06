#!/usr/bin/env python3
"""
moving_throat_pde_stage088_loading_ratio_from_minimal_module_sympy_audit.py

SymPy audit for Stage 088.

Checks:
1. exact contact-plus-pole decomposition in terms of rho_alpha = alpha_req/alpha_mix;
2. exact inverse map from (c0,c1) to rho_alpha and zeta_req;
3. match to the minimal isotropic quadrupole module c0=3/4, c1=1/4;
4. exact regime classification Pi_tr = (4/3) C_mix.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 088 — LOADING-RATIO EXTRACTION FROM THE MINIMAL ISOTROPIC MODULE")

rho_alpha, omega, Omega_Q = sp.symbols("rho_alpha omega Omega_Q", positive=True, real=True)
alpha_req, alpha_mix = sp.symbols("alpha_req alpha_mix", positive=True, real=True)
c0, c1 = sp.symbols("c0 c1", positive=True, real=True)

Y_loading = sp.simplify(alpha_mix/alpha_req + (alpha_req-alpha_mix)/alpha_req/(1 - omega**2/Omega_Q**2))
Y_rho = sp.simplify(1/rho_alpha + (rho_alpha-1)/rho_alpha/(1 - omega**2/Omega_Q**2))

print("Y_loading(omega) =", Y_loading)
print("Y_rho(omega)     =", Y_rho)

expect_zero(
    "loading form - rho form",
    Y_loading.subs(alpha_req, rho_alpha*alpha_mix) - Y_rho
)

# Recover (c0, c1) from Y_rho by independent paths, rather than defining them.
# Strategy: extract c1 as the pole residue at u = omega^2/Omega_Q^2 -> 1, then
# recover c0 by subtraction at u = 0. This avoids the tautological
# "define c0 := 1/rho_alpha then assert c0 == 1/rho_alpha" pattern.
#
# Substitute omega**2 = u * Omega_Q**2 (Y_rho depends on omega only through
# omega**2 so this gives a clean rational function of u). Subbing on the
# combined ratio `omega**2/Omega_Q**2 -> u` fails after sp.simplify reshapes
# the denominator.
u = sp.symbols("u", real=True)
Y_rho_u = sp.simplify(Y_rho.subs(omega**2, u * Omega_Q**2))
c1_extracted = sp.simplify(sp.limit((1 - u) * Y_rho_u, u, 1))
c0_extracted = sp.simplify(Y_rho_u.subs(u, 0) - c1_extracted)

print("c0 extracted from Y_rho (subtract pole at u=0) =", c0_extracted)
print("c1 extracted from Y_rho (pole residue at u=1)  =", c1_extracted)

expect_zero("c0_extracted - 1/rho_alpha", c0_extracted - 1/rho_alpha)
expect_zero(
    "c1_extracted - (rho_alpha-1)/rho_alpha",
    c1_extracted - (rho_alpha - 1)/rho_alpha,
)

# With c0, c1 now derived (not defined), the contact-plus-pole reconstruction
# is a real check.
expect_zero(
    "contact-plus-pole reconstruction",
    Y_rho - (c0_extracted + c1_extracted/(1 - omega**2/Omega_Q**2))
)
expect_zero("c0 + c1 - 1", c0_extracted + c1_extracted - 1)

rho_from_c0 = sp.simplify(1/c0)
rho_from_c1 = sp.simplify(1/(1-c1))
zeta_from_c = sp.simplify(c1/c0)

print("rho_alpha from c0 =", rho_from_c0)
print("rho_alpha from c1 =", rho_from_c1)
print("zeta_req from (c0,c1) =", zeta_from_c)

expect_zero("rho(c0(rho)) - rho", rho_from_c0.subs(c0, c0_extracted) - rho_alpha)
expect_zero("rho(c1(rho)) - rho", rho_from_c1.subs(c1, c1_extracted) - rho_alpha)
expect_zero(
    "zeta(c(rho)) - (rho-1)",
    zeta_from_c.subs({c0: c0_extracted, c1: c1_extracted}) - (rho_alpha - 1),
)

# Minimal isotropic quadrupole module: paper Y_Q^cons = 3/4 + (1/4)/(1-omega^2/Omega_Q^2).
# Extract (c0, c1) from the paper form by the same independent path; do NOT
# substitute fractions into pre-defined symbols.
Y_Q_paper = sp.Rational(3, 4) + sp.Rational(1, 4)/(1 - omega**2/Omega_Q**2)
Y_Q_paper_u = sp.simplify(Y_Q_paper.subs(omega**2, u * Omega_Q**2))
c1_paper = sp.simplify(sp.limit((1 - u) * Y_Q_paper_u, u, 1))
c0_paper = sp.simplify(Y_Q_paper_u.subs(u, 0) - c1_paper)
expect_zero("c0_paper - 3/4", c0_paper - sp.Rational(3, 4))
expect_zero("c1_paper - 1/4", c1_paper - sp.Rational(1, 4))

rho_min = sp.simplify(rho_from_c0.subs(c0, c0_paper))
zeta_min = sp.simplify(zeta_from_c.subs({c0: c0_paper, c1: c1_paper}))

print("rho_alpha(minimal isotropic module) =", rho_min)
print("zeta_req(minimal isotropic module) =", zeta_min)

expect_zero("rho_min - 4/3", rho_min - sp.Rational(4, 3))
expect_zero("zeta_min - 1/3", zeta_min - sp.Rational(1, 3))

# Product-language regime: Pi_tr = rho_alpha * C_mix is the stage-085 product
# identity (verified upstream in scripts/moving_throat_pde_stage085_*). Substitute
# rho_min (derived above from c_contact = 3/4) and confirm Pi_tr = (4/3) C_mix.
C_mix = sp.symbols("C_mix", positive=True, real=True)
Pi_tr_from_rho = rho_min * C_mix
expect_zero("Pi_tr_from_rho - (4/3) C_mix", Pi_tr_from_rho - sp.Rational(4, 3)*C_mix)

assert rho_min > 1, "rho_min must exceed 1 for the mixed-baseline-not-enough regime"
assert rho_min < 2, "rho_min must lie below 2 to stay in the symmetric-lowest-twin regime"

print("\nRegime classification:")
print("  Pi_tr = rho_min * C_mix = (4/3) C_mix")
print("  therefore C_mix < Pi_tr < 2 C_mix")
print("  and zeta_req = 1/3 < 1, so the symmetric lowest twin suffices.")
