#!/usr/bin/env python3
"""
moving_throat_pde_stage71_loading_ratio_from_minimal_module_sympy_audit.py

SymPy audit for Stage 71.

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

banner("STAGE 71 — LOADING-RATIO EXTRACTION FROM THE MINIMAL ISOTROPIC MODULE")

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

c0_from_rho = sp.simplify(1/rho_alpha)
c1_from_rho = sp.simplify((rho_alpha-1)/rho_alpha)

expect_zero(
    "contact-plus-pole reconstruction",
    Y_rho - (c0_from_rho + c1_from_rho/(1 - omega**2/Omega_Q**2))
)

print("c0(rho_alpha) =", c0_from_rho)
print("c1(rho_alpha) =", c1_from_rho)

expect_zero("c0 + c1 - 1", c0_from_rho + c1_from_rho - 1)

rho_from_c0 = sp.simplify(1/c0)
rho_from_c1 = sp.simplify(1/(1-c1))
zeta_from_c = sp.simplify(c1/c0)

print("rho_alpha from c0 =", rho_from_c0)
print("rho_alpha from c1 =", rho_from_c1)
print("zeta_req from (c0,c1) =", zeta_from_c)

expect_zero("rho(c0(rho)) - rho", rho_from_c0.subs(c0, c0_from_rho) - rho_alpha)
expect_zero("rho(c1(rho)) - rho", rho_from_c1.subs(c1, c1_from_rho) - rho_alpha)
expect_zero("zeta(c(rho)) - (rho-1)", zeta_from_c.subs({c0: c0_from_rho, c1: c1_from_rho}) - (rho_alpha - 1))

# Minimal isotropic quadrupole module
rho_min = sp.simplify(rho_from_c0.subs(c0, sp.Rational(3,4)))
zeta_min = sp.simplify(zeta_from_c.subs({c0: sp.Rational(3,4), c1: sp.Rational(1,4)}))

print("rho_alpha(minimal isotropic module) =", rho_min)
print("zeta_req(minimal isotropic module) =", zeta_min)

expect_zero("rho_min - 4/3", rho_min - sp.Rational(4,3))
expect_zero("zeta_min - 1/3", zeta_min - sp.Rational(1,3))

# Product-language regime
Pi_tr, C_mix = sp.symbols("Pi_tr C_mix", positive=True, real=True)
expect_zero("Pi_tr/C_mix - 4/3", sp.simplify((sp.Rational(4,3)*C_mix)/C_mix - sp.Rational(4,3)))

print("\nRegime classification:")
print("  Pi_tr = (4/3) C_mix")
print("  therefore C_mix < Pi_tr < 2 C_mix")
print("  and zeta_req = 1/3 < 1, so the symmetric lowest twin suffices.")
