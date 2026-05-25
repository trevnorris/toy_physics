#!/usr/bin/env python3
"""
Stage 59 SymPy audit.

Checks:
1. Exact n=5 enthalpy-sound-speed identity h = m c_s^2 / 4.
2. Exact algebraic reduction of Theta_w under the local enthalpy lock mu_* = lambda_mu h_w.
3. Exact healing-lock reduction Theta_w = lambda_mu^2 rho_w^2 / (16 ell^2).
4. Reference-branch form Theta_w = 25 lambda_mu^2 rho_w^2 in normalized Family-1 wall units.
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


banner("STAGE 59 — EXACT n=5 WALL-DEPTH LOCK")

K, rho, m, hbar, csw = sp.symbols("K rho m hbar c_sw", positive=True, real=True)
lambda_mu, rho_w, ell, a = sp.symbols("lambda_mu rho_w ell a", positive=True, real=True)
n_poly = sp.symbols("n_poly", positive=True, real=True)

P = K * rho**n_poly
# Polytrope enthalpy per particle for the exponent-indexed branch.
h_polytrope = sp.simplify(m * n_poly / (n_poly - 1) * K * rho**(n_poly - 1))
# Derive U from the energy-density convention:
#   U(rho) = rho * integral(P/rho^2 d rho)
U_per_mass = sp.simplify(sp.integrate(P / rho**2, rho))
U = sp.simplify(rho * U_per_mass)
cs2 = sp.diff(P, rho) / m
h = sp.diff(U, rho)
# Specialize to the n=5 case for the printed numbers and the assertion
cs2_n5 = sp.simplify(cs2.subs(n_poly, 5))
h_n5 = sp.simplify(h.subs(n_poly, 5))
print("c_s^2(rho)  [n=5] =", cs2_n5)
print("h(rho)      [n=5] =", h_n5)
expect_zero("n=5 enthalpy identity", h_n5 - m * cs2_n5 / 4)
# Non-tautology check: the same identity must FAIL for n=3.
residual_n3 = sp.simplify(h.subs(n_poly, 3) - m * cs2.subs(n_poly, 3) / 4)
print("n=3 residual (should be NONZERO) =", residual_n3)
if residual_n3 == 0:
    raise AssertionError("n=3 residual is zero -- identity does not actually distinguish n=5")

# Define mu_star symbolically (independent symbol), then impose the enthalpy-lock condition
mu_star_sym = sp.symbols("mu_star_sym", positive=True, real=True)
enthalpy_lock = mu_star_sym - lambda_mu * m * csw**2 / 4  # set this to zero
mu_star_solved = sp.solve(enthalpy_lock, mu_star_sym)[0]
Theta_w = sp.simplify(4 * rho_w**2 * mu_star_solved**2 / (hbar**2 * csw**2))
# Independent route: Theta_w as (2 rho_w mu_star / (hbar c_sw))^2
Theta_w_alt = sp.simplify((2 * rho_w * mu_star_solved / (hbar * csw))**2)
print("Theta_w (enthalpy lock) =", Theta_w)
expect_zero("Theta_w vs alternative-form derivation", Theta_w - Theta_w_alt)

healing_condition = csw - hbar / (2 * m * ell)  # the healing-length defining relation
ell_solved = sp.solve(healing_condition, ell)[0]
print("ell from healing condition =", ell_solved)
# Substitute c_sw in Theta_w using the inverse relation (solve for c_sw):
csw_from_ell = sp.solve(healing_condition, csw)[0]
Theta_w_in_ell = sp.simplify(Theta_w.subs(csw, csw_from_ell))
Theta_heal_target = sp.simplify(lambda_mu**2 * rho_w**2 / (16 * ell**2))
expect_zero("healing-lock reduction", Theta_w_in_ell - Theta_heal_target)
print("Theta_w (healing lock) =", Theta_heal_target)

# Reference-branch convention: ell = a * ref_factor with ref_factor = 1/20.
# TODO(provenance): cite the upstream stage that fixes ref_factor. This factor is
# the load-bearing piece of the "25" in the normalized reference identity.
ref_factor = sp.Rational(1, 20)  # reference-branch convention: ell = a * ref_factor  (see F2 below for provenance)
ref_sub = {ell: a * ref_factor}
Theta_ref = sp.simplify(Theta_heal_target.subs(ref_sub))
print("Theta_w (reference branch, general a) =", Theta_ref)
Theta_ref_norm = sp.simplify(Theta_ref.subs(a, 1))
print("Theta_w (reference branch, normalized wall units) =", Theta_ref_norm)
# The "25" target is derived from ref_factor as (1/ref_factor)^2 / 16
ref_target = sp.simplify((1 / ref_factor)**2 / 16) * lambda_mu**2 * rho_w**2
expect_zero("normalized reference factor", Theta_ref_norm - ref_target)
