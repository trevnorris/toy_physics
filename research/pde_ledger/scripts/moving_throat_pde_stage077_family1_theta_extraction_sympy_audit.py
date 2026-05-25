#!/usr/bin/env python3
"""
Stage 60 SymPy/mpmath audit.

Checks:
1. Exact cut point xi_* for the Family-1 radial wall profile.
2. Exact canonical support normalization If = 1/3.
3. Numerical shell-weighted moments <rho>_chi and <rho^2>_chi on the alpha_r=10 branch.
4. Numerical effective wall-depth datum Theta_w^(chi) and conservative Jensen floor Theta_w^(J).
"""
from __future__ import annotations

import mpmath as mp
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


banner("STAGE 60 — SHELL-WEIGHTED THETA EXTRACTION")

xi, alpha_r, lambda_mu = sp.symbols("xi alpha_r lambda_mu", positive=True, real=True)
S = (1 + sp.tanh(xi)) / 2
chi = sp.diff(S, xi)
If = sp.integrate(sp.simplify(chi**2), (xi, -sp.oo, sp.oo))
print("chi_phi(xi) =", sp.simplify(chi))
print("I_f         =", sp.simplify(If))
expect_zero("I_f - 1/3", If - sp.Rational(1, 3))

xi_star = sp.simplify(sp.atanh(2 / sp.sqrt(alpha_r) - 1))
print("xi_* =", xi_star)
print("xi_*(alpha_r=10) =", sp.N(xi_star.subs(alpha_r, 10), 20))

S_at_star = (1 + sp.tanh(xi_star)) / 2
rho_quartic_at_star = sp.simplify(1 - alpha_r * S_at_star**2)
expect_zero("1 - alpha_r * S(xi_*)**2", rho_quartic_at_star)

# Numerical branch evaluation at alpha_r = 10.
mp.mp.dps = 50
alpha_num = mp.mpf('10')

def S_num(x):
    return mp.mpf('0.5') * (1 + mp.tanh(x))

def chi_num(x):
    return mp.mpf('0.5') / (mp.cosh(x)**2)

def rho_num(x):
    val = 1 - alpha_num * (S_num(x)**2)
    return val**mp.mpf('0.25') if val > 0 else mp.mpf('0')

xi_cut = mp.atanh(2 / mp.sqrt(alpha_num) - 1)
print("numeric xi_* =", xi_cut)

den = mp.quad(lambda x: chi_num(x)**2, [-mp.inf, -4, 0, mp.inf])
num1 = mp.quad(lambda x: chi_num(x)**2 * rho_num(x), [-mp.inf, -4, xi_cut])
num2 = mp.quad(lambda x: chi_num(x)**2 * rho_num(x)**2, [-mp.inf, -4, xi_cut])
R1 = num1 / den
R2 = num2 / den

print("<rho>_chi    =", R1)
print("<rho^2>_chi  =", R2)
print("denominator  =", den)

Theta_chi = mp.mpf('25') * R2
Theta_J = mp.mpf('25') * (R1**2)

print("Theta_w^(chi) / lambda_mu^2 =", Theta_chi)
print("Theta_w^(J)   / lambda_mu^2 =", Theta_J)

def expect_close(name: str, value, target, tol) -> None:
    diff = abs(value - target)
    print(f"{name} diff = {diff}")
    if diff > tol:
        raise AssertionError(f"{name} exceeds tol {tol}")

expect_close(
    "<rho>_chi",
    R1,
    mp.mpf('0.19261900555649309777068139356018510792903510747507'),
    mp.mpf('1e-28'),
)
expect_close(
    "<rho^2>_chi",
    R2,
    mp.mpf('0.16274529400326462037087418498629868328210821103971'),
    mp.mpf('1e-28'),
)
expect_close(
    "Theta_w^(chi)",
    Theta_chi,
    mp.mpf('4.0686323500816155092718546246574670820527052759928'),
    mp.mpf('1e-26'),
)
expect_close(
    "Theta_w^(J)",
    Theta_J,
    mp.mpf('0.92755203253930797183993260663904217023332624032789'),
    mp.mpf('1e-27'),
)
if not (Theta_chi >= Theta_J > 0):
    raise AssertionError("Expected Theta_w^(chi) >= Theta_w^(J) > 0")
