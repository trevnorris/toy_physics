#!/usr/bin/env python3
"""
moving_throat_pde_stage54_tanh_wall_branch_sympy_audit.py

SymPy audit for Stage 54:
canonical tanh-wall branch and natural local mouth closure.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 54 — CANONICAL TANH-WALL BRANCH")

xi = sp.symbols("xi", real=True)
a, L, ell = sp.symbols("a L ell", positive=True, real=True)
rho_w, c_sw, V0 = sp.symbols("rho_w c_sw V0", positive=True, real=True)
m, hbar = sp.symbols("m hbar", positive=True, real=True)
pi = sp.pi

f = sp.Rational(1, 2) * (1 + sp.tanh(xi))
fp = sp.simplify(sp.diff(f, xi))
fpp = sp.simplify(sp.diff(fp, xi))

If = sp.simplify(sp.integrate(fp**2, (xi, -sp.oo, sp.oo)))
t = sp.symbols("t", real=True)
If_sub = sp.simplify(sp.integrate(sp.Rational(1, 4) * (1 - t**2), (t, -1, 1)))
Ig = sp.simplify(sp.integrate(fpp**2, (xi, -sp.oo, sp.oo)))
Ig_sub = sp.simplify(sp.integrate(t**2 * (1 - t**2), (t, -1, 1)))  # tanh-substitution form

print("f'(xi) =", fp)
print("f''(xi) =", fpp)
print("I_f =", If)
print("I_f_sub =", If_sub)
print("I_g =", Ig)
print("I_g_sub =", Ig_sub)
expect_zero("I_f - 1/3", If - sp.Rational(1, 3))
expect_zero("I_f direct - substitution", If - If_sub)
expect_zero("I_g - 4/15", Ig - sp.Rational(4, 15))
expect_zero("I_g direct - substitution", Ig - Ig_sub)
expect_zero("I_g/I_f - 4/5", sp.simplify(Ig / If - sp.Rational(4, 5)))

Hw = sp.simplify(m * c_sw**2 / rho_w)
Tx = sp.simplify(sp.pi * a**2 * ell * If * hbar**2 / (m * rho_w))
Kx = sp.simplify(4 * sp.pi * a**2 * ell * If * Hw + sp.pi * a**2 * Ig * hbar**2 / (m * rho_w * ell))
J1 = sp.simplify(If / Hw)
Wwall = sp.simplify(4 * rho_w**2 * V0**2 * L**2 / (hbar**2 * c_sw**2 * ell**2))

print("T_X =", Tx)
print("K_X =", Kx)
print("J_1 =", J1)
print("W_wall =", Wwall)

Km = sp.simplify(Tx / ell)
eta = sp.simplify(Km * L / Tx)
print("K_m =", Km)
print("eta =", eta)
Km_expected = sp.pi * a**2 * hbar**2 / (3 * m * rho_w)
expect_zero("K_m - pi a^2 hbar^2 / (3 m rho_w)", Km - Km_expected)
expect_zero("eta - L/ell (from closed-form K_m)", (Km_expected * L / Tx) - L / ell)

chi_s, Lambda_ell, Upsilon_w = sp.symbols("chi_s Lambda_ell Upsilon_w", positive=True, real=True)
kappa = sp.simplify(Kx * L**2 / Tx)
kappa_red = sp.simplify(kappa.subs({
    m * c_sw * L / hbar: chi_s,
    L / ell: Lambda_ell
}))
print("kappa =", kappa)
print("kappa reduced =", kappa_red)
expect_zero("kappa_reduced - [4 chi_s^2 + 4/5 Lambda_ell^2]",
            kappa_red - (4 * chi_s**2 + sp.Rational(4, 5) * Lambda_ell**2))

W_red = sp.simplify(Wwall.subs({
    4 * rho_w**2 * V0**2 / (hbar**2 * c_sw**2): Upsilon_w,
    L / ell: Lambda_ell
}))
print("W_wall reduced =", W_red)
expect_zero("W_wall_reduced - Upsilon_w Lambda_ell^2", W_red - Upsilon_w * Lambda_ell**2)

banner("STAGE 54 THEOREM LEDGER")
print("I_f = 1/3")
print("I_g = 4/15")
print("I_g/I_f = 4/5")
print("eta = L/ell  under K_m = T_X/ell")
print("kappa = 4 chi_s^2 + (4/5) Lambda_ell^2")
print("W_wall = Upsilon_w Lambda_ell^2")
