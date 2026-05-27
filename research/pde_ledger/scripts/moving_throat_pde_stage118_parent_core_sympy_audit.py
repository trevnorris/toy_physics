#!/usr/bin/env python3
"""
moving_throat_pde_stage118_parent_core_sympy_audit.py

SymPy audit for the first parent-action extraction of
(K_s, K_q, lambda, g_s, g_q) from an explicit GNLS + localized-Maxwell throat-core ansatz.
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

banner("I. Tanh-wall shell moments")

y = sp.symbols("y", real=True)
chi_s = sp.Rational(1,2) * sp.sech(y)**2
t = sp.symbols("t", real=True)
# Use t = tanh(y), dt = sech(y)^2 dy, to avoid special-function ambiguity.
I_f = sp.simplify(sp.integrate((sp.Rational(1,4))*(1 - t**2), (t, -1, 1)))
I_g = sp.simplify(sp.integrate(t**2*(1 - t**2), (t, -1, 1)))
print("I_f =", I_f)
print("I_g =", I_g)
if sp.simplify(I_f - sp.Rational(1,3)) != 0 or sp.simplify(I_g - sp.Rational(4,15)) != 0:
    raise AssertionError("Unexpected tanh-wall moments.")

banner("II. D/N half-wave normalization data")

z, L_W = sp.symbols("z L_W", positive=True, real=True)
chi = sp.sqrt(2/L_W) * sp.sin(sp.pi*z/(2*L_W))
chi_norm = sp.simplify(sp.integrate(chi**2, (z, 0, L_W)))
chi_grad = sp.simplify(sp.integrate(sp.diff(chi,z)**2, (z, 0, L_W)))
chi_int = sp.simplify(sp.integrate(chi, (z, 0, L_W)))
chi_prime_0 = sp.simplify(sp.diff(chi,z).subs(z,0))
print("∫ chi^2 dz =", chi_norm)
print("∫ (chi')^2 dz =", chi_grad)
print("∫ chi dz =", chi_int)
print("chi'(0) =", chi_prime_0)
expect_zero("D/N norm check", chi_norm - 1)
expect_zero("D/N stiffness check", chi_grad - sp.pi**2/(4*L_W**2))

banner("III. Shell stiffness K_s")

a, ell, H_w, hbar, mpsi, rho_w = sp.symbols("a ell H_w hbar mpsi rho_w", positive=True, real=True)
K_s = 4*sp.pi*a**2*(H_w*ell*I_f + (hbar**2/(4*mpsi*rho_w))*(I_g/ell))
K_s_expected = 4*sp.pi*a**2*(H_w*ell/sp.Integer(3) + hbar**2/(15*mpsi*rho_w*ell))
print("K_s =", sp.simplify(K_s))
expect_zero("K_s closed form", K_s - K_s_expected)

c_sw = sp.symbols("c_sw", positive=True, real=True)
healing_sub = {H_w: mpsi*c_sw**2/rho_w, ell: hbar/(2*mpsi*c_sw)}
K_s_heal = sp.simplify(K_s.subs(healing_sub))
print("K_s on healing lock =", K_s_heal)
expect_zero(
    "healing-lock K_s",
    K_s_heal - 3*sp.pi*a**2*hbar**2/(5*mpsi*rho_w*(hbar/(2*mpsi*c_sw)))
)

banner("IV. Bilinear GNLS shell/mixed hybridization")

rho0, varrho_s, v0, qstar, A_q, s, q, m = sp.symbols(
    "rho0 varrho_s v0 qstar A_q s q m", real=True
)
expr = sp.expand(sp.Rational(1,2)*m*(rho0 + s*varrho_s)*(v0 - (qstar/m)*q*A_q)**2)
sq_coeff = sp.expand(expr).coeff(s,1).coeff(q,1)
print("sq coefficient =", sq_coeff)
expect_zero("bilinear sq coefficient", sq_coeff + qstar*varrho_s*v0*A_q)

banner("V. Effective mixed stiffness and mouth couplings")

mu0, Zq, c_s, Tm = sp.symbols("mu0 Zq c_s Tm", positive=True, real=True)
K_q = sp.simplify((Zq/mu0) * (sp.pi**2*c_s**2/(4*L_W**2)))
g_q = sp.simplify((Zq/mu0) * chi_prime_0)
J_s = sp.simplify(4*sp.pi*a**2*ell*I_f)
g_s = sp.simplify(Tm * J_s)
I_q = chi_int
I_sq_uniform = sp.simplify(J_s * I_q)
lam_uniform = sp.simplify(-qstar * v0 * I_sq_uniform)

print("K_q =", K_q)
print("g_q =", g_q)
print("J_s =", J_s)
print("g_s =", g_s)
print("I_sq (uniform-core closure) =", I_sq_uniform)
print("lambda (uniform-core closure) =", lam_uniform)

expect_zero("K_q closed form", K_q - (Zq/mu0) * sp.pi**2 * c_s**2 / (4*L_W**2))
expect_zero("g_q closed form", g_q - (Zq/mu0) * sp.pi/(sp.sqrt(2)*L_W**sp.Rational(3,2)))
expect_zero("J_s closed form", J_s - 4*sp.pi*a**2*ell/3)
expect_zero("g_s closed form", g_s - Tm * 4*sp.pi*a**2*ell/3)
expect_zero("I_q closed form", I_q - 2*sp.sqrt(2*L_W)/sp.pi)
expect_zero("lambda uniform closure",
            lam_uniform + 8*sp.sqrt(2)*qstar*v0*a**2*ell*sp.sqrt(L_W)/3)
expect_zero("lambda from bilinear", lam_uniform + qstar * v0 * J_s * I_q)

print("\nAll Stage 118 symbolic checks passed.")
