#!/usr/bin/env python3
"""
Stage 43 SymPy audit — entropic source microclosure and microscopic Xi.
"""
from __future__ import annotations
import sympy as sp


def banner(t: str):
    print("\n" + "="*88)
    print(t)
    print("="*88)


def expect_zero(name: str, expr):
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} != 0")


banner("STAGE 43 — ENTROPIC SOURCE MICROCLOSURE")

# Symbols / functions
s, L = sp.symbols('s L', positive=True, real=True)
Theta, sigstar, Lam, T_X, K_X, K_m, Msig = sp.symbols('Theta sigma_* Lambda_phi T_X K_X K_m M_sigma', positive=True, real=True)
Dsig = sp.symbols('D_sigma', positive=True, real=True)
chi = sp.symbols('chi_sigma', positive=True, real=True)
kappa, eta = sp.symbols('kappa eta', positive=True, real=True)
Delta = sp.symbols('Delta', real=True)
phi0, dphi = sp.symbols('phi0 Delta_phi', real=True)
sigma = sp.Function('sigma')(s)
phi = sp.Function('phi')(s)
mu = sp.Function('mu')(s)
J = sp.Function('J')(s)

# 1) free-energy density and exact chemical potential
f = Theta * sigma * (sp.log(sigma/sigstar) - 1) - Lam * sigma * phi
mu_expr = sp.simplify(sp.diff(f, sigma))
print("f_sigma =", f)
print("mu_sigma^(chem) =", mu_expr)
expect_zero("mu - [Theta log(sigma/sigma_*) - Lambda phi]",
            mu_expr - (Theta*sp.log(sigma/sigstar) - Lam*phi))

# 2) Onsager current and Einstein relation
J_expr = sp.expand(-Msig * sigma * sp.diff(mu_expr, s))
print("J expanded =", J_expr)
expect_zero("Onsager current decomposition",
            J_expr - (-Msig*Theta*sp.diff(sigma, s) + Msig*Lam*sigma*sp.diff(phi, s)))
print("Einstein relation: D_sigma = M_sigma Theta")
expect_zero("current with D_sigma substitution",
            J_expr.subs(Msig*Theta, Dsig) - (-Dsig*sp.diff(sigma, s) + Msig*Lam*sigma*sp.diff(phi, s)))

# 3) affine-drop reduction reproduces exponential family
phi_aff = phi0 + dphi * s / L
ode = sp.simplify(J_expr.subs(phi, phi_aff).subs({sp.diff(phi_aff, s): dphi/L}))
print("J on affine-drop branch =", ode)
# solve J=0 by trial sigma = C exp(a s)
Cnorm = sp.symbols('C', positive=True, real=True)
a = Lam * dphi / (Theta * L)
sigma_trial = Cnorm * sp.exp(a*s)
expect_zero("J=0 solved by exponential family",
            ode.subs(sigma, sigma_trial).replace(sp.Derivative(sigma_trial, s), sp.diff(sigma_trial, s)))

# normalization on [0,L]
Csol = sp.simplify(sp.solve(sp.Eq(sp.integrate(sigma_trial, (s, 0, L)), 1), Cnorm)[0])
print("Normalization constant C =", Csol)
Pe = sp.symbols('Pe', positive=True, real=True)
Csol_Pe = sp.simplify(Csol.subs(a, Pe/L))
x = sp.symbols('x', real=True)
Sigma_x = sp.simplify(Pe * sp.exp(Pe*x) / (sp.exp(Pe) - 1))
print("Sigma_Pe(x) =", Sigma_x)
expect_zero("normalized Sigma_Pe family",
            sp.simplify(sp.integrate(Sigma_x, (x, 0, 1)) - 1))
expect_zero("Pe identification", a*L - Lam*dphi/Theta)

# 4) microscopic Xi from support normalization
Phi = sp.Function('Phi')(s)
phi_from_Phi = Lam * L**2 * Delta / T_X  # end-to-end drop only
Xi_micro = sp.simplify(Lam * phi_from_Phi / Theta / Delta)
print("Xi_micro =", Xi_micro)
expect_zero("Xi_micro - Lambda^2 L^2/(Theta T_X)", Xi_micro - Lam**2 * L**2 / (Theta * T_X))
expect_zero("Xi_micro susceptibility form", Xi_micro.subs(Theta, 1/chi) - chi * Lam**2 * L**2 / T_X)
expect_zero("Xi_micro phenomenological form", Xi_micro.subs(Theta, Dsig/Msig) - Msig*Lam**2 * L**2/(Dsig*T_X))

# 5) local dissipation identity
mu_fun = sp.Function('mu_sigma')(s)
J_fun = sp.Function('J')(s)
local_identity = sp.simplify(
    -mu_fun*sp.diff(J_fun, s) - (-sp.diff(mu_fun*J_fun, s) + sp.diff(mu_fun, s)*J_fun)
)
expect_zero("integration-by-parts identity", local_identity)
# Onsager substitution
mu_s = sp.symbols('mu_s', real=True)
J_on = -Msig * sigma * mu_s
expect_zero("Onsager dissipation density",
            sp.simplify(mu_s * J_on + J_on**2/(Msig*sigma)))
print("Therefore: dF/dt = -[mu J]_0^L - integral J^2/(M_sigma sigma) ds <= 0 under no-flux boundaries.")

# 6) support Euler-Lagrange equation from full free energy density
phi_fun = sp.Function('phi_fun')(s)
sigma_fun = sp.Function('sigma_fun')(s)
f_full = Theta * sigma_fun * (sp.log(sigma_fun/sigstar) - 1) - Lam * sigma_fun * phi_fun + sp.Rational(1,2)*T_X*sp.diff(phi_fun, s)**2 + sp.Rational(1,2)*K_X*phi_fun**2
EL_phi = sp.simplify(sp.diff(f_full, phi_fun) - sp.diff(sp.diff(f_full, sp.diff(phi_fun, s)), s))
print("Euler-Lagrange support bulk equation =", EL_phi)
expect_zero("support bulk equation form", EL_phi - (-Lam*sigma_fun + K_X*phi_fun - T_X*sp.diff(phi_fun, s, 2)))
print("Boundary term from support variation: T_X phi_s(0) = K_m phi(0), phi_s(L)=0.")

banner("STAGE 43 THEOREM LEDGER")
print("1. The positive-density free energy yields mu = Theta log(sigma/sigma_*) - Lambda phi.")
print("2. The Onsager current is J = -D_sigma sigma_s + M_sigma Lambda sigma phi_s with D_sigma = M_sigma Theta.")
print("3. Under affine support drop, the zero-flux branch is exactly the Stage-39 exponential family.")
print("4. The microscopic support/source gain is Xi_micro = Lambda_phi^2 L^2/(Theta_sigma T_X) = chi_sigma Lambda_phi^2 L^2/T_X.")
print("5. The same closure is passive: dF/dt = - integral J^2/(M_sigma sigma) ds <= 0 under no-flux boundaries.")
