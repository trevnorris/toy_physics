#!/usr/bin/env python3
"""
5pn_stage43_entropic_microclosure.py

Stage 43 audit: entropic source microclosure and the microscopic support/source gain.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr_s = sp.simplify(sp.together(sp.expand(expr)))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")


banner("STAGE 43 — ENTROPIC SOURCE MICROCLOSURE AND THE MICROSCOPIC GAIN")

s, L = sp.symbols("s L", positive=True, real=True)
Theta_sigma, sigma_star = sp.symbols("Theta_sigma sigma_star", positive=True, real=True)
Lambda_phi, T_X, K_X, K_m = sp.symbols("Lambda_phi T_X K_X K_m", positive=True, real=True)
M_sigma, D_sigma = sp.symbols("M_sigma D_sigma", positive=True, real=True)
Delta_phi = sp.symbols("Delta_phi", real=True)
Delta = sp.symbols("Delta", real=True)

sigma = sp.Function("sigma")(s)
phi = sp.Function("phi")(s)
mu = Theta_sigma * sp.log(sigma / sigma_star) - Lambda_phi * phi

subbanner("43.1 — Exact free-energy variations")
print("mu_sigma^(chem) =")
sp.pprint(mu)

phi_density = T_X * sp.diff(phi, s)**2 / 2 + K_X * phi**2 / 2 - Lambda_phi * sigma * phi
phi_euler = sp.simplify(sp.diff(phi_density, phi) - sp.diff(sp.diff(phi_density, sp.diff(phi, s)), s))
expect_zero("support Euler–Lagrange equation", phi_euler - (K_X * phi - Lambda_phi * sigma - T_X * sp.diff(phi, s, 2)))
print("Robin/Neumann boundary conditions: T_X phi_s(0)=K_m phi(0),  phi_s(L)=0")

subbanner("43.2 — Onsager current and Einstein relation")
J = -M_sigma * sigma * sp.diff(mu, s)
J_expanded = sp.expand(J.doit())
expect_zero(
    "drift-diffusion current form",
    J_expanded - (-M_sigma * Theta_sigma * sp.diff(sigma, s) + M_sigma * Lambda_phi * sigma * sp.diff(phi, s)),
)
print("Thus D_sigma = M_sigma Theta_sigma.")

subbanner("43.3 — Stationary zero-flux branch")
ode_expr = sp.simplify(sp.diff(sigma, s) - (Lambda_phi * Delta_phi) / (Theta_sigma * L) * sigma)
print("J=0 with affine phi gives ODE:")
sp.pprint(ode_expr)
print("Solution family: sigma(s) = C * exp[(Lambda_phi Delta_phi)/(Theta_sigma L) * s]")
print("Normalized on x=s/L, this becomes Sigma_Pe(x) = Pe exp(Pe x)/(exp(Pe)-1)")
print("with exact Peclet number  Pe = (Lambda_phi/Theta_sigma) Delta_phi.")

subbanner("43.4 — Exact microscopic coupling Xi_micro")
Delta_phi_support = Lambda_phi * L**2 * Delta / T_X
Pe_from_support = sp.simplify((Lambda_phi / Theta_sigma) * Delta_phi_support)
Xi_formula = sp.simplify(Pe_from_support / Delta)
expect_zero("Xi_micro formula", Xi_formula - Lambda_phi**2 * L**2 / (Theta_sigma * T_X))
expect_zero(
    "phenomenological/microscopic Xi match after Einstein relation",
    Xi_formula.subs(Theta_sigma, D_sigma / M_sigma) - M_sigma * Lambda_phi**2 * L**2 / (D_sigma * T_X),
)

subbanner("43.5 — Exact passivity identity")
mu_s = sp.diff(mu, s)
expect_zero("local entropy-production identity", mu_s * J_expanded + J_expanded**2 / (M_sigma * sigma))
print("So under no-flux boundaries, dF/dt = -∫ ds J^2/(M_sigma sigma) <= 0.")

banner("STAGE 43 FINAL LEDGER")
print("On the first explicit positive-density Onsager closure, the support/source strength is")
print("  Xi_micro = Lambda_phi^2 L^2 / (Theta_sigma T_X) = chi_sigma Lambda_phi^2 L^2 / T_X,")
print("and the source branch is automatically passive.")
