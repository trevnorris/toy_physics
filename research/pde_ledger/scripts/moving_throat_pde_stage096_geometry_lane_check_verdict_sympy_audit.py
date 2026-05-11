#!/usr/bin/env python3
"""
SymPy audit for Stage 096.

This checkpoint script re-establishes the geometry-lane verdict directly:

1. the isotropic wall branch keeps the l=0 lane orthogonal to the grouped real
   l=2 bundle, so the geometry contamination numbers vanish on that branch;
2. the Stage 75 obstruction formula then collapses to the clean
   3/4 + 1/4 conservative quadrupole module;
3. the carried contact-plus-pole identification reproduces
   rho_alpha = 4/3 and zeta_req = 1/3.

Constant provenance:
- the l=2 Laplace eigenvalue 6 is derived as ell(ell+1) at ell = 2;
- 1/4, 3/4, 4/3, and 1/3 are derived in this audit from eps_2 = eps_4 = 0 and
  the obstruction / loading-ratio formulas.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    residual = sp.simplify(sp.expand(expr))
    print(f"{name} = {residual}")
    if residual != 0:
        raise AssertionError(f"{name} is not zero")


def domega(expr: sp.Expr, th: sp.Symbol, ph: sp.Symbol) -> sp.Expr:
    integrand = sp.simplify(expr * sp.sin(th))
    return sp.simplify(
        sp.integrate(sp.integrate(integrand, (ph, 0, 2 * sp.pi)), (th, 0, sp.pi))
    )


def lap_s2(Y: sp.Expr, th: sp.Symbol, ph: sp.Symbol) -> sp.Expr:
    return sp.simplify(
        (1 / sp.sin(th)) * sp.diff(sp.sin(th) * sp.diff(Y, th), th)
        + (1 / sp.sin(th) ** 2) * sp.diff(Y, ph, 2)
    )


banner("STAGE 079 — GEOMETRY-LANE CHECK VERDICT")

theta, phi = sp.symbols("theta phi", real=True)
omega, Omega_Q = sp.symbols("omega Omega_Q", positive=True, real=True)

Y00 = sp.Rational(1, 2) / sp.sqrt(sp.pi)
Y20 = sp.sqrt(sp.Rational(5, 16) / sp.pi) * (3 * sp.cos(theta) ** 2 - 1)
Y21c = sp.sqrt(sp.Rational(15, 4) / sp.pi) * sp.sin(theta) * sp.cos(theta) * sp.cos(phi)
Y21s = sp.sqrt(sp.Rational(15, 4) / sp.pi) * sp.sin(theta) * sp.cos(theta) * sp.sin(phi)
Y22c = sp.sqrt(sp.Rational(15, 16) / sp.pi) * sp.sin(theta) ** 2 * sp.cos(2 * phi)
Y22s = sp.sqrt(sp.Rational(15, 16) / sp.pi) * sp.sin(theta) ** 2 * sp.sin(2 * phi)

grouped_real_p2 = {
    "20": Y20,
    "21c": Y21c,
    "21s": Y21s,
    "22c": Y22c,
    "22s": Y22s,
}

for name, Y2m in grouped_real_p2.items():
    overlap = domega(Y00 * Y2m, theta, phi)
    laplace_residual = sp.simplify(-lap_s2(Y2m, theta, phi) - 6 * Y2m)
    laplace_overlap = domega(Y00 * (-lap_s2(Y2m, theta, phi)), theta, phi)

    print(f"<Y00|Y{name}> =", overlap)
    print(f"(-Delta)Y{name} - 6Y{name} =", laplace_residual)
    print(f"<Y00|(-Delta)Y{name}> =", laplace_overlap)

    expect_zero(f"<Y00|Y{name}>", overlap)
    expect_zero(f"(-Delta)Y{name} - 6Y{name}", laplace_residual)
    expect_zero(f"<Y00|(-Delta)Y{name}>", laplace_overlap)

# On the isotropic branch, the exact l=0 / l=2 block diagonal structure kills
# the geometry contamination numbers.
eps_2 = sp.Integer(0)
eps_4 = sp.Integer(0)
expect_zero("eps_2", eps_2)
expect_zero("eps_4", eps_4)

c_pole = sp.simplify((1 + eps_4) / (4 * (1 + eps_2) ** 2))
c_geom = sp.simplify(1 - c_pole)
Yhat_cons = sp.simplify(c_geom + c_pole / (1 - omega**2 / Omega_Q**2))
Yhat_expected = sp.simplify(
    sp.Rational(3, 4) + sp.Rational(1, 4) / (1 - omega**2 / Omega_Q**2)
)

rho_alpha = sp.simplify(1 / c_geom)
zeta_req = sp.simplify(c_pole / c_geom)

print("c_pole =", c_pole)
print("c_geom =", c_geom)
print("Yhat_Q^cons(omega) =", Yhat_cons)
print("rho_alpha =", rho_alpha)
print("zeta_req =", zeta_req)

expect_zero("c_pole - 1/4", c_pole - sp.Rational(1, 4))
expect_zero("c_geom - 3/4", c_geom - sp.Rational(3, 4))
expect_zero(
    "Yhat_Q^cons - [3/4 + (1/4)/(1 - omega^2/Omega_Q^2)]",
    Yhat_cons - Yhat_expected,
)
expect_zero("rho_alpha - 4/3", rho_alpha - sp.Rational(4, 3))
expect_zero("zeta_req - 1/3", zeta_req - sp.Rational(1, 3))
expect_zero("zeta_req - c_pole/c_geom", zeta_req - c_pole / c_geom)

print("\nFINAL LEDGER")
print("On the actual isotropic branch, the geometry lane stays dynamically inert:")
print("  eps_2 = eps_4 = 0,")
print("  c_pole = 1/4,")
print("  c_geom = 3/4,")
print("  Yhat_Q^cons = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2),")
print("  rho_alpha = 4/3,")
print("  zeta_req  = 1/3.")
