#!/usr/bin/env python3
"""
moving_throat_pde_stage159_hybrid_outlet_projection_sympy_audit.py

SymPy-backed audit for Stage 159: linear projection of the co-evolving Family-1
defect onto the compensated Robin–mixed outlet.

Checks:
1. Linearized hybrid-outlet defects (delta E2, delta E4, Delta_Q) around the
   compensated canonical branch.
2. Mouth-gain to outlet-loading map and the exact cancellation of delta Sigma0
   in delta C = delta rho_R - 4 delta sigma_W.
3. Canonical-even preservation implies delta C = delta kappa_W = 0.
4. On that branch the last linear 2.5PN defect collapses to pure delta gamma_W.
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

banner("STAGE 142 — HYBRID OUTLET LINEARIZATION")

# General hybrid outlet data.
rho, sigma, kappa, gamma = sp.symbols("rho sigma kappa gamma", real=True)
sigma0 = sp.symbols("sigma0", real=True)
drho, dsigma, dkappa, dgamma = sp.symbols("drho dsigma dkappa dgamma", real=True)

L0 = -3 + rho - sigma
L2 = sp.Rational(1, 3) - sigma * kappa
L4 = sp.Rational(1, 9) - sigma * kappa**2
chi = sp.simplify(3 * (1 - 9 * sigma * gamma) / (3 - rho + sigma))
E2 = sp.simplify(-L2 / L0 - sp.Rational(1, 9))
E4 = sp.simplify(L2**2 / L0**2 - L4 / L0 - sp.Rational(4, 81))

subs = {
    rho: 4 * sigma0 + drho,
    sigma: sigma0 + dsigma,
    kappa: sp.Rational(1, 3) + dkappa,
    gamma: sp.Rational(1, 9) + dgamma,
}

def linearize(expr: sp.Expr, vars: list[sp.Symbol]) -> sp.Expr:
    expr = sp.expand(expr)
    for v in vars:
        expr = sp.series(expr, v, 0, 2).removeO()
    expr = sp.expand(expr)
    keep = 0
    poly = sp.Poly(expr, *vars)
    for monom, coeff in zip(poly.monoms(), poly.coeffs()):
        if sum(monom) <= 1:
            term = coeff
            for v, p in zip(vars, monom):
                term *= v**p
            keep += term
    return sp.expand(keep)

chi_lin = linearize(chi.subs(subs), [drho, dsigma, dkappa, dgamma])
E2_lin = linearize(E2.subs(subs), [drho, dsigma, dkappa, dgamma])
E4_lin = linearize(E4.subs(subs), [drho, dsigma, dkappa, dgamma])

delta_chi = sp.expand(chi_lin - 1)

delta_chi_expected = sp.expand((drho - 4 * dsigma - 27 * sigma0 * dgamma) / (3 * (1 - sigma0)))
E2_expected = sp.expand((drho - 4 * dsigma - 9 * sigma0 * dkappa) / (27 * (1 - sigma0)))
E4_expected = sp.expand((5 * drho - 20 * dsigma - 72 * sigma0 * dkappa) / (243 * (1 - sigma0)))

expect_zero("delta chi formula", delta_chi - delta_chi_expected)
expect_zero("delta E2 formula", E2_lin - E2_expected)
expect_zero("delta E4 formula", E4_lin - E4_expected)

banner("MOUTH-GAIN -> HYBRID-LOADING TRANSPORT")

Xi = sp.symbols("Xi", real=True)
Sigma0_can, dSigma0, dR = sp.symbols("Sigma0_can dSigma0 dR", real=True)
dMs = dSigma0
dMq = -sp.Rational(1, 4) * dSigma0 - Sigma0_can * dR
drho_expr = Xi * dMs
dsigma_expr = -Xi * dMq
deltaC_expr = sp.expand(drho_expr - 4 * dsigma_expr)
deltaC_expected = sp.expand(-4 * Xi * Sigma0_can * dR)
expect_zero("deltaC mouth transport", deltaC_expr - deltaC_expected)

sigma_star = sp.symbols("sigma_star", real=True)
expect_zero(
    "sigma_star substitution",
    deltaC_expected.subs(Xi, 4 * sigma_star / Sigma0_can) + 16 * sigma_star * dR
)

banner("CANONICAL-EVEN PRESERVATION")

deltaC, dk = sp.symbols("deltaC dk", real=True)
eq1 = sp.Eq(deltaC - 9 * sigma0 * dk, 0)
eq2 = sp.Eq(5 * deltaC - 72 * sigma0 * dk, 0)
sol = sp.solve([eq1, eq2], [deltaC, dk], dict=True)
print("solution =", sol)
if sol != [{deltaC: 0, dk: 0}]:
    raise AssertionError("Canonical-even preservation did not collapse to deltaC = dkappa = 0.")

det = sp.Matrix([[1, -9 * sigma0], [5, -72 * sigma0]]).det()
print("determinant =", sp.factor(det))

banner("FINAL REDUCED DEFECT")
final_defect = sp.simplify((deltaC - 27 * sigma0 * dgamma) / (3 * (1 - sigma0))).subs({deltaC: 0})
expect_zero("final Delta_Q + 9 sigma* dgamma /(1-sigma*)", final_defect + 9 * sigma0 * dgamma / (1 - sigma0))

print("\nCarry-forward formulas:")
print("  delta C := delta rho_R - 4 delta sigma_W")
print("  delta E2 = (delta C - 9 sigma_* delta kappa_W)/(27(1-sigma_*))")
print("  delta E4 = (5 delta C - 72 sigma_* delta kappa_W)/(243(1-sigma_*))")
print("  Delta_Q  = (delta C - 27 sigma_* delta gamma_W)/(3(1-sigma_*))")
print("  Canonical-even preservation => delta C = delta kappa_W = 0")
print("  Hence Delta_Q = -9 sigma_* delta gamma_W/(1-sigma_*)")
