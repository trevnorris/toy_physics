#!/usr/bin/env python3
"""
g2_step54_explicit_outlet_deformation_targets.py

Step 54 of the g-2 chain.

What this script does
---------------------
1. Builds the exact isotropic DtN deformation law around the canonical outgoing l=2 branch.
2. Imposes the canonical even fingerprint and extracts chi_Q as a function of the
   deformation data (S, beta, Sigma0, Sigma5).
3. Solves the electron target chi_e = 1/(1+delta) in a few explicit outlet classes:
   - pure additive / Robin-like core,
   - compensated Robin-mixed outlet,
   - linearized tangent branch.
"""

from __future__ import annotations

import sympy as sp

I = sp.I


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


banner("STEP 54 — EXPLICIT OUTLET DEFORMATION TARGETS")

# ---------------------------------------------------------------------------
# 1. General isotropic DtN deformation preserving the canonical even fingerprint
# ---------------------------------------------------------------------------

z = sp.symbols("z", real=True)
S, beta = sp.symbols("S beta", positive=True, real=True)
Sigma0, Sigma2, Sigma4, Sigma5 = sp.symbols("Sigma0 Sigma2 Sigma4 Sigma5", real=True)

Lambda_out = -3 + z**2/sp.Integer(3) + z**4/sp.Integer(9) + I*z**5/sp.Integer(9)
Lambda_def = sp.expand(S * Lambda_out.subs(z, beta*z) + Sigma0 + Sigma2*z**2 + Sigma4*z**4 + I*Sigma5*z**5)

L0 = sp.expand(Lambda_def).coeff(z, 0)
L2 = sp.expand(Lambda_def).coeff(z, 2)
L4 = sp.expand(Lambda_def).coeff(z, 4)
L5 = sp.simplify(sp.expand(Lambda_def).coeff(z, 5) / I)

print("L0 =", L0)
print("L2 =", L2)
print("L4 =", L4)
print("L5 =", L5)

# Canonical even fingerprint.
sol_even = sp.solve(
    [
        sp.Eq(-L2 / L0, sp.Rational(1, 9)),
        sp.Eq(L2**2 / L0**2 - L4 / L0, sp.Rational(4, 81)),
    ],
    [Sigma2, Sigma4],
    simplify=True,
    dict=True,
)[0]

Sigma2_sol = sp.simplify(sol_even[Sigma2])
Sigma4_sol = sp.simplify(sol_even[Sigma4])

print("Sigma2 =", Sigma2_sol)
print("Sigma4 =", Sigma4_sol)

chi_Q = sp.simplify((-L5 / L0) / sp.Rational(1, 27))
chi_Q_even = sp.simplify(chi_Q.subs({Sigma2: Sigma2_sol, Sigma4: Sigma4_sol}))
print("chi_Q =", chi_Q_even)

chi_expected = sp.simplify(3 * (S * beta**5 + 9 * Sigma5) / (3 * S - Sigma0))
expect_zero("chi_Q - expected", chi_Q_even - chi_expected)

# ---------------------------------------------------------------------------
# 2. Electron target branch
# ---------------------------------------------------------------------------

delta = sp.symbols("delta", positive=True, real=True)
chi_e = sp.simplify(1 / (1 + delta))
print("chi_e =", chi_e)

# Pure additive / Robin-like core: beta=1, Sigma5=0.
Sigma0_add = sp.solve(sp.Eq(chi_Q_even.subs({beta: 1, Sigma5: 0}), chi_e), Sigma0)[0]
print("Sigma0 (pure additive core) =", sp.simplify(Sigma0_add))

# Pure Robin outlet chi_Q^R = 3/(3-rho_R)
rho_R = sp.symbols("rho_R", real=True)
chi_R = sp.simplify(3 / (3 - rho_R))
rho_R_sol = sp.solve(sp.Eq(chi_R, chi_e), rho_R)[0]
print("rho_R (pure Robin outlet) =", sp.simplify(rho_R_sol))

# Compensated Robin-mixed outlet.
sigma_W, gamma_W = sp.symbols("sigma_W gamma_W", positive=True, real=True)
chi_hyb = sp.simplify((1 - 9 * sigma_W * gamma_W) / (1 - sigma_W))
gamma_hyb_sol = sp.solve(sp.Eq(chi_hyb, chi_e), gamma_W)[0]
print("gamma_W (compensated hybrid outlet) =", sp.simplify(gamma_hyb_sol))

# ---------------------------------------------------------------------------
# 3. Linearized branch-selection triple
# ---------------------------------------------------------------------------

eps = sp.symbols("eps", real=True)
s, b, a0, a5 = sp.symbols("s b a0 a5", real=True)

chi_lin = sp.simplify(
    chi_Q_even.subs(
        {
            S: 1 + eps * s,
            beta: 1 + eps * b,
            Sigma0: eps * a0,
            Sigma5: eps * a5,
        }
    )
)
chi_lin_series = sp.expand(sp.series(chi_lin, eps, 0, 2).removeO())
print("chi_lin =", chi_lin_series)

coeff_eps = sp.simplify(sp.expand(chi_lin_series).coeff(eps, 1))
print("linear coefficient =", coeff_eps)
expect_zero("linear coefficient - expected", coeff_eps - (5*b + a0/sp.Integer(3) + 9*a5))

banner("FINAL LEDGER")
print("Exact isotropic outgoing deformation law:")
print("  chi_Q = 3 (S beta^5 + 9 Sigma5) / (3S - Sigma0)")
print("Electron target:")
print("  chi_e = 1 / (1 + delta)")
print("Representative exact realizations:")
print("  pure additive core : Sigma0 = {}".format(sp.simplify(Sigma0_add)))
print("  pure Robin outlet  : rho_R  = {}".format(sp.simplify(rho_R_sol)))
print("  compensated hybrid : gamma_W = {}".format(sp.simplify(gamma_hyb_sol)))
print("Linearized tangent law:")
print("  delta chi_Q = eps * (5 b + a0/3 + 9 a5) + O(eps^2)")
