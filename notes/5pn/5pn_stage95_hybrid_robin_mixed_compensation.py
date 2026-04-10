#!/usr/bin/env python3
"""
5pn_stage95_hybrid_robin_mixed_compensation.py

Stage 95 audit: exact Robin-mixed compensation law.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr) -> None:
    expr_s = sp.simplify(sp.factor(sp.together(sp.expand(expr))))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 95 — EXACT ROBIN–MIXED COMPENSATION LAW")

rho_R, sigma_W, kappa_W, gamma_W, z = sp.symbols("rho_R sigma_W kappa_W gamma_W z", real=True)

Lambda_out = -3 + z**2 / 3 + z**4 / 9 + sp.I * z**5 / 9
hidden_inv = sp.series(1 / (1 - kappa_W * z**2 - sp.I * gamma_W * z**5), z, 0, 6).removeO()
Lambda_hyb = sp.expand(Lambda_out + rho_R - sigma_W * hidden_inv)

print("Lambda_2^hyb(z) =")
sp.pprint(Lambda_hyb)

L0 = sp.expand(Lambda_hyb).coeff(z, 0)
L2 = sp.expand(Lambda_hyb).coeff(z, 2)
L4 = sp.expand(Lambda_hyb).coeff(z, 4)
L5 = sp.simplify(sp.expand(Lambda_hyb).coeff(z, 5) / sp.I)

expect_zero("L0 - (-3 + rho_R - sigma_W)", L0 - (-3 + rho_R - sigma_W))
expect_zero("L2 - (1/3 - sigma_W kappa_W)", L2 - (sp.Rational(1, 3) - sigma_W * kappa_W))
expect_zero("L4 - (1/9 - sigma_W kappa_W^2)", L4 - (sp.Rational(1, 9) - sigma_W * kappa_W**2))
expect_zero("L5 - (1/9 - sigma_W gamma_W)", L5 - (sp.Rational(1, 9) - sigma_W * gamma_W))

eq2 = sp.simplify(-L2 / L0 - sp.Rational(1, 9))
eq4 = sp.simplify(L2**2 / L0**2 - L4 / L0 - sp.Rational(4, 81))

sol_branches = sp.solve([sp.Eq(eq2, 0), sp.Eq(eq4, 0)], [rho_R, kappa_W], dict=True)
print("canonical-even branches =", sol_branches)

branch_trivial = {rho_R: sigma_W, kappa_W: 0}
branch_nontriv = {rho_R: 4 * sigma_W, kappa_W: sp.Rational(1, 3)}

# Verify both branches.
expect_zero("trivial branch eq2", eq2.subs(branch_trivial))
expect_zero("trivial branch eq4", eq4.subs(branch_trivial))
expect_zero("nontrivial branch eq2", eq2.subs(branch_nontriv))
expect_zero("nontrivial branch eq4", eq4.subs(branch_nontriv))

Lambda_nontriv = sp.expand(Lambda_hyb.subs(branch_nontriv))
print("nontrivial compensated branch =")
sp.pprint(Lambda_nontriv)

expect_zero(
    "nontrivial branch form",
    Lambda_nontriv
    - (-3 * (1 - sigma_W) + (1 - sigma_W) * z**2 / 3 + (1 - sigma_W) * z**4 / 9
       + sp.I * (sp.Rational(1, 9) - sigma_W * gamma_W) * z**5),
)

chi_hyb = sp.simplify((-(sp.expand(Lambda_nontriv).coeff(z, 5) / sp.I) / sp.expand(Lambda_nontriv).coeff(z, 0)) / sp.Rational(1, 27))
print("chi_Q^hyb =", sp.factor(chi_hyb))
expect_zero("chi_Q^hyb formula", chi_hyb - (1 - 9 * sigma_W * gamma_W) / (1 - sigma_W))

gamma_pres = sp.solve(sp.Eq(chi_hyb, 1), gamma_W)[0]
print("gamma_W preserving canonical normalization =", gamma_pres)
expect_zero("gamma_W preservation value", gamma_pres - sp.Rational(1, 9))

Lambda_scale = sp.expand(Lambda_nontriv.subs(gamma_W, sp.Rational(1, 9)))
expect_zero("collapse to pure-scale class", Lambda_scale - (1 - sigma_W) * Lambda_out)

chi_hyb_series = sp.series(chi_hyb, sigma_W, 0, 2).removeO()
print("small-sigma_W expansion =", chi_hyb_series)
expect_zero("linearized preservation law", chi_hyb_series - (1 + sigma_W * (1 - 9 * gamma_W)))

banner("STAGE 95 FINAL LEDGER")
print("The hybrid outlet admits exactly two canonical-even branches.")
print("The nontrivial compensated branch is")
print("  rho_R = 4 sigma_W,   kappa_W = 1/3,")
print("with outgoing normalization")
print("  chi_Q^hyb = (1 - 9 sigma_W gamma_W)/(1 - sigma_W).")
print("It preserves the canonical outgoing branch exactly when gamma_W = 1/9, in which case the whole outlet collapses to a harmless pure-scale deformation.")
