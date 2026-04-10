#!/usr/bin/env python3
"""
5pn_stage93_robin_outlet_model.py

Stage 93 audit: explicit isotropic Robin outlet model.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr) -> None:
    expr_s = sp.simplify(sp.together(sp.expand(expr)))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 93 — EXPLICIT ISOTROPIC ROBIN OUTLET MODEL")

rho_R, z = sp.symbols("rho_R z", real=True)
Lambda_out = -3 + z**2 / 3 + z**4 / 9 + sp.I * z**5 / 9
Lambda_R = sp.expand(Lambda_out + rho_R)
Y_R = sp.simplify((rho_R - 3) / Lambda_R)
Y_R_series = sp.series(Y_R, z, 0, 6).removeO()

print("Lambda_2^R(z) =")
sp.pprint(Lambda_R)
print("Yhat_2^R(z) =")
sp.pprint(Y_R_series)

coef2 = sp.simplify(sp.expand(Y_R_series).coeff(z, 2))
coef4 = sp.simplify(sp.expand(Y_R_series).coeff(z, 4))
coef5 = sp.simplify(sp.expand(Y_R_series).coeff(z, 5) / sp.I)

expect_zero("z^2 coefficient", coef2 - 1 / (9 - 3 * rho_R))
expect_zero("z^4 coefficient", coef4 - (4 - rho_R) / (9 * (3 - rho_R)**2))
expect_zero("z^5 coefficient", coef5 - 1 / (27 - 9 * rho_R))

chi_R = sp.simplify(coef5 / sp.Rational(1, 27))
print("chi_Q^R =", sp.factor(chi_R))
expect_zero("chi_Q^R formula", chi_R - 3 / (3 - rho_R))

chi_R_series = sp.series(chi_R, rho_R, 0, 3).removeO()
print("small-rho_R expansion =", chi_R_series)
expect_zero(
    "small-rho_R expansion check",
    chi_R_series - (1 + rho_R / 3 + rho_R**2 / 9),
)

banner("STAGE 93 FINAL LEDGER")
print("The raw isotropic Robin core deforms both the even fingerprint and the odd normalization.")
print("Its exact normalization factor is chi_Q^R = 3/(3-rho_R), so the linearized branch-selection triple is")
print("  (b, a_0, a_5) = (0, rho_R, 0).")
