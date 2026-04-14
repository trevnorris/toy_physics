#!/usr/bin/env python3
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


def expect_zero(name: str, expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


banner("SAME-CHARGE BARRIER AUDIT — STAGE 024")
banner("EXACT PRIMITIVE RANKING ON THE SELECTED TWIN-SUPPORT BRANCH")

varrho, beta = sp.symbols("varrho beta", positive=True, real=True)

subbanner("I. Exact selected twin-support curve")
epsilon_star = sp.simplify(1 - sp.Rational(3, 2) * varrho)
sigma = sp.simplify(sp.Rational(4, 1) / (3 * varrho) - 2)
print(f"epsilon_* = {epsilon_star}")
print(f"sigma = {sigma}")

mixed_only_ceiling = sp.simplify(1 / varrho - 2)
twin_ceiling = sp.simplify(2 / varrho - 2)
print(f"mixed-only ceiling = {mixed_only_ceiling}")
print(f"symmetric-twin ceiling = {twin_ceiling}")

expect_zero(
    "selected branch above mixed-only lower edge",
    sp.simplify(sigma - mixed_only_ceiling - 1 / (3 * varrho)),
)
expect_zero(
    "selected branch below twin upper edge",
    sp.simplify(twin_ceiling - sigma - 2 / (3 * varrho)),
)

subbanner("II. Surviving threshold 1: q_W versus q_Lambda")
varrho_WL = sp.simplify(2 * (1 + beta**2) / (3 * (2 + beta**2)))
print(f"varrho_(WΛ) = {varrho_WL}")
# Recover from epsilon_* = 1/(2+beta^2)
expect_zero(
    "selected-curve inversion for varrho_(WΛ)",
    sp.simplify((1 - 3 * varrho_WL / 2) - 1 / (2 + beta**2)),
)

subbanner("III. Surviving threshold 2: |q_U| versus q_Lambda")
varrho_UL = sp.simplify(2 * (1 + beta**2) / (3 * (1 + beta + beta**2)))
print(f"varrho_(UΛ) = {varrho_UL}")
# Recover from epsilon_* = beta/(1+beta+beta^2)
expect_zero(
    "selected-curve inversion for varrho_(UΛ)",
    sp.simplify((1 - 3 * varrho_UL / 2) - beta / (1 + beta + beta**2)),
)

subbanner("IV. Threshold ordering")
ordering_gap = sp.simplify(varrho_UL - varrho_WL)
top_gap = sp.simplify(sp.Rational(2, 3) - varrho_UL)
print(f"varrho_(UΛ) - varrho_(WΛ) = {ordering_gap}")
print(f"2/3 - varrho_(UΛ) = {top_gap}")

expect_zero(
    "exact threshold ordering gap formula",
    ordering_gap
    - 2 * (1 + beta**2) * (1 - beta) / (3 * (1 + beta + beta**2) * (2 + beta**2)),
)
expect_zero(
    "upper selected-curve gap formula",
    top_gap - 2 * beta / (3 * (1 + beta + beta**2)),
)

subbanner("V. Exact endpoint values on the constructive coherent bound beta = 2/11")
beta_max = sp.Rational(2, 11)
varrho_WL_beta0 = sp.simplify(varrho_WL.subs(beta, 0))
varrho_WL_betamax = sp.simplify(varrho_WL.subs(beta, beta_max))
varrho_UL_beta0 = sp.simplify(varrho_UL.subs(beta, 0))
varrho_UL_betamax = sp.simplify(varrho_UL.subs(beta, beta_max))

print(f"varrho_(WΛ)(beta=0) = {varrho_WL_beta0}")
print(f"varrho_(WΛ)(beta=2/11) = {varrho_WL_betamax}")
print(f"varrho_(UΛ)(beta=0) = {varrho_UL_beta0}")
print(f"varrho_(UΛ)(beta=2/11) = {varrho_UL_betamax}")

expect_zero("varrho_(WΛ)(0) - 1/3", varrho_WL_beta0 - sp.Rational(1, 3))
expect_zero("varrho_(WΛ)(2/11) - 125/369", varrho_WL_betamax - sp.Rational(125, 369))
expect_zero("varrho_(UΛ)(0) - 2/3", varrho_UL_beta0 - sp.Rational(2, 3))
expect_zero("varrho_(UΛ)(2/11) - 250/441", varrho_UL_betamax - sp.Rational(250, 441))

subbanner("VI. Regime meaning")
print("Region I:  0 < varrho < varrho_(WΛ)      => q_chi > q_Z > q_Lambda > q_W > |q_U|")
print("Region II: varrho_(WΛ) < varrho < varrho_(UΛ) => q_chi > q_Z > q_W > q_Lambda > |q_U|")
print("Region III: varrho_(UΛ) < varrho < 2/3  => q_chi > q_Z > q_W > |q_U| > q_Lambda")

print("\nAll Stage 024 checks passed.")
