#!/usr/bin/env python3
"""
5pn_stage91_robustness_classes.py

Stage 91 audit: robustness classes for chi_Q.
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

banner("STAGE 91 — ROBUSTNESS CLASSES FOR chi_Q")

S, beta = sp.symbols("S beta", real=True)
Sigma0, Sigma2, Sigma4, Sigma5 = sp.symbols("Sigma_0 Sigma_2 Sigma_4 Sigma_5", real=True)
z = sp.symbols("z", real=True)

# Stage-90 general isotropic DtN deformation data.
chi_Q = sp.simplify(3 * (S * beta**5 + 9 * Sigma5) / (3 * S - Sigma0))
print("chi_Q =", sp.factor(chi_Q))

# Canonical outgoing fingerprint.
Y_out = 1 + z**2 / 9 + 4 * z**4 / 81 + sp.I * z**5 / 27

# ------------------------------------------------------------------
# Class A: pure scale deformation
# ------------------------------------------------------------------
banner("CLASS A — PURE SCALE DEFORMATION")

Y_A = sp.simplify((-3 * S) / (S * (-3 + z**2 / 3 + z**4 / 9 + sp.I * z**5 / 9)))
Y_A_series = sp.series(Y_A, z, 0, 6).removeO()
print("Yhat_A(z) =")
sp.pprint(Y_A_series)
expect_zero("Class-A normalized response minus canonical", Y_A_series - Y_out)
expect_zero("Class-A chi_Q - 1", chi_Q.subs({beta: 1, Sigma0: 0, Sigma5: 0}) - 1)

# ------------------------------------------------------------------
# Class B: pure scale + argument deformation
# ------------------------------------------------------------------
banner("CLASS B — PURE SCALE + ARGUMENT DEFORMATION")

Y_B = sp.simplify((-3 * S) / (S * (-3 + beta**2 * z**2 / 3 + beta**4 * z**4 / 9 + sp.I * beta**5 * z**5 / 9)))
Y_B_series = sp.series(Y_B, z, 0, 6).removeO()
print("Yhat_B(z) =")
sp.pprint(Y_B_series)

coef2 = sp.simplify(sp.expand(Y_B_series).coeff(z, 2))
coef4 = sp.simplify(sp.expand(Y_B_series).coeff(z, 4))
coef5 = sp.simplify(sp.expand(Y_B_series).coeff(z, 5) / sp.I)
print("z^2 coefficient =", coef2)
print("z^4 coefficient =", coef4)
print("z^5 coefficient =", coef5)

# Keeping the even fingerprint canonical forces beta = 1 on the positive branch.
sol_beta = sp.solve([sp.Eq(coef2, sp.Rational(1, 9)), sp.Eq(coef4, sp.Rational(4, 81))], [beta], dict=True)
print("even-preserving beta branches =", sol_beta)
expect_zero("beta=1 preserves Class-B response", Y_B_series.subs(beta, 1) - Y_out)
expect_zero("Class-B chi_Q(beta=1) - 1", chi_Q.subs({beta: 1, Sigma0: 0, Sigma5: 0}) - 1)

# ------------------------------------------------------------------
# Class C: additive isotropic throat-core channel
# ------------------------------------------------------------------
banner("CLASS C — ADDITIVE ISOTROPIC THROAT-CORE CHANNEL")

# Hold beta = 1 but allow additive core deformations.
chi_C = sp.simplify(chi_Q.subs(beta, 1))
# Even-preserving slots from Stage 90 at beta=1.
Sigma2_even = sp.simplify(-Sigma0 / 9)
Sigma4_even = sp.simplify(-Sigma0 / 27)
print("Sigma_2 =", Sigma2_even)
print("Sigma_4 =", Sigma4_even)

chi_C_even = sp.simplify(chi_C)
print("chi_Q^(Class C) =", sp.factor(chi_C_even))
expect_zero(
    "Class-C chi_Q formula",
    chi_C_even - 3 * (S + 9 * Sigma5) / (3 * S - Sigma0),
)

chi_C_pure_even = sp.simplify(chi_C_even.subs(Sigma5, 0))
print("pure-even additive chi_Q =", sp.factor(chi_C_pure_even))
expect_zero(
    "pure-even additive chi_Q formula",
    chi_C_pure_even - 3 * S / (3 * S - Sigma0),
)

# Exact preservation submanifold.
Sigma5_pres = sp.solve(sp.Eq(chi_Q, 1), Sigma5)[0]
print("Sigma_5 preservation law =", sp.factor(Sigma5_pres))
expect_zero(
    "preservation submanifold",
    Sigma5_pres - (S * (1 - beta**5) / 9 - Sigma0 / 27),
)

banner("STAGE 91 FINAL LEDGER")
print("Class A: pure overall mouth normalization is invisible to the normalized outgoing fingerprint.")
print("Class B: pure scale+argument deformation preserves the canonical even branch only on the natural positive branch beta=1, hence chi_Q=1.")
print("Class C: a genuine additive isotropic throat-core channel can move chi_Q while leaving the lower even moments canonical.")
print("Exact preservation submanifold:")
print("  Sigma_5 = S(1 - beta^5)/9 - Sigma_0/27.")
