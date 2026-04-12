
#!/usr/bin/env python3
"""
5pn_stage90_general_isotropic_dtn_deformation.py

Stage 90 audit: general isotropic l=2 DtN deformation algebra.
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

banner("STAGE 90 — GENERAL ISOTROPIC l=2 DtN DEFORMATION ALGEBRA")

S, beta = sp.symbols("S beta", real=True)
Sigma0, Sigma2, Sigma4, Sigma5 = sp.symbols("Sigma_0 Sigma_2 Sigma_4 Sigma_5", real=True)
z = sp.symbols("z", real=True)

L0 = -3 * S + Sigma0
L2 = S * beta**2 / 3 + Sigma2
L4 = S * beta**4 / 9 + Sigma4
L5 = S * beta**5 / 9 + Sigma5

Ydef = sp.simplify(L0 / (L0 + L2 * z**2 + L4 * z**4 + sp.I * L5 * z**5))
Ydef_series = sp.series(Ydef, z, 0, 6).removeO()
print("Yhat_2^def(z) =")
sp.pprint(Ydef_series)

coef2 = sp.simplify(sp.expand(Ydef_series).coeff(z, 2))
coef4 = sp.simplify(sp.expand(Ydef_series).coeff(z, 4))
coef5 = sp.simplify(sp.expand(Ydef_series).coeff(z, 5) / sp.I)

expect_zero("z^2 coefficient + L2/L0", coef2 + L2 / L0)
expect_zero("z^4 coefficient - (L2^2/L0^2 - L4/L0)", coef4 - (L2**2 / L0**2 - L4 / L0))
expect_zero("z^5 coefficient + L5/L0", coef5 + L5 / L0)

sol_even = sp.solve(
    [
        sp.Eq(-L2 / L0, sp.Rational(1, 9)),
        sp.Eq(L2**2 / L0**2 - L4 / L0, sp.Rational(4, 81)),
    ],
    [Sigma2, Sigma4],
    dict=True,
)[0]
print("Sigma_2 =", sp.factor(sol_even[Sigma2]))
print("Sigma_4 =", sp.factor(sol_even[Sigma4]))

expect_zero(
    "Sigma_2 formula",
    sol_even[Sigma2] + (3 * S * beta**2 - 3 * S + Sigma0) / 9,
)
expect_zero(
    "Sigma_4 formula",
    sol_even[Sigma4] + (3 * S * beta**4 - 3 * S + Sigma0) / 27,
)

chi_Q = sp.simplify((-L5 / L0) / sp.Rational(1, 27))
chi_Q_even = sp.simplify(chi_Q.subs(sol_even))
print("chi_Q =", sp.factor(chi_Q_even))
expect_zero(
    "chi_Q - 3(S beta^5 + 9 Sigma_5)/(3S - Sigma_0)",
    chi_Q_even - 3 * (S * beta**5 + 9 * Sigma5) / (3 * S - Sigma0),
)
expect_zero(
    "chi_Q - 1 formula",
    sp.simplify(chi_Q_even - 1)
    - (3 * S * (beta**5 - 1) + Sigma0 + 27 * Sigma5) / (3 * S - Sigma0),
)

banner("STAGE 90 FINAL LEDGER")
print("The first explicit isotropic moving-throat DtN deformation family gives an exact map")
print("from throat-core deformation data to the retarded quadrupole normalization scalar:")
print("  chi_Q = 3(S beta^5 + 9 Sigma_5)/(3S - Sigma_0).")
print("So the only isotropic branch data that can move the canonical value chi_Q = 1 are")
print("  beta, Sigma_0, Sigma_5,")
print("with Sigma_2 and Sigma_4 fixed by the requirement that the canonical even moments stay unchanged.")
