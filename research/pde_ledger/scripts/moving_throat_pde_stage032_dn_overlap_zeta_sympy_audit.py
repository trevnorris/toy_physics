#!/usr/bin/env python3
"""
Stage 32 SymPy audit.

Checks:
1. Exact D/N half-wave momentum and stiffness data.
2. Exact uniform-source overlap integral on the D/N ladder.
3. Exact ratio I_n / I_0 = 1 / (2n+1).
4. Exact microscopic coherent-support law zeta_n^(phys).
5. Same-operator twin-family stiffness relation.
6. Exact twin-lane support ratio zeta_n^(twin).
7. Lowest twin half-wave value zeta_0^(twin) = 1.
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


def dn_halfwave_momentum(n: sp.Expr, L: sp.Expr) -> sp.Expr:
    return sp.simplify((n + sp.Rational(1, 2)) * sp.pi / L)


def uniform_dn_overlap(n: sp.Expr, L: sp.Expr) -> sp.Expr:
    return sp.simplify(sp.sqrt(2 * L) / ((n + sp.Rational(1, 2)) * sp.pi))


def overlap_ratio(n: sp.Expr) -> sp.Expr:
    return sp.simplify(sp.Rational(1, 1) / (2 * n + 1))


def twin_support_ratio(n: sp.Expr, x: sp.Expr) -> sp.Expr:
    return sp.simplify(1 / ((2 * n + 1) ** 2 * (1 + x * n * (n + 1))))


def main() -> None:
    banner("STAGE 32 — EXPLICIT D/N OVERLAP EXTRACTION OF THE PHYSICAL SUPPORT RATIO")

    n = sp.symbols("n", integer=True, nonnegative=True)
    s = sp.symbols("s", real=True)
    L = sp.symbols("L", positive=True, real=True)
    lambda_star = sp.symbols("lambda_star", nonzero=True, real=True)
    KW_eff, Kphi_eff = sp.symbols("K_W_eff K_phi_eff", positive=True, real=True)
    KX, TX = sp.symbols("K_X T_X", positive=True, real=True)

    k_n = dn_halfwave_momentum(n, L)
    chi_n = sp.sqrt(2 / L) * sp.sin(k_n * s)
    print("chi_n(s) =", chi_n)
    print("k_n =", k_n)

    expect_zero(
        "k_n - (n+1/2) pi / L",
        k_n - (n + sp.Rational(1, 2)) * sp.pi / L,
    )

    overlap_from_integral = sp.simplify(sp.integrate(chi_n, (s, 0, L)))
    overlap_formula = uniform_dn_overlap(n, L)
    print("I_n from direct integral =", overlap_from_integral)
    print("I_n closed form =", overlap_formula)
    expect_zero("uniform overlap integral", overlap_from_integral - overlap_formula)

    I0 = uniform_dn_overlap(sp.Integer(0), L)
    ratio = sp.simplify(overlap_formula / I0)
    print("I_0 =", I0)
    print("I_n / I_0 =", ratio)
    expect_zero("overlap ratio hierarchy", ratio - overlap_ratio(n))

    lambda_W = sp.simplify(lambda_star * I0)
    lambda_phi_n = sp.simplify(lambda_star * overlap_formula)
    zeta_phys = sp.simplify(lambda_phi_n**2 * KW_eff / (lambda_W**2 * Kphi_eff))
    zeta_phys_expected = sp.simplify(KW_eff / Kphi_eff * overlap_ratio(n) ** 2)
    print("zeta_n^(phys) =", zeta_phys_expected)
    expect_zero("microscopic coherent-support law", zeta_phys - zeta_phys_expected)

    KW_twin = sp.simplify(KX + sp.pi**2 * TX / (4 * L**2))
    Kphi_twin = sp.simplify(KX + (n + sp.Rational(1, 2)) ** 2 * sp.pi**2 * TX / L**2)
    twin_gap = sp.simplify(sp.pi**2 * TX * n * (n + 1) / L**2)
    expect_zero(
        "same-operator twin stiffness relation",
        Kphi_twin - (KW_twin + twin_gap),
    )

    x_expr = sp.simplify(sp.pi**2 * TX / (L**2 * KW_twin))
    zeta_twin = sp.simplify(KW_twin / Kphi_twin * overlap_ratio(n) ** 2)
    zeta_twin_expected = twin_support_ratio(n, x_expr)
    print("x =", x_expr)
    print("zeta_n^(twin) =", zeta_twin_expected)
    expect_zero("exact twin-lane support ratio", zeta_twin - zeta_twin_expected)
    expect_zero("lowest twin half-wave", zeta_twin_expected.subs(n, 0) - 1)

    print("\nCarry-forward formulas:")
    print("  I_n = sqrt(2L) / ((n+1/2) pi)")
    print("  I_n / I_0 = 1 / (2n+1)")
    print("  zeta_n^(phys) = (K_W_eff / K_phi,n_eff) (I_n / I_0)^2")
    print("  zeta_n^(twin) = 1 / ((2n+1)^2 (1 + x n(n+1)))")
    print("  zeta_0^(twin) = 1")


if __name__ == "__main__":
    main()
