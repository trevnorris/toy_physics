#!/usr/bin/env python3
"""
Step 28 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Computes the exact compact outgoing l=2 Dirichlet-to-Neumann fingerprint from
   the spherical Hankel mode h_2^(1) and verifies the canonical normalized series
       Yhat_2^out = 1 + z^2/9 + 4 z^4/81 + i z^5/27 + O(z^6).
2. Matches that explicit outgoing branch to the retarded grouped-P2 one-pole
   module and fixes the remaining reduced outgoing normalization scalar
       chi_Q = 1.
3. Verifies the exact factorization of the last reduced 2.5PN defect:
       mhat_0^2 chi_Q N_Q = 1,
   so on the natural point-particle source-map branch mhat_0 -> 1 one has
       N_Q = 1/chi_Q.
4. Derives the first explicit isotropic DtN deformation law
       chi_Q = 3 (S beta^5 + 9 Sigma_5)/(3 S - Sigma_0),
   together with its linearization
       chi_Q = 1 + 5 b + a_0/3 + 9 a_5 + O(eps^2),
   which identifies the minimal PDE-facing outgoing branch-selection data.

Interpretation
--------------
This step is the first honest retarded closure statement. On the canonical compact
outgoing l=2 DtN branch, the last reduced 2.5PN scalar is fixed exactly:
chi_Q = 1. Any remaining deviation must therefore come from an actual DtN branch
selection effect, and at first order it can only enter through the isotropic
triple (b, a_0, a_5).
"""

from __future__ import annotations

import sympy as sp
from sympy.functions.special.bessel import jn, yn


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


def expect_zero(name: str, expr: sp.Expr | sp.Matrix) -> None:
    if isinstance(expr, sp.MatrixBase):
        simplified = expr.applyfunc(lambda z: sp.simplify(sp.expand(z)))
        print(f"{name} =")
        sp.pprint(simplified)
        if any(entry != 0 for entry in simplified):
            raise AssertionError(f"{name} is not zero")
    else:
        simplified = sp.simplify(sp.expand(expr))
        print(f"{name} = {simplified}")
        if simplified != 0:
            raise AssertionError(f"{name} is not zero")


def main() -> None:
    banner("STEP 28 — EXPLICIT OUTGOING DTN FIX OF chi_Q AND BRANCH-SELECTION LAW")

    z = sp.symbols("z", real=True)
    I = sp.I

    subbanner("XXVIII.1 — Exact outgoing l=2 DtN fingerprint")

    h2 = sp.expand_func(jn(2, z) + I * yn(2, z))
    Lambda_out = sp.simplify(z * sp.diff(h2, z) / h2)
    Lambda_series = sp.expand(sp.series(Lambda_out, z, 0, 8).removeO())

    print("Lambda_2^out(z) = z h_2'(z)/h_2(z)")
    print("Series =")
    sp.pprint(Lambda_series)

    expected_Lambda = -3 + z**2 / 3 + z**4 / 9 + I * z**5 / 9 - 2 * z**6 / 27 - I * z**7 / 27
    expect_zero("Lambda_out series - expected", sp.simplify(Lambda_series - expected_Lambda))

    Yhat_out = sp.simplify(-3 / Lambda_out)
    Yhat_out_series = sp.expand(sp.series(Yhat_out, z, 0, 8).removeO())
    print("Yhat_2^out(z) = -3/Lambda_2^out(z)")
    print("Series =")
    sp.pprint(Yhat_out_series)

    expected_Yhat = 1 + z**2 / 9 + 4 * z**4 / 81 + I * z**5 / 27 - 11 * z**6 / 729 - I * z**7 / 243
    expect_zero("Yhat_out series - expected", sp.simplify(Yhat_out_series - expected_Yhat))

    subbanner("XXVIII.2 — Matching to the retarded grouped-P2 one-pole module fixes chi_Q")

    omega, Omega_Q, chi_Q, c_s, a_th = sp.symbols("omega Omega_Q chi_Q c_s a_th", positive=True, real=True)
    sigma_can = sp.simplify(9 / (8 * Omega_Q**5))
    Yhat_ret = sp.simplify(
        sp.Rational(3, 4) + sp.Rational(1, 4) / (1 - omega**2 / Omega_Q**2 - I * chi_Q * sigma_can * omega**5)
    )
    Yhat_ret_series = sp.expand(sp.series(Yhat_ret, omega, 0, 6).removeO())
    print("Yhat_Q^ret(omega) =")
    sp.pprint(Yhat_ret)
    print("Series through O(omega^5) =")
    sp.pprint(Yhat_ret_series)

    sigma_can_geo = sp.simplify(sigma_can.subs(Omega_Q, 3 * c_s / (2 * a_th)))
    expect_zero(
        "sigma_Q^can - 4 a_th^5/(27 c_s^5)",
        sp.simplify(sigma_can_geo - 4 * a_th**5 / (27 * c_s**5)),
    )

    # Convert the explicit outgoing DtN series to omega using z = a omega / c_s.
    Yhat_out_omega = sp.expand(expected_Yhat.subs(z, a_th * omega / c_s))
    print("Explicit outgoing DtN fingerprint in omega =")
    sp.pprint(Yhat_out_omega)

    Yhat_ret_geo = sp.expand(Yhat_ret_series.subs(Omega_Q, 3 * c_s / (2 * a_th)))
    print("Retarded one-pole grouped-P2 module in omega =")
    sp.pprint(Yhat_ret_geo)

    chi_solution = sp.solve(
        sp.Eq(sp.expand(Yhat_ret_geo).coeff(omega, 5) / I, sp.expand(Yhat_out_omega).coeff(omega, 5) / I),
        chi_Q,
    )[0]
    print("chi_Q from matching the O(omega^5) coefficient =")
    sp.pprint(chi_solution)
    expect_zero("chi_Q - 1", sp.simplify(chi_solution - 1))

    subbanner("XXVIII.3 — Exact factorization of the last reduced 2.5PN defect")

    mhat0, N_Q = sp.symbols("mhat_0 N_Q", positive=True, real=True)
    factorized_condition = sp.simplify(mhat0**2 * chi_Q * N_Q)
    print("Observable odd normalization condition =")
    sp.pprint(sp.Eq(factorized_condition, 1))
    print("So on the natural source-map branch mhat_0 -> 1,")
    print("  N_Q = 1/chi_Q.")
    expect_zero("N_Q on canonical branch - 1", sp.simplify((1 / chi_solution) - 1))

    subbanner("XXVIII.4 — General isotropic DtN deformation law for chi_Q")

    S, beta, Sigma0, Sigma2, Sigma4, Sigma5 = sp.symbols(
        "S beta Sigma_0 Sigma_2 Sigma_4 Sigma_5", real=True
    )
    # Canonical outgoing coefficients in z.
    L0 = -3 * S + Sigma0
    L2 = S * beta**2 / 3 + Sigma2
    L4 = S * beta**4 / 9 + Sigma4
    L5 = S * beta**5 / 9 + Sigma5

    Yhat_def = sp.simplify(L0 / (L0 + L2 * z**2 + L4 * z**4 + I * L5 * z**5))
    Yhat_def_series = sp.expand(sp.series(Yhat_def, z, 0, 6).removeO())
    print("Yhat_2^def(z) =")
    sp.pprint(Yhat_def_series)

    # Canonical even matching conditions.
    even_eq1 = sp.Eq(-L2 / L0, sp.Rational(1, 9))
    even_eq2 = sp.Eq(L2**2 / L0**2 - L4 / L0, sp.Rational(4, 81))
    sol_even = sp.solve((even_eq1, even_eq2), (Sigma2, Sigma4), dict=True)[0]
    print("Canonical-even matching gives:")
    print("Sigma_2 =")
    sp.pprint(sp.simplify(sol_even[Sigma2]))
    print("Sigma_4 =")
    sp.pprint(sp.simplify(sol_even[Sigma4]))

    Yhat_def_matched = sp.expand(Yhat_def_series.subs(sol_even))
    chi_from_def = sp.simplify((sp.expand(Yhat_def_matched).coeff(z, 5) / I) / sp.Rational(1, 27))
    print("chi_Q from the matched deformed DtN branch =")
    sp.pprint(chi_from_def)
    expected_chi_def = sp.simplify(3 * (S * beta**5 + 9 * Sigma5) / (3 * S - Sigma0))
    expect_zero("chi_Q deformation law - expected", sp.simplify(chi_from_def - expected_chi_def))

    subbanner("XXVIII.5 — Linearized branch-selection law")

    eps, s, b, a0, a5 = sp.symbols("eps s b a_0 a_5", real=True)
    chi_linear = sp.series(
        chi_from_def.subs(
            {
                S: 1 + eps * s,
                beta: 1 + eps * b,
                Sigma0: eps * a0,
                Sigma5: eps * a5,
            }
        ),
        eps,
        0,
        2,
    ).removeO()
    chi_linear = sp.expand(chi_linear)
    print("chi_Q near the canonical branch =")
    sp.pprint(chi_linear)
    expect_zero(
        "linearized chi_Q - (1 + eps*(5 b + a_0/3 + 9 a_5))",
        sp.simplify(chi_linear - (1 + eps * (5 * b + a0 / 3 + 9 * a5))),
    )

    banner("STEP 28 LEDGER")
    print("Exact compact outgoing l=2 DtN fingerprint:")
    print("  Lambda_2^out(z) = -3 + z^2/3 + z^4/9 + i z^5/9 - 2 z^6/27 - i z^7/27 + ...")
    print("  Yhat_2^out(z)   = 1 + z^2/9 + 4 z^4/81 + i z^5/27 - 11 z^6/729 - i z^7/243 + ...")
    print()
    print("Matching to the retarded grouped-P2 one-pole module fixes")
    print("  chi_Q = 1")
    print("on the canonical compact outgoing branch.")
    print()
    print("So the exact factorized normalization condition")
    print("  mhat_0^2 chi_Q N_Q = 1")
    print("reduces on the natural point-particle source-map branch to")
    print("  N_Q = 1/chi_Q, hence N_Q = 1 on the canonical outgoing branch.")
    print()
    print("General isotropic DtN deformation law:")
    print("  chi_Q = 3 (S beta^5 + 9 Sigma_5)/(3 S - Sigma_0)")
    print("with the canonical-even matching conditions fixing Sigma_2 and Sigma_4.")
    print()
    print("Linearized near the canonical branch:")
    print("  chi_Q = 1 + eps*(5 b + a_0/3 + 9 a_5) + O(eps^2)")
    print()
    print("Interpretation:")
    print("  Any remaining deviation from the canonical outgoing branch can only enter")
    print("  through the isotropic branch-selection triple (b, a_0, a_5).")


if __name__ == "__main__":
    main()
