#!/usr/bin/env python3
"""
Step 27 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Starts from the actual isotropic grouped-P2 one-pole conservative module
       Yhat_Q^cons(omega) = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2),
   and derives its exact low-frequency coefficients.
2. Combines that conservative module with the carried minimal outgoing grouped-P2
   branch to show that the odd Burke-Thorne coefficient is not independent:
       Gammabar_5 = 9 Kbar_0 / (32 Omega_Q^5).
3. Defines the GR target normalization
       Kbar_0^target = 64 G Omega_Q^5 / (45 c^5)
   and proves that every low-frequency coefficient on the actual isotropic
   branch scales with the same scalar defect
       N_Q = Kbar_0 / Kbar_0^target.
4. Shows that the full reduced 2.5PN/4PN theorem gap on that branch is therefore
   equivalent to the single condition
       N_Q = 1.

Interpretation
--------------
This is the sharpest conservative normalization statement reached so far. Once the
actual isotropic grouped-P2 one-pole branch is accepted, there is no remaining
multi-parameter conservative ambiguity: the even moments and the odd Burke-Thorne
coefficient are all controlled by the one scalar normalization defect N_Q.
"""

from __future__ import annotations

import sympy as sp


PI = sp.pi


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


def coeff(series_expr: sp.Expr, symbol: sp.Symbol, power: int) -> sp.Expr:
    return sp.simplify(sp.expand(series_expr).coeff(symbol, power))


def main() -> None:
    banner("STEP 27 — THE ACTUAL ISOTROPIC BRANCH HAS ONE NORMALIZATION DEFECT")

    omega, Omega_Q = sp.symbols("omega Omega_Q", positive=True, real=True)
    Kbar0 = sp.symbols("Kbar_0", positive=True, real=True)
    G, c, cs, a_th = sp.symbols("G c c_s a_th", positive=True, real=True)
    N_Q = sp.symbols("N_Q", positive=True, real=True)

    subbanner("XXVII.1 — Exact low-frequency coefficients of the actual conservative module")

    Yhat_cons = sp.simplify(sp.Rational(3, 4) + sp.Rational(1, 4) / (1 - omega**2 / Omega_Q**2))
    Yhat_cons_series = sp.expand(sp.series(Yhat_cons, omega, 0, 6).removeO())
    print("Yhat_Q^cons(omega) =")
    sp.pprint(Yhat_cons)
    print("Series through O(omega^4) =")
    sp.pprint(Yhat_cons_series)

    u2 = coeff(Yhat_cons_series, omega, 2)
    u4 = coeff(Yhat_cons_series, omega, 4)
    print("u_2 =")
    sp.pprint(u2)
    print("u_4 =")
    sp.pprint(u4)
    expect_zero("u_2 - 1/(4 Omega_Q^2)", sp.simplify(u2 - 1 / (4 * Omega_Q**2)))
    expect_zero("u_4 - 1/(4 Omega_Q^4)", sp.simplify(u4 - 1 / (4 * Omega_Q**4)))

    Kbar2 = sp.simplify(Kbar0 * u2)
    Kbar4 = sp.simplify(Kbar0 * u4)
    print("Kbar_2 =")
    sp.pprint(Kbar2)
    print("Kbar_4 =")
    sp.pprint(Kbar4)
    expect_zero("Kbar_2 - Kbar_0/(4 Omega_Q^2)", sp.simplify(Kbar2 - Kbar0 / (4 * Omega_Q**2)))
    expect_zero("Kbar_4 - Kbar_0/(4 Omega_Q^4)", sp.simplify(Kbar4 - Kbar0 / (4 * Omega_Q**4)))

    subbanner("XXVII.2 — The odd Burke–Thorne coefficient is fixed by the same branch")

    Gammabar5 = sp.simplify(9 * Kbar2**sp.Rational(5, 2) / Kbar0**sp.Rational(3, 2))
    print("Gammabar_5 from the minimal isotropic outgoing branch =")
    sp.pprint(Gammabar5)
    expect_zero(
        "Gammabar_5 - 9 Kbar_0/(32 Omega_Q^5)",
        sp.simplify(Gammabar5 - 9 * Kbar0 / (32 * Omega_Q**5)),
    )

    subbanner("XXVII.3 — Exact GR target normalization")

    Kbar0_target = sp.simplify(64 * G * Omega_Q**5 / (45 * c**5))
    Kbar2_target = sp.simplify(Kbar0_target / (4 * Omega_Q**2))
    Kbar4_target = sp.simplify(Kbar0_target / (4 * Omega_Q**4))
    Gammabar5_target = sp.simplify(2 * G / (5 * c**5))

    print("Kbar_0^target =")
    sp.pprint(Kbar0_target)
    print("Kbar_2^target =")
    sp.pprint(Kbar2_target)
    print("Kbar_4^target =")
    sp.pprint(Kbar4_target)
    print("Gammabar_5^target =")
    sp.pprint(Gammabar5_target)

    Kbar0_target_geo = sp.simplify(Kbar0_target.subs(Omega_Q, 3 * cs / (2 * a_th)))
    print("With Omega_Q = 3 c_s/(2 a_th):")
    print("Kbar_0^target =")
    sp.pprint(Kbar0_target_geo)
    expect_zero(
        "Kbar_0^target - 54 G c_s^5/(5 a_th^5 c^5)",
        sp.simplify(Kbar0_target_geo - 54 * G * cs**5 / (5 * a_th**5 * c**5)),
    )

    subbanner("XXVII.4 — The full actual-branch mismatch is one scalar N_Q")

    Kbar0_actual = sp.simplify(N_Q * Kbar0_target)
    Kbar2_actual = sp.simplify(Kbar0_actual / (4 * Omega_Q**2))
    Kbar4_actual = sp.simplify(Kbar0_actual / (4 * Omega_Q**4))
    Gammabar5_actual = sp.simplify(9 * Kbar0_actual / (32 * Omega_Q**5))

    print("Actual isotropic branch with normalization defect N_Q:")
    print("Kbar_0 =")
    sp.pprint(Kbar0_actual)
    print("Kbar_2 =")
    sp.pprint(Kbar2_actual)
    print("Kbar_4 =")
    sp.pprint(Kbar4_actual)
    print("Gammabar_5 =")
    sp.pprint(Gammabar5_actual)

    expect_zero("Kbar_0/Kbar_0^target - N_Q", sp.simplify(Kbar0_actual / Kbar0_target - N_Q))
    expect_zero("Kbar_2/Kbar_2^target - N_Q", sp.simplify(Kbar2_actual / Kbar2_target - N_Q))
    expect_zero("Kbar_4/Kbar_4^target - N_Q", sp.simplify(Kbar4_actual / Kbar4_target - N_Q))
    expect_zero(
        "Gammabar_5/Gammabar_5^target - N_Q",
        sp.simplify(Gammabar5_actual / Gammabar5_target - N_Q),
    )

    R0 = sp.simplify(Kbar0_actual / Kbar0_target - 1)
    R2 = sp.simplify(Kbar2_actual / Kbar2_target - 1)
    R4 = sp.simplify(Kbar4_actual / Kbar4_target - 1)
    R5 = sp.simplify(Gammabar5_actual / Gammabar5_target - 1)
    expect_zero("R0 - (N_Q-1)", sp.simplify(R0 - (N_Q - 1)))
    expect_zero("R2 - (N_Q-1)", sp.simplify(R2 - (N_Q - 1)))
    expect_zero("R4 - (N_Q-1)", sp.simplify(R4 - (N_Q - 1)))
    expect_zero("R5 - (N_Q-1)", sp.simplify(R5 - (N_Q - 1)))

    banner("STEP 27 LEDGER")
    print("Actual isotropic grouped-P2 one-pole conservative module:")
    print("  Yhat_Q^cons(omega) = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)")
    print()
    print("Exact low-frequency coefficients on that branch:")
    print("  Kbar_2 = Kbar_0/(4 Omega_Q^2)")
    print("  Kbar_4 = Kbar_0/(4 Omega_Q^4)")
    print("  Gammabar_5 = 9 Kbar_0/(32 Omega_Q^5)")
    print()
    print("GR target normalization:")
    print("  Kbar_0^target = 64 G Omega_Q^5/(45 c^5)")
    print("  = 54 G c_s^5/(5 a_th^5 c^5) after Omega_Q = 3 c_s/(2 a_th)")
    print()
    print("Single actual-branch defect:")
    print("  N_Q := Kbar_0 / Kbar_0^target")
    print()
    print("And then automatically:")
    print("  Kbar_2/Kbar_2^target = N_Q")
    print("  Kbar_4/Kbar_4^target = N_Q")
    print("  Gammabar_5/Gammabar_5^target = N_Q")
    print()
    print("So the entire reduced 2.5PN/4PN theorem gap on the actual isotropic branch")
    print("is equivalent to the one scalar condition N_Q = 1.")


if __name__ == "__main__":
    main()
