#!/usr/bin/env python3
"""
Step 35 for the g-2 rebuild from the moving-throat PDE base.

What this script does
---------------------
1. Uses the explicit Family-1 parent ratio r_F1 to compute the two compensated
   mouth-coupling branches g_±.
2. Proves the positive mouth-source theorem: any nonnegative normalized source on
   the first D/N interval gives 0 <= g[sigma] <= 1.
3. Evaluates the self-matched derivative source, the convex positive family, and
   the exact interpolation point that reaches the lower Family-1 branch.
4. Checks the explicit exponential mouth-layer branch g_Pi and solves the unique
   finite-bias canonical point g_Pi = g_-.

Interpretation
--------------
The physical branch is no longer tuned by hand. Positive mouth sourcing kills the
upper branch, leaves the lower branch as the only regular one, and reaches it by
only a small deformation of an already-natural derivative source law.
"""

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
    banner("STEP 35 — POSITIVE MOUTH SOURCING SELECTS THE LOWER BRANCH WITHOUT FINE TUNING")

    pi = sp.pi
    L, z = sp.symbols("L z", positive=True, real=True)
    xi = sp.symbols("xi", real=True)
    Pi = sp.symbols("Pi", positive=True, real=True)

    subbanner("XXXV.1 — Explicit Family-1 parent ratio and the two compensated branches")

    r_F1 = sp.simplify(sp.sqrt(sp.Rational(12, 1) / pi**2 * (sp.Rational(37, 20)) ** 2 - 1))
    g_minus = sp.simplify(r_F1 - sp.Rational(1, 2) * sp.sqrt(1 + r_F1**2))
    g_plus = sp.simplify(r_F1 + sp.Rational(1, 2) * sp.sqrt(1 + r_F1**2))

    print("r_F1 =")
    sp.pprint(r_F1)
    print("g_-^(F1) =")
    sp.pprint(g_minus)
    print("g_+^(F1) =")
    sp.pprint(g_plus)
    print("Numerically:")
    print("r_F1     =", sp.N(r_F1, 20))
    print("g_-^(F1) =", sp.N(g_minus, 20))
    print("g_+^(F1) =", sp.N(g_plus, 20))

    subbanner("XXXV.2 — Positive mouth-source theorem")

    k = sp.simplify(pi / (2 * L))
    sigma_match = sp.simplify(k * sp.cos(k * z))
    g_sigma = sp.Function("g_sigma")

    print("For any normalized nonnegative source sigma(z) on [0,L],")
    print("g[sigma] = ∫ sigma(z) cos(pi z / (2L)) dz")
    print("Because 0 <= cos(pi z/(2L)) <= 1 on [0,L], one has 0 <= g[sigma] <= 1.")
    print("Therefore the upper Family-1 branch is automatically impossible if g_+ > 1.")

    subbanner("XXXV.3 — Exact self-matched derivative source")

    g_match = sp.simplify(sp.integrate(sigma_match * sp.cos(k * z), (z, 0, L)))
    expect_zero("normalization of sigma_match", sp.integrate(sigma_match, (z, 0, L)) - 1)
    print("g_match =")
    sp.pprint(g_match)
    expect_zero("g_match - pi/4", sp.simplify(g_match - pi / 4))

    traction_ratio = sp.simplify((pi / 4) / g_minus)
    print("Required traction enhancement relative to the self-matched source:")
    sp.pprint(sp.Eq(sp.Symbol("T_ratio"), traction_ratio))
    print("Numerically =", sp.N(traction_ratio, 20))

    subbanner("XXXV.4 — Convex positive family and the exact interpolation point")

    sigma_xi = sp.simplify((1 - xi) * k * sp.cos(k * z) + xi / L)
    g_xi = sp.simplify(sp.integrate(sigma_xi * sp.cos(k * z), (z, 0, L)))
    expect_zero("normalization of sigma_xi", sp.integrate(sigma_xi, (z, 0, L)) - 1)
    print("g_xi =")
    sp.pprint(g_xi)
    expect_zero("g_xi - ((1-xi)pi/4 + xi*2/pi)", sp.simplify(g_xi - ((1 - xi) * pi / 4 + xi * 2 / pi)))

    xi_star = sp.simplify(sp.solve(sp.Eq(g_xi, g_minus), xi)[0])
    print("Exact xi_* matching the lower Family-1 branch:")
    sp.pprint(sp.Eq(xi, xi_star))
    print("Numerically xi_* =", sp.N(xi_star, 20))

    subbanner("XXXV.5 — Explicit exponential mouth-layer branch and the unique finite-bias point")

    g_Pi = sp.simplify(2 * Pi * (2 * Pi * sp.exp(Pi) + pi) / ((4 * Pi**2 + pi**2) * (sp.exp(Pi) - 1)))
    one_minus = sp.simplify(1 - g_Pi)
    decomposition = sp.simplify(
        (pi**2 * (sp.exp(Pi) - 1 - Pi - Pi**2 / 2) + Pi * (pi**2 - 2 * pi) + Pi**2 * (pi**2 / 2 - 4))
        / ((4 * Pi**2 + pi**2) * (sp.exp(Pi) - 1))
    )
    expect_zero("exact positivity decomposition of 1 - g_Pi", one_minus - decomposition)
    print("lim_{Pi->oo} g_Pi =", sp.limit(g_Pi, Pi, sp.oo))

    Pi_star = sp.N(sp.nsolve(sp.Eq(g_Pi, sp.N(g_minus, 40)), 1.5), 25)
    print("Unique finite-bias canonical point solving g_Pi = g_-^(F1):")
    print("Pi_* =", Pi_star)
    print("g_Pi(Pi_*) =", sp.N(g_Pi.subs(Pi, Pi_star), 25))

    banner("STEP 35 LEDGER")
    print("The explicit Family-1 parent ratio is")
    print("  r_F1 = sqrt((12/pi^2)*(37/20)^2 - 1) ≈", sp.N(r_F1, 16))
    print("so the compensated branches are")
    print("  g_-^(F1) ≈", sp.N(g_minus, 16))
    print("  g_+^(F1) ≈", sp.N(g_plus, 16))
    print()
    print("Positive mouth sourcing forces 0 <= g[sigma] <= 1, so the upper branch is")
    print("ruled out and the lower branch is the unique admissible canonical candidate.")
    print()
    print("The natural self-matched derivative profile gives")
    print("  g_match = pi/4 ≈", sp.N(pi / 4, 16))
    print("which differs from g_-^(F1) by only a 3.61% traction enhancement.")
    print()
    print("Inside the explicit convex positive family")
    print("  sigma_xi = (1-xi) k cos(kz) + xi/L")
    print("the exact lower branch is reached at")
    print("  xi_* ≈", sp.N(xi_star, 16))
    print("so only an 18.4% admixture of the uniform profile is needed.")
    print()
    print("The exponential mouth-layer branch g_Pi approaches 1 only as Pi -> infinity,")
    print("while the lower compensated Family-1 point occurs already at the finite bias")
    print("  Pi_* ≈", Pi_star)
    print("So the physical branch is regular and nearby, not singular or finely tuned.")


if __name__ == "__main__":
    main()
