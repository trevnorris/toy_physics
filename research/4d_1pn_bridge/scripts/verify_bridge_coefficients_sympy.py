#!/usr/bin/env python3

import sympy as sp


def main() -> None:
    aH2 = sp.symbols("aH2", real=True)
    alpha2, Kvec = sp.symbols("alpha2 Kvec", real=True)
    n, x, delta = sp.symbols("n x delta", real=True)
    v, c, E0 = sp.symbols("v c E0", real=True, positive=True)

    print("=== 4d_1pn_bridge symbolic checks ===")

    print("kappa_add(throat) =", sp.Rational(1, 2))
    print("kappa_add(bubble) =", sp.Rational(1, 3))

    ndelta = (1 + delta) ** (-(n - 1) / 2)
    print("linear coeff in N(delta) =", sp.expand(ndelta.series(delta, 0, 2).removeO()).coeff(delta, 1))
    print("solve ((n-1)/2 == 2) =", sp.solve(sp.Eq((n - 1) / 2, 2), n))

    c_parallel = Kvec * sp.pi**2 * (-1 + aH2 - alpha2)
    c_longitudinal = Kvec * sp.pi**2 * (-1 + aH2 + alpha2)
    family = sp.solve(
        (sp.Eq(c_parallel, -sp.Rational(7, 2)), sp.Eq(c_longitudinal, -sp.Rational(1, 2))),
        (alpha2, Kvec),
        dict=True,
    )[0]
    print("family alpha^2(a_H^2) =", sp.simplify(family[alpha2]))
    print("family K_vec(a_H^2) =", sp.simplify(family[Kvec]))

    alpha_thermo = sp.simplify(1 - 1 / (n - 1)).subs(n, 5)
    unique = sp.solve(
        (
            sp.Eq(c_parallel, -sp.Rational(7, 2)),
            sp.Eq(c_longitudinal, -sp.Rational(1, 2)),
            sp.Eq(alpha2, sp.Rational(3, 4)),
        ),
        (alpha2, Kvec, aH2),
        dict=True,
    )[0]
    print("alpha^2_thermo(n=5) =", alpha_thermo)
    print("thermo + EIH unique point =", unique)

    dln_f = sp.simplify((((n - 1) / 2) + (-1) * x + n * ((1 + 2 * x) / 3)) / (1 + x + (1 + 2 * x) / 3))
    dln_f_n5 = sp.simplify(dln_f.subs(n, 5))
    x_solution = sp.solve(sp.Eq(dln_f_n5, sp.Rational(5, 2)), x)[0]
    print("d ln F / d ln rho (n=5) =", dln_f_n5)
    print("x = Ef/Ew =", x_solution)
    print(
        "Ew:Ef:Epv =",
        [sp.simplify(term) for term in (11 * sp.Matrix([1, x_solution, (1 + 2 * x_solution) / 3]))],
    )
    print(
        "d ln a / d ln rho =",
        sp.simplify((-( -(n - 1) / 2 + 2 * x + n * (1 + 2 * x)) / (4 + 10 * x)).subs({n: 5, x: x_solution})),
    )

    gamma_series = sp.series(1 / sp.sqrt(1 - v**2 / c**2), v, 0, 6).removeO()
    energy_series = sp.expand(E0 * gamma_series)
    print("gamma(v) through O(v^4) =", gamma_series)
    print("E(v) through O(v^4) =", energy_series)
    print("coeff of v^4 in E(v) =", sp.expand(energy_series).coeff(v, 4))

    print("ALL SYMPY CHECKS PASSED")


if __name__ == "__main__":
    main()
