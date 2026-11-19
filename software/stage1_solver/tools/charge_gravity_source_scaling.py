"""Test proportional electric and drain source strengths across two physical mouths."""

import sympy as sp


def main():
    m1, m2, e, x1, x2, a, b, k = sp.symbols("M1 M2 e x1 x2 a b k", positive=True)
    # Columns are two matched-ambient preparations. The second row contains
    # electric magnitudes so reflecting a mouth does not change the rank test.
    shared = sp.Matrix([[a*x1, a*x2], [b*x1, b*x2]])
    required = sp.Matrix([[k*m1, k*m2], [e, e]])
    shared_minor = sp.factor(shared.det())
    required_minor = sp.factor(required.det())
    print("proportional source rows, determinant =", shared_minor, flush=True)
    print("target source rows, determinant =", required_minor, flush=True)
    print("equal-inertia control =", required_minor.subs(m2, m1), flush=True)
    # Relax only the assumption that electric strength is proportional to x.
    # Find the response coefficient ratio needed to keep one electric unit.
    b2 = sp.symbols("b2", positive=True)
    ratio, = sp.solve(b*x1-b2*x2, b2)
    electric_gain_ratio = sp.factor((ratio/b).subs(x2, x1*m2/m1))
    print("electric gain ratio required if drain strength tracks x =", electric_gain_ratio, flush=True)
    assert shared_minor == 0
    assert sp.simplify(required_minor - k*e*(m1-m2)) == 0
    assert required_minor.subs(m2, m1) == 0
    assert sp.simplify(electric_gain_ratio-m1/m2) == 0


if __name__ == "__main__":
    main()
