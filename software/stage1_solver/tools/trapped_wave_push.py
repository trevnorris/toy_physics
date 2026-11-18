"""Independent wave and rim coordinates driven by an external force on the rim."""

import sympy as sp


def show(label, value, expected):
    result = sp.factor(sp.trigsimp(value))
    print(f"{label} = {result}", flush=True)
    if sp.simplify(result - expected) != 0:
        raise ArithmeticError(label)


def main():
    # z is time in units of the relative-mode timescale. Both field inertias
    # remain free; alpha is the wave fraction of their sum.
    s, z = sp.symbols("s z", positive=True)
    alpha, f = sp.symbols("alpha f", real=True)
    xhat, yhat = sp.symbols("xhat yhat")
    initial_rest = [s**2 * xhat + alpha * (xhat - yhat) - f / s,
                    s**2 * yhat - (1 - alpha) * (xhat - yhat)]
    sol = sp.solve(initial_rest, (xhat, yhat), dict=True)[0]
    x, y = [sp.inverse_laplace_transform(sol[q], s, z) for q in (xhat, yhat)]
    print("normalized rim trajectory =", x)
    print("normalized envelope trajectory =", y)
    show("total momentum integral", (1-alpha)*x + alpha*y, (1-alpha)*f*z**2/2)
    show("envelope acceleration from slip", sp.diff(y, z, 2) - (1-alpha)*(x-y), 0)
    show("rim second derivative at onset", sp.diff(x,z,2).subs(z,0), f)
    show("envelope second derivative at onset", sp.diff(y,z,2).subs(z,0), 0)
    show("envelope fourth derivative at onset", sp.diff(y,z,4).subs(z,0), (1-alpha)*f)
    for drive in [0, 1, -1, sp.Rational(1,2), -sp.Rational(1,2)]:
        gap = (x - y).subs(f, drive)
        show(f"rim minus envelope at pi for drive {drive}", gap.subs(z,sp.pi), 2*drive)
        show(f"rim minus envelope at 2pi for drive {drive}", gap.subs(z,2*sp.pi), 0)
    # Independent trajectory samples can resolve the two field inertias in
    # the twice-integrated total momentum balance if this determinant is nonzero.
    sample = sp.Matrix([[u.subs({z:phase,f:1}) for u in (x,y)]
                        for phase in (sp.pi, 2*sp.pi)])
    show("full drive two-time position determinant", sample.det(), 4*sp.pi**2*(1-alpha))


if __name__ == "__main__":
    main()
