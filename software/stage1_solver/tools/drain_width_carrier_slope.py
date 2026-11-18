"""Imposed broad thickness changes at fixed radial throughput, local far-field limit."""

import sympy as sp


def main():
    n, width, current = sp.symbols("n W I", positive=True)
    pressure = (n**5 - 1) / 5
    speed = current / (width*n)
    bernoulli = (n**4-1)/4 + speed**2/2
    slope = sp.factor(-sp.diff(pressure, n)*sp.diff(bernoulli, width)/sp.diff(bernoulli, n))
    predicted = n**5*speed**2/(width*(n**4-speed**2))
    print("carrier dP/dW =", slope)
    assert sp.simplify(slope - predicted) == 0

    # Solve independently at each width with the same distant density,
    # chemical potential and intake. Choose the branch connected to n=1.
    for radius in (8, 16):
        flux = 1/(16*sp.pi*radius**2)  # Q=1/4 divided by the 3D sphere area
        def state(w):
            equation = bernoulli.subs({width: w, current: flux})
            value = sp.nsolve(equation, n, 1, prec=55)
            if not (value > 0 and value**4 > speed.subs({current: flux, width: w, n: value})**2
                    and abs(equation.subs(n, value).evalf(50)) < 1e-48):
                raise ValueError("Unresolved subsonic Bernoulli state")
            return value, pressure.subs(n, value)
        n0, p0 = state(1)
        exact = slope.subs({n: n0, width: 1, current: flux}).evalf(35)
        print(f"r={radius} baseline_excess_pressure={float(p0):+.12g} width_slope={float(exact):+.12g}")
        for delta in (sp.Rational(1, 50), sp.Rational(1, 100)):
            _, wider = state(1+delta)
            _, narrower = state(1-delta)
            central = (wider - narrower)/(2*delta)
            print(f"  +/-delta={float(delta):g} wider_dP={float(wider-p0):+.12g} "
                  f"narrower_dP={float(narrower-p0):+.12g} finite_slope={float(central):+.12g}")
            if abs(central/exact-1) > 3*delta**2:
                raise ValueError("Width pairs disagree with the linear slope")


if __name__ == "__main__":
    main()
