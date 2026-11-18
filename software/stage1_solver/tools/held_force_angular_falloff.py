"""Exact integral of the steady momentum tensor on the midplane between mouths."""
import sympy as sp


def area(n):
    return 2*sp.pi**(sp.Rational(n, 2))/sp.gamma(sp.Rational(n, 2))


def check(label, value, expected):
    residual = sp.factor(value-expected)
    print(f"{label}: {sp.collect(sp.expand(value), distance)} | residual {residual}", flush=True)
    if residual != 0:
        raise ArithmeticError(label)


distance, rho = sp.symbols("d rho", positive=True)
qt, qs, bt, bs = sp.symbols("q_t q_s b_t b_s", real=True)


def midplane(n):
    # On x=d/2 use t=transverse radius/(d/2). Factor the same radial
    # denominator out of both velocity components before the stress integral.
    t = sp.Symbol("t", positive=True)
    half = distance/2
    common = 1+t**2

    def velocity_numerators(q, b, side):
        mono = q*half*common/(rho*area(n))
        return side*mono+b*(common-n), t*(mono-b*n*side)

    def integrate(numerator):
        poly = sp.Poly(numerator, t)
        terms = []
        for (power,), coeff in poly.terms():
            arg = sp.Rational(power+n-1, 2)
            # integral of t^(power+n-2)/(1+t^2)^(n+2) from 0 to infinity
            terms.append(coeff*sp.gamma(arg)*sp.gamma(n+2-arg)/(2*sp.gamma(n+2)))
        return sp.factor(rho*area(n-1)*sum(terms)/half**(n+1))

    tx, ty = velocity_numerators(qt, bt, 1)
    sx, sy = velocity_numerators(qs, bs, -1)
    # Pi_xx = rho*(vx^2-vy^2)/2, with the common far pressure removed.
    self_load = integrate((tx**2-ty**2)/2)
    check(f"{n}D one-mouth plane integral at arbitrary b_t", self_load, 0)
    return integrate(tx*sx-ty*sy)


def main():
    print("SymPy version =", sp.__version__, flush=True)
    for n in (3,4):
        load = midplane(n)
        mono = -qt*qs/(rho*area(n)*distance**(n-1))
        linear = (n-1)*(qs*bt-qt*bs)/distance**n
        quadratic = n*(n-1)*rho*area(n)*bt*bs/distance**(n+1)
        check(f"{n}D finite paired-minus-alone holder", load, mono+linear+quadratic)
        reciprocal = load.xreplace({qt:qs, qs:qt, bt:-bs, bs:-bt})
        check(f"{n}D same attraction in reversed center coordinates", reciprocal-load, 0)
        limit = sp.limit(distance**2*load, distance, sp.oo)
        check(f"{n}D inverse-square holder coefficient for fixed moments",
              limit, -qt*qs/(rho*area(n)) if n==3 else 0)
        # Probe an angular redistribution that itself responds as 1/d^2.
        kt, ks = sp.symbols("k_t k_s", real=True)
        induced = (load-mono).subs({bt:kt/distance**2, bs:ks/distance**2})
        check(f"{n}D rescaled leading induced correction",
              sp.limit(distance**(n+2)*induced, distance, sp.oo), (n-1)*(qs*kt-qt*ks))


if __name__=="__main__":
    main()
