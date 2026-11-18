"""What a rim drive can determine about a bound envelope's follow rate."""

import sympy as sp


def check(label, value, target):
    result = sp.factor(value)
    print(f"{label} = {result}", flush=True)
    if sp.simplify(result - target) != 0:
        raise ArithmeticError(label)


def main():
    mr, mw, k, s, nu, eps = sp.symbols("Mr Mw k s nu eps", positive=True)
    x, y, force = sp.symbols("X Y F")
    sol = sp.solve([mr*s**2*x + k*(x-y) - force, mw*s**2*y - k*(x-y)],
                   (x, y), dict=True)[0]
    follow = sp.factor(sol[y]/sol[x])
    app = sp.factor(force/(s**2*sol[x]))
    check("envelope/rim for zero initial disturbance", follow, k/(k+mw*s**2))
    check("rim force per acceleration", app, mr+mw*follow)
    gap2 = sp.Symbol("g", positive=True)  # squared clamped-rim follow frequency
    harmonic_error = (follow.subs({s: sp.I*nu, k: mw*gap2})-1).factor()
    check("signed harmonic excess below pole", harmonic_error, nu**2/(gap2-nu**2))
    check("squared follow rate at fractional tolerance",
          sp.solve(harmonic_error-eps, gap2)[0], nu**2*(1+1/eps))
    m0 = app.subs(s, 0)
    ratio = sp.factor(((m0-app)/(m0-app.subs(s,2*s))).subs(k,mw*gap2))
    r = sp.Symbol("R", positive=True)
    inferred = sp.solve(r-ratio, gap2)[0]
    check("ratio of two apparent-inertia deficits", ratio, (gap2+4*s**2)/(4*(gap2+s**2)))
    check("squared rate from mouth readings", inferred, 4*s**2*(1-r)/(4*r-1))
    check("slope of squared rate with measured ratio", sp.diff(inferred,r), -12*s**2/(4*r-1)**2)
    # Step response from the forced two-coordinate dynamics, normalized to
    # the displacement of their common center of inertia.
    t = sp.Symbol("t", positive=True)
    omega = sp.sqrt(k*(1/mr+1/mw))
    diff = sp.inverse_laplace_transform((sol[x]-sol[y])/s, s, t)
    cm = force*t**2/(2*(mr+mw))
    check("relative displacement divided by center displacement",
          diff/cm, 2*(1-sp.cos(omega*t))/(k/mw*t**2))
    # Change relative stiffness without changing centered radial equilibrium.
    a, d, lam = sp.symbols("a delta lam", real=True)
    energy = sp.Function("E0")(a) + lam*sp.Function("b")(a)*d**2/2
    check("rest radial derivative", sp.diff(energy,a).subs(d,0), sp.diff(sp.Function("E0")(a),a))
    check("relative curvature", sp.diff(energy,d,2).subs(d,0), lam*sp.Function("b")(a))


if __name__ == "__main__":
    main()
