"""Axial momentum export by a rounded 4D receiver, at fixed mass throughput."""

import sympy as sp


def check(name, expression):
    residual = sp.factor(sp.simplify(expression))
    print(f"{name} = {residual}", flush=True)
    if residual != 0:
        raise ArithmeticError(name)


def hemisphere():
    # In four space dimensions dOmega = 4*pi*sqrt(1-mu**2) dmu.
    # Use exact beta moments to integrate the velocity and stress on the cut.
    mu = sp.symbols("mu", real=True)
    r, a, mass, rho = sp.symbols("R a D_a rho", positive=True)
    e = sp.symbols("epsilon", real=True)

    def integrate(value):
        terms = sp.Poly(value, mu).terms()
        total = 0
        for (power,), coeff in terms:
            arg = sp.Rational(power+1, 2)
            moment = sp.gamma(arg)*sp.gamma(sp.Rational(3, 2))/(2*sp.gamma(arg+sp.Rational(3, 2)))
            total += coeff*moment
        return sp.factor(4*sp.pi*r**3*total)

    q = mass/(sp.pi**2*rho)
    potential = -q/(2*r**2)-q*e*a**2*(4*mu**2-1)/(4*r**4)
    radial = potential.diff(r)
    tangent = potential.diff(mu)/r  # coefficient along e_w - mu*e_r
    axial = mu*radial + (1-mu**2)*tangent
    speed_squared = radial**2+(1-mu**2)*tangent**2
    carriage = integrate(rho*axial*radial)
    pressure = integrate(-rho*mu*speed_squared/2)
    flux = sp.factor(carriage+pressure)
    print("outward mass rate =", integrate(rho*radial))
    print("upper hemisphere carriage =", carriage)
    print("upper hemisphere excess pressure =", pressure)
    print("upper hemisphere total =", sp.collect(sp.expand(flux), e))
    check("mass integral residual", integrate(rho*radial)-mass)
    # Reflecting equatorial plane, with normal out of the fluid along -w.
    plane_pressure = (-rho*speed_squared/2).subs(mu, 0)
    check("hemisphere plus plane-annulus conservation",
          flux.diff(r)-4*sp.pi*r**2*plane_pressure)
    check("far axial flux", sp.limit(flux,r,sp.oo))
    uniform = flux.subs({r:a,e:0})
    print("uniform inner-cut flux =", uniform)
    print("same-rate shaped / uniform inner flux =",
          sp.factor(flux.subs({r:a,e:sp.Rational(1,4)})/uniform))


def carrier():
    # Weak stationary radial discharge. Use v=O(delta), eta=n-n0=O(delta^2)
    # and uniform chemical potential at a quiescent infinity.
    r, mass, m, n0, k, hp, xi = sp.symbols("R D_a m n0 k hprime xi", positive=True)
    eta = sp.Function("eta")(r)
    v = mass/(m*n0*sp.pi**2*r**3)
    bernoulli = hp*eta+m*v**2/2-k*(eta.diff(r,2)+3*eta.diff(r)/r)/n0
    rr = m*n0*v**2+n0*hp*eta-k*eta.diff(r,2)
    tangential = n0*hp*eta-k*eta.diff(r)/r
    curvature = sp.solve(bernoulli,eta.diff(r,2))[0]
    flux = sp.expand((4*sp.pi*r**3/3)*rr.subs(eta.diff(r,2),curvature))
    print("weak radial carrier hemisphere flux =",flux)
    check("carrier curved-cut and flat-annulus balance",
          (flux.diff(r)-4*sp.pi*r**2*tangential).subs(eta.diff(r,2),curvature))
    # For two runs at the same D_a, the inhomogeneous Bernoulli term cancels.
    # b''+b'/x-(1+1/x^2)b=0 is the order-one modified Bessel equation.
    x, amp = sp.symbols("x A", positive=True)
    b = sp.Function("b")(x)
    mode = b/x
    ode = b + b/x**2-b.diff(x)/x
    check("decaying density mode equation",
          (mode.diff(x,2)+3*mode.diff(x)/x-mode).subs(b.diff(x,2),ode))
    mode = amp*sp.besselk(1,r/xi)/r
    shift = 4*sp.pi*k*r**2*mode.diff(r)
    shift = sp.simplify(sp.expand_func(shift+4*sp.pi*k*amp*r*sp.besselk(2,r/xi)/xi))
    check("density-mode export minus closed form", shift)
    print("same-rate flux shift =", -4*sp.pi*k*amp*r*sp.besselk(2,r/xi)/xi)


if __name__ == "__main__":
    hemisphere()
    carrier()
