"""Constrained rim reaction and a field-only ambiguity on a prescribed path."""

import sympy as sp


def reaction():
    # A radial variation changes volume by R^2 dR dOmega and its x moment
    # by mu R^3 dR dOmega. Differentiate the centroid before constraining it.
    moment, volume, r, mu, e = sp.symbols("I V R mu e", real=True)
    centroid = (moment + e*mu*r**3)/(volume + e*r**2)
    gradient = sp.factor(centroid.diff(e).subs(e, 0))
    a = sp.Symbol("a", positive=True)
    sphere = gradient.subs({moment: 0, volume: 4*sp.pi*a**3/3, r: a})
    unit_shift = sp.integrate(2*sp.pi*mu*sphere, (mu,-1,1))
    # A complete residual includes every conservative and open-system
    # contribution except the work of the commanded actuator.
    gx, gy, lam, h = sp.symbols("G_x G_y lambda h", nonzero=True)
    g = sp.Matrix([gx, gy])
    lift = sp.Matrix([(1-h*gy)/gx, h])
    virtual_work = sp.expand(lift.dot(lam*g))
    print("centroid gradient per solid angle =", gradient)
    print("rigid translation normalization, work for a unit lift =", unit_shift, virtual_work)
    if unit_shift != 1 or virtual_work != lam:
        raise ArithmeticError("Geometric translation force is misnormalized")


def alternative_reactions():
    # y is a brane field coordinate normalized to kinetic energy ydot^2/2.
    # The interface options below are illustrative, not fitted to the drain.
    t, s, omega = sp.symbols("t s omega", positive=True)
    A, c = sp.symbols("A c", real=True)
    x, y, ell = [sp.Function(name)(t) for name in ("X", "y", "ell")]
    kinetic = y.diff(t)**2/2
    rigid = kinetic + ell*(y-c*x)
    compliant = kinetic - omega**2*(y-c*x)**2/2
    def residual(lag, coord):
        return sp.diff(lag.diff(coord.diff(t)), t) - lag.diff(coord)
    multiplier, = sp.solve(residual(rigid, y), ell)
    locked_force = residual(rigid, x).subs(ell, multiplier).subs(y, c*x).doit()
    # Laplace transfer from the center to the field, starting from rest.
    xs, ys = sp.symbols("X_s y_s")
    zero_state = {y.diff(t, 2): s**2*ys, y: ys, x: xs}
    solved_y, = sp.solve(residual(compliant, y).xreplace(zero_state), ys)
    force_s = residual(compliant, x).xreplace(zero_state).subs(ys, solved_y)
    per_accel = sp.factor(force_s/(s**2*xs))
    spring_step = sp.simplify(sp.inverse_laplace_transform(A*per_accel/s, s, t))
    locked_step = locked_force.subs(x.diff(t, 2), A)
    energy = kinetic + omega**2*(y-c*x)**2/2
    power = sp.expand(energy.diff(t)-residual(compliant, x)*x.diff(t)
                      -residual(compliant, y)*y.diff(t))
    print("locked and compliant extra force per specified acceleration =",
          sp.factor(locked_step/A), sp.factor(spring_step/A))
    print("compliant response per acceleration, slow and fast Laplace limits =",
          per_accel, sp.limit(per_accel,s,0), sp.limit(per_accel,s,sp.oo))
    print("compliant power residual =", power)
    for drive in (sp.Rational(1,50), -sp.Rational(1,50),
                  sp.Rational(1,100), -sp.Rational(1,100)):
        interval = sp.simplify((locked_step-spring_step).subs(A, drive))
        print("A =", drive, "difference between the two required forces =", interval)
    if (power != 0 or locked_step != c**2*A
            or sp.simplify(locked_step-spring_step-c**2*A*sp.cos(omega*t)) != 0):
        raise ArithmeticError("Field constraint or energy transfer failed")


if __name__ == "__main__":
    reaction()
    alternative_reactions()
