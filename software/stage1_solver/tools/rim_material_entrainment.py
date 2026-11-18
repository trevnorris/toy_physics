"""Normal dipole traction on a spherical mouth in a weakly advecting viscous brane."""

import sympy as sp


def check(label, actual, expected):
    residual = sp.factor(actual - expected)
    print(f"{label}: {sp.factor(actual)} | {sp.factor(expected)} | residual {residual}", flush=True)
    if residual != 0:
        raise ArithmeticError(label)


def main():
    r, a, eta, rho, q, gamma = sp.symbols("r a eta rho q Gamma", positive=True)
    w, b = sp.symbols("W b", real=True)
    # In the mouth frame W is the far transverse material speed, b the
    # normal-velocity dipole at the rim; the independent base flow is -q/r^2.
    # Radial velocity is f*cos(theta), polar velocity is -g*sin(theta).
    def polar(f):
        return f + r * sp.diff(f, r) / 2

    def momentum(f, pressure):
        g = polar(f)
        radial = eta * (f.diff(r, 2) + 2*f.diff(r)/r + 4*(g-f)/r**2) - pressure.diff(r)
        theta = eta * (-g.diff(r, 2) - 2*g.diff(r)/r + 2*(g-f)/r**2) + pressure/r
        return sp.factor(radial), sp.factor(theta)

    def shear(f):
        g = polar(f)
        return eta * (-g.diff(r) + (g-f)/r)

    def normal_force(f, pressure):
        return 4*sp.pi*a**2/3 * (2*eta*f.diff(r) - pressure).subs(r, a)

    # Uniform far data and the two decaying homogeneous Stokes dipoles.
    c, h = sp.symbols("c h")
    f0 = w + c/r + h/r**3
    p0 = eta*c/r**2
    at_rim = sp.solve([f0.subs(r, a)-b, shear(f0).subs(r, a)], (c, h), dict=True)[0]
    f0, p0 = f0.subs(at_rim), p0.subs(at_rim)
    for component in momentum(f0, p0):
        check("base transverse momentum residual", component, 0)
    z = 4*sp.pi*eta*a
    force0 = normal_force(f0, p0)
    check("leading mechanical force", force0, z*(w-b))

    # Linearize convective momentum transport in the transverse drive
    # keeping both v0.grad(dv) and dv.grad(v0), in spherical components.
    v0 = -q/r**2
    source_r = sp.factor(rho*(v0*f0.diff(r) + f0*v0.diff(r)))
    source_theta = sp.factor(-rho*v0*(polar(f0).diff(r) + polar(f0)/r))
    print("suction momentum sources (radial, polar) =", source_r, source_theta, flush=True)

    # A 1/r^2 velocity term accommodates the dipolar 1/r^4 source.
    # Allow the pressure dipole to vary independently. Unchanged far and
    # rim velocities at this order leave zero boundary data at fixed b.
    c1, c2, c3, p2, p3 = sp.symbols("c1 c2 c3 p2 p3")
    f1 = c1/r + c2/r**2 + c3/r**3
    p1 = p2/r**2 + p3/r**3
    radial, theta = momentum(f1, p1)
    equations = [f1.subs(r, a), shear(f1).subs(r, a)]
    for residual in (radial-source_r, theta-source_theta):
        equations.extend(sp.Poly(sp.cancel(r**5*residual), r).all_coeffs())
    solution, = sp.solve(equations, (c1, c2, c3, p2, p3), dict=True)
    print("first-order radial profile =", sp.factor(f1.subs(solution)), flush=True)
    print("first-order pressure profile =", sp.factor(p1.subs(solution)), flush=True)
    corrected_f, corrected_p = f1.subs(solution), p1.subs(solution)
    for actual, expected in zip(momentum(corrected_f, corrected_p), (source_r, source_theta)):
        check("first-order exterior momentum", actual, expected)
    check("first-order normal rim velocity", corrected_f.subs(r, a), 0)
    check("first-order tangential traction", shear(corrected_f).subs(r, a), 0)
    mass_rate = 4*sp.pi*rho*q
    force1 = normal_force(corrected_f, corrected_p)
    check("suction correction to mechanical force", force1, mass_rate*(w/6-b/2))

    # gamma is the angularly integrated settled normal contact coefficient.
    # Expand the solved contact law consistently to first order in suction.
    rim_velocity, = sp.solve(gamma*b-force0-force1, b)
    drive_coefficient = sp.factor(gamma * rim_velocity / w)
    leading = drive_coefficient.subs(q, 0)
    suction = sp.factor(q*drive_coefficient.diff(q).subs(q, 0))
    check("settled force per relative speed, leading", leading, gamma*z/(gamma+z))
    check("settled force per relative speed, suction correction", suction,
          mass_rate*gamma*(gamma-2*z)/(6*(gamma+z)**2))
    check("strong-contact leading force per speed", sp.limit(leading, gamma, sp.oo), z)
    check("strong-contact suction correction per speed", sp.limit(suction, gamma, sp.oo), mass_rate/6)
    check("no-contact leading force", leading.subs(gamma, 0), 0)
    check("no-contact suction correction", suction.subs(gamma, 0), 0)


if __name__ == "__main__":
    main()
