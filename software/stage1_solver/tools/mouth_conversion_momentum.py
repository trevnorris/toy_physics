"""Signed material momentum through a moving sphere in the draining exterior."""

import sympy as sp


def check(label, actual, expected):
    residual = sp.factor(actual - expected)
    print(f"{label}: {sp.factor(actual)} | {sp.factor(expected)} | residual {residual}", flush=True)
    if residual != 0:
        raise ArithmeticError(label)


def main():
    a, r, q, rho, eta = sp.symbols("a r q rho eta", positive=True)
    mu, w, b, v, shift = sp.symbols("mu W b V shift", real=True)
    drive = sp.Symbol("drive")
    rate = 4 * sp.pi * rho * q
    base = -q / r**2

    def surface(expr):
        return sp.factor(2 * sp.pi * r**2 * sp.integrate(expr, (mu, -1, 1)))

    def polar(f):
        return f + r * f.diff(r) / 2

    # Supplied material profiles at fixed b from rim_material_entrainment.
    f = w + a*(b-w)/r
    g = polar(f)
    pressure = eta*a*(b-w)/r**2
    f1 = -q*rho*(b-w)*(r-a)*(2*r-a)/(4*eta*r**3)
    p1 = -q*rho*(3*w*a-3*w*r-3*a*b+b*r)/(2*r**3)

    # Relative velocity sets the throughput. Absolute velocity sets the
    # lab momentum. O(q*drive) transport needs only the leading dipole.
    inward_mass = -rho*(base + drive*f*mu)
    axial_lab = base*mu + drive*(v + f*mu**2 + g*(1-mu**2))
    incoming = surface((inward_mass*axial_lab).diff(drive).subs(drive, 0))
    redistributed = surface(inward_mass.diff(drive)*axial_lab.subs(drive, 0))
    sped_up = surface(inward_mass.subs(drive, 0)*axial_lab.diff(drive))
    check("mass intake", surface(inward_mass), rate)
    check("undriven axial inflow", surface((inward_mass*axial_lab).subs(drive, 0)), 0)
    check("rim mass redistribution contribution", redistributed.subs(r, a), rate*b/3)
    check("rim velocity contribution", sped_up.subs(r, a), rate*(v+w/3+2*b/3))
    check("rim lab axial momentum rate", incoming.subs(r, a), rate*(v+w/3+b))
    check("held strong-contact carriage", incoming.subs({r:a, v:0, b:0}), rate*w/3)
    check("lab shift of carried rate", incoming.subs(v, v+shift)-incoming, rate*shift)

    def stress(f_r, pressure_r):
        g_r = polar(f_r)
        shear = eta*(-g_r.diff(r)+(g_r-f_r)/r)
        return surface((2*eta*f_r.diff(r)-pressure_r)*mu**2-shear*(1-mu**2))

    mechanical = stress(f, pressure) + stress(f1, p1)
    print("arbitrary-radius mechanical rate =", sp.collect(sp.expand(mechanical), q), flush=True)
    print("arbitrary-radius carried rate =", incoming, flush=True)
    net = sp.factor(mechanical + incoming)
    check("net flux independent of radius", net, net.subs(r, a))
    check("combined inward relative momentum rate", net-rate*v,
          4*sp.pi*eta*a*(w-b) + rate*(w+b)/2)


if __name__ == "__main__":
    main()
