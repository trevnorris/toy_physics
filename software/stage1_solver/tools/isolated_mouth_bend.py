"""Compute the normal load and midpoint tail from the two outer fluid cuts."""

import sympy as sp


def check(label, actual, expected):
    print(f"{label} = {actual}", flush=True)
    if sp.simplify(actual - expected) != 0:
        raise ArithmeticError(label)


def main():
    r, a, k, D, rho, rho_f, K = sp.symbols("r a k D rho rho_f K", positive=True)
    z, A, B, e, T = sp.symbols("z A B e T", real=True)
    # z runs from lower to upper face, H=pi/k. Only the zero mode carries net mass.
    feed = (A + B * sp.exp(-k*r) * sp.cos(k*z)) / r
    vr, vw = sp.diff(feed, r), sp.diff(feed, z)
    def section(integrand):
        return sp.simplify(4 * sp.pi * r**2 * sp.integrate(integrand, (z, 0, sp.pi/k)))

    check("inward feed mass rate", section(-rho_f * vr), 4*sp.pi**2*rho_f*A/k)
    feed_flux = section(rho_f * vr * vw)
    check("feed outward w flux", feed_flux, 8*sp.pi*rho_f*A*B*sp.exp(-k*r)/r)

    # On the 4D hemisphere t=sin(theta)^2. After taking out cos(theta),
    # the axial stress is polynomial in t, and the angular weight is sqrt(t)/2.
    t = sp.symbols("t", nonnegative=True)
    v0 = D / (sp.pi**2 * rho * r**3)
    shape = e*(a/r)**2
    radial = v0*(1 + shape*(3-4*t))
    polar_factor = 2*v0*shape  # v_theta = polar_factor * sin(theta)*cos(theta)
    axial_over_cos = radial - polar_factor*t
    speed2 = radial**2 + polar_factor**2*t*(1-t)
    stress_over_cos = rho*(radial*axial_over_cos - speed2/2)
    bulk_flux = sp.expand(2*sp.pi*r**3*sp.integrate(
        sp.expand(stress_over_cos)*sp.sqrt(t), (t, 0, 1)))
    prefactor = 2*D**2/(sp.pi**3*rho)
    check("receiver outward w flux", bulk_flux,
          prefactor*(1/(3*r**3) - 2*e*a**2/(5*r**5) + e**2*a**4/(7*r**7)))

    # A single material budget encloses both cuts and all compact contacts.
    # Static, no independent normal support, all return acts on the midpoint.
    load = -bulk_flux-feed_flux
    check("enclosed common load at infinity", sp.limit(load, r, sp.oo), 0)
    bulk_height = sp.expand(sp.integrate(bulk_flux/(4*sp.pi*K*r**2), r))
    check("receiver part of midpoint displacement", bulk_height,
          -prefactor/(4*sp.pi*K)*(1/(12*r**4) - e*a**2/(15*r**6) + e**2*a**4/(56*r**8)))
    check("receiver height at infinity", sp.limit(bulk_height, r, sp.oo), 0)
    # m(infinity)=0, m'=-load/(4*pi*K*r^2). Use l'Hopital for r*m and r^4*m.
    slope = -load/(4*sp.pi*K*r**2)
    check("coefficient of 1/r in full displacement", sp.limit(-r**2*slope, r, sp.oo), 0)
    fourth = sp.limit(-r**5*slope/4, r, sp.oo)
    check("limit of r^4 times full displacement", fourth, -D**2/(24*sp.pi**4*rho*K))
    # Positive control: a separate sector exports signed w momentum T per second.
    check("1/r coefficient with independent outward transfer T",
          sp.limit((load-T)/(4*sp.pi*K), r, sp.oo), -T/(4*sp.pi*K))


if __name__ == "__main__":
    main()
