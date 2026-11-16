"""Companion-directed outward rates on two fixed flat face annuli.

Use the uniform face supply of drain_on_replenishing_slab.py. Superpose two
prescribed inward monopoles in the in-plane constituent current and assign
transferred material the local material velocity. The native core is excluded.
"""
import sympy as sp


def main():
    a, b, d, r, distance = sp.symbols("a b d r t", positive=True)
    h, n, ell, mass = sp.symbols("h n ell m", positive=True)
    q_target, q_source = sp.symbols("q_target q_source", positive=True)
    background_gradient = sp.Symbol("gamma", real=True)
    mu, w = sp.symbols("mu w", real=True)

    # At fixed r < d, distance from the companion ranges from d-r to d+r.
    cosine = (d**2 + r**2 - distance**2) / (2*d*r)
    radial_speed = background_gradient*r - q_target/(4*sp.pi*h*n*r**2)
    control_axial = radial_speed*mu
    companion_axial = q_source*(d-r*mu)/(4*sp.pi*h*n*distance**3)
    measure = 2*sp.pi*r*distance/d

    def face_flow(axial_speed, mass_flux):
        integrand = sp.expand(mass_flux*axial_speed.subs(mu, cosine)*measure)
        primitive = sp.integrate(integrand, distance)
        shell = sp.factor(primitive.subs(distance, d+r) - primitive.subs(distance, d-r))
        return sp.factor(sp.integrate(shell, (r, a, b)))

    number_current_w = -ell*w
    normal_pressure = sp.Symbol("p", real=True)
    direction = sp.Matrix([1, 0, 0, 0])
    sample = {a: 1, b: 3, d: 10, h: 2, n: 3, ell: sp.Rational(3, 25),
              mass: 1, q_source: 24*sp.pi}
    for label, sign in (("upper", 1), ("lower", -1)):
        normal = sp.Matrix([0, 0, 0, sign])
        mass_flux = sign*mass*number_current_w.subs(w, sign*h/2)
        isolated = face_flow(control_axial, mass_flux)
        paired = face_flow(control_axial + companion_axial, mass_flux)
        delta = sp.factor(paired - isolated)
        pressure_projection = direction.dot(normal_pressure*normal)
        print(label, "outward mass flux per true area =", mass_flux)
        print(label, "flow [isolated, paired, change] =", [isolated, paired, delta])
        print(label, "sample change in flow =", delta.subs(sample))
        print(label, "outward normal stress along axis per area =", pressure_projection)

    # The spherical intake probe uses volume throughput rather than number rate.
    rho, volume_rate, k, incident = sp.symbols("rho Q k U", real=True)
    outlet = rho*volume_rate*k*incident
    for fraction in (0, 1):
        scale = {rho: 3/(4*sp.pi), volume_rate: 4*sp.pi, k: fraction}
        control = outlet.subs(scale).subs(incident, 0)
        driven = outlet.subs(scale).subs(incident, sp.Rational(1, 100))
        print("probe k =", fraction, "outlet flow [control, driven, change] =",
              [control, driven, driven-control])


if __name__ == "__main__":
    main()
