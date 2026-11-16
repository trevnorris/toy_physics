"""Reversal and force-probe experiment on a spherical, uniformly draining intake.

The ideal-fluid exterior is three-dimensional. The outlet and mechanical mass
are selectable benchmark inputs. Forces point along the imposed flow axis.
"""
import json

import sympy as s


def main():
    a, r, rho, Q, mass = s.symbols("a r rho Q M", positive=True)
    mu, U, V, A, force, hold, k = s.symbols("mu U V A F H k", real=True)
    monopole, dipole = s.symbols("monopole dipole")
    # At the instant X=0, coordinates are centered on the mouth. Determine the
    # l=0,1 exterior harmonics from the normal speed relative to the surface.
    potential = U*r*mu + monopole/r + dipole*mu/r**2
    normal_error = s.diff(potential, r).subs(r, a) - V*mu + Q/(4*s.pi*a**2)
    harmonics = s.solve(s.Poly(normal_error, mu).all_coeffs(), (monopole, dipole))
    phi = potential.subs(harmonics)
    radial = s.diff(phi, r)
    angular = s.diff(phi, mu)/r
    axial = radial*mu + angular*(1-mu**2)
    speed_sq = radial**2 + (1-mu**2)*angular**2
    # Bernoulli uses the time derivative at fixed lab position. Here V=0
    # at the observation, so the coordinate shift contributes zero, even for A!=0.
    pressure = -rho*(s.diff(phi, V)*A + (speed_sq-U**2)/2)

    def inward(integrand):
        at_rest = integrand.subs(V, 0).expand()
        return s.factor(-2*s.pi*r**2*s.integrate(at_rest, (mu, -1, 1)))

    net_pressure = inward(pressure*mu)
    net_carried = inward(rho*axial*radial)
    incoming = s.factor(net_pressure + net_carried)
    # Same control volume and outlet for holding and probe. Constant inventory,
    # instantaneous outlet mean speed V+k(U-V), zero additional outlet traction.
    outlet = rho*Q*(V+k*(U-V))
    balance = mass*A - (incoming - outlet.subs(V, 0) + force + hold)
    holding = s.factor(s.solve(balance.subs({A: 0, force: 0}), hold)[0])
    acceleration = s.factor(s.solve(balance.subs(hold, 0), A)[0])

    print("normal_bc_residual =", s.simplify(normal_error.subs(harmonics)))
    print("potential =", phi)
    print("inward_pressure_at_r =", net_pressure)
    print("inward_carried_at_r =", net_carried)
    print("inward_total_at_r =", incoming)
    print("external_hold =", holding)
    print("released_acceleration =", acceleration)
    print("outward_outlet_x_at_release =", outlet.subs(V, 0))

    # Units: a, suction speed Q/(4 pi a^2), displaced mass. These values
    # set units, not native throat properties. The benchmark body has M=1.
    units = {a: 1, rho: 3/(4*s.pi), Q: 4*s.pi, mass: 1}
    probe = s.Rational(3, 200)
    measurements = {}
    for fraction in (0, 1):
        for speed in (s.Rational(1, 100), 0, s.Rational(-1, 100)):
            state = {**units, k: fraction, U: speed}
            key = f"k={fraction},U={speed}"
            h = holding.subs(state)
            released = [acceleration.subs({**state, force: f}) for f in (-probe, 0, probe)]
            measurements[key] = [str(s.cancel(v)) for v in (h, *released)]
    print("measurements =", json.dumps(measurements, sort_keys=True))


if __name__ == "__main__":
    main()
