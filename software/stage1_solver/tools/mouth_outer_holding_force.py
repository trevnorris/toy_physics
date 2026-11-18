"""Signed momentum on external cuts of a held, steady two-channel ideal fluid.

The brane and receiver are distinct effective exteriors with a common mass
throughput. The perturbation is the axial dipole of a scalar velocity potential.
"""

import sympy as sp


def cut(dimension):
    radius, density = sp.symbols("R rho", positive=True)
    rate, speed, dipole, load = sp.symbols("q V B p1", real=True)
    angle = sp.symbols("theta", real=True)
    mu = sp.cos(angle)
    perturbation = sp.Function("f")(radius)
    area_unit = 2 * sp.pi**(sp.Rational(dimension, 2)) / sp.gamma(sp.Rational(dimension, 2))
    area = area_unit * radius**(dimension-1)
    angular_weight = sp.sin(angle)**(dimension-2)
    normalization = sp.integrate(angular_weight, (angle, 0, sp.pi))

    def average(value):
        return sp.simplify(sp.integrate(sp.expand(value)*angular_weight,
                                        (angle, 0, sp.pi)) / normalization)

    base_vr = rate/(density*area)
    dvr = mu * perturbation.diff(radius)
    dvx = mu*dvr + (1-mu**2)*perturbation/radius
    adv = average(area*density*base_vr*(dvx+mu*dvr))
    pressure = average(-area*density*base_vr*dvr*mu)
    total = sp.factor(adv + pressure)
    profile = speed*radius + dipole/radius**(dimension-1)
    def evaluate(expr):
        return sp.expand(expr.subs(perturbation, profile).doit())

    # A separate boundary diagnostic for an unspecified pressure harmonic.
    # No solution with this extra pressure is assumed.
    pressure_mode = average(area * load*mu/radius**(dimension-1) * mu)
    print(f"D={dimension} general outward flux: {total}")
    print(f"D={dimension} harmonic [advective, pressure, total]:",
          [evaluate(term) for term in (adv, pressure, total)])
    print(f"D={dimension} pressure-only tail integral: {pressure_mode}")
    return evaluate(total), rate, speed


def main():
    brane, qb, ub = cut(3)
    receiver, qr, ur = cut(4)
    intake = sp.Rational(1, 4)
    for flow in (sp.Rational(1, 50), -sp.Rational(1, 50),
                 sp.Rational(1, 100), -sp.Rational(1, 100)):
        incoming_cut = brane.subs({qb: -intake, ub: flow})
        values = [incoming_cut + receiver.subs({qr: intake, ur: bulk})
                  for bulk in (0, flow/2, flow)]
        print(f"U={flow} holder [W=0, W=U/2, W=U]: {values}")


if __name__ == "__main__":
    main()
