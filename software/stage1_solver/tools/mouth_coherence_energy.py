"""Variational sector bounds for a spherical mouth with fixed uniform curvature.

Energies have unit B = 2*pi*R**2*alpha**2/u.  A squared amplitude
has unit alpha/u.  Averages use the uniform probability measure on the
sphere.  The supplied magnetic-Laplacian floor for degree m is |m|/(2R**2).
"""

import sympy as sp


def angular_average(expression, z):
    return sp.factor(sp.integrate(sp.cancel(expression), (z, 0, 1)))


def interval_coefficients(expression, radius, lower, upper):
    """Exact Bernstein coefficients on the closed radius interval."""
    polynomial = sp.Poly(expression, radius)
    degree = polynomial.degree()
    v = sp.Symbol("v")
    pulled = sp.Poly(sum(
        coefficient * (lower + upper*v)**i * (1 + v)**(degree - i)
        for (i,), coefficient in polynomial.terms()
    ), v)
    return [pulled.nth(i) / sp.binomial(degree, i) for i in range(degree + 1)]


def main():
    z = sp.Symbol("z", positive=True)
    t, mass = sp.symbols("t mass", real=True)
    x, y, radius = sp.symbols("x y r", positive=True)
    n, h, s = sp.symbols("n h s", integer=True)
    d = sp.Symbol("d", nonnegative=True)
    checks = []

    # For both helicities take d = min_h |n-C_h|.  The kinetic lower bound
    # is additive in the component masses. Jensen uses their SUM, keeping
    # the cross term in the quartic energy of an arbitrary mixed field.
    total_floor = y*n**2 + mass**2 - 2*(1 - x*d)*mass
    sector_floor = y*n**2 - sp.Max(0, 1 - x*d)**2
    print("two-component mass polynomial:", total_floor)
    print("sector floor over nonnegative mass:", sector_floor)

    degree = n + 2*h*s
    same_s_opposite_n = degree.subs({h: -h, n: -n}, simultaneous=True)
    opposite_s_same_n = degree.subs({h: -h, s: -s}, simultaneous=True)
    opposite_s_opposite_n = degree.subs({s: -s, n: -n}, simultaneous=True)
    print("degrees (original, conjugate+helicity swap, mirror+helicity swap,"
          " mirror+conjugation):", degree, same_s_opposite_n, opposite_s_same_n,
          opposite_s_opposite_n)
    checks.extend([same_s_opposite_n == -degree, opposite_s_same_n == degree,
                   opposite_s_opposite_n == -degree])

    # z=(1+cos(theta))/2. The north potential is (1-z) dphi up to sign.
    # In the south chart the conjugate phase winds once. The factor sqrt(z)
    # is a smooth simple zero in that chart after its transition phase.
    f = sp.sqrt(z) * (1 + t*(1-z))
    sin2 = 4*z*(1-z)
    kinetic_density = sin2/4 * sp.diff(f, z)**2 + (1-z)**2/sin2 * f**2
    gradient = angular_average(kinetic_density, z)
    quadratic = angular_average(f**2, z)
    quartic = angular_average(f**4, z)
    gain = quadratic - 2*x*gradient
    trial = y - gain**2/quartic  # conditional on gain > 0
    print("shape averages (R^2 kinetic, quadratic, quartic):",
          gradient, quadratic, quartic)
    print("unit-mismatch energy polynomial:", y - 2*gain*mass + quartic*mass**2)
    print("optimum squared amplitude and energy for positive gain:",
          sp.factor(gain/quartic), sp.factor(trial))

    # The matched degree-zero component at |n|=2 can be covariantly constant.
    # For every |n|>=2 the local energy floor -1 and n^2>=4 suffice.
    nonunit_floor = 4*y - 1
    neutral_floor = y*n**2 - 1
    print("large-sector excess above |n|=2 from Coulomb alone:",
          sp.factor((y*n**2 - 1) - nonunit_floor))
    print("C_h=0 control (all-sector rough floor, n=0 constant, excess):",
          neutral_floor, total_floor.subs({d: 0, n: 0, mass: 1}),
          sp.factor(neutral_floor + 1))
    print("C_h=+-2, |n|=2 constant / floor:",
          total_floor.subs({d: 0, n: 2, mass: 1}), nonunit_floor)
    checks.append(total_floor.subs({d: 0, n: 2, mass: 1}) == nonunit_floor)

    # r=R/R0.  Stiffness, condensation coefficients and k_e are shared.
    x_r, y_r = 1/(3*radius**2), 2/(9*radius**3)
    law = {x: x_r, y: y_r}
    shape = sp.Rational(2, 5)
    lower, upper = sp.Rational(97, 100), sp.Rational(101, 100)
    print("radius laws, trial shape, radius interval:", x_r, y_r, shape, lower, upper)
    print("fixed shape averages:", *(p.subs(t, shape)
                                    for p in (gradient, quadratic, quartic)))
    unit_upper = sp.factor(trial.subs(t, shape).subs(law))
    zero_floor = -(1 - 2*x_r)**2  # conditional on 1 - 2*x_r > 0
    two_floor = nonunit_floor.subs(law)
    normal_floor = sp.Integer(0)
    print("radius expressions (unit upper, zero lower, |n|>=2 lower):",
          unit_upper, sp.factor(zero_floor), two_floor)
    gaps = [zero_floor - unit_upper, two_floor - unit_upper, normal_floor - unit_upper]
    for point in (lower, sp.Integer(1), upper):
        energies = [p.subs(radius, point) for p in
                    (unit_upper, zero_floor, two_floor, normal_floor)]
        margins = [p.subs(radius, point) for p in gaps]
        print(f"r={point} (x,y):", x_r.subs(radius, point), y_r.subs(radius, point))
        print(f"r={point} (unit upper, zero lower, |n|>=2 lower, normal lower):",
              *energies)
        print(f"r={point} margins (zero, |n|>=2, normal):", *margins)
        checks.extend(p > 0 for p in margins)

    # On the interval positivity of all Bernstein coefficients suffices,
    # since their nonnegative basis polynomials sum to one. Powers of r
    # clear only positive denominators. This includes both interval ends.
    witnesses = (
        ("trial gain", radius**2*gain.subs(t, shape).subs(law)),
        ("zero-sector positive part", radius**2*(1 - 2*x_r)),
        ("gap zero", radius**4*gaps[0]),
        ("gap |n|>=2", radius**4*gaps[1]),
        ("gap normal", radius**4*gaps[2]),
    )
    for label, scaled in witnesses:
        polynomial = sp.cancel(scaled)
        coefficients = interval_coefficients(polynomial, radius, lower, upper)
        print(f"interval {label} (positive-scaled polynomial, min Bernstein):",
              polynomial, min(coefficients))
        checks.append(min(coefficients) > 0)

    radial_base = {x: x_r.subs(radius, 1), y: y_r.subs(radius, 1)}
    unshaped = sp.factor(trial.subs(t, 0).subs(radial_base))
    print("r=1 undeformed trial, zero lower, matched double:",
          unshaped, zero_floor.subs(radius, 1), two_floor.subs(radius, 1))

    # Optional material response: the helicity f retains alpha and its
    # partner has eta*alpha. The original trial is in the retained branch.
    # Weakening a branch raises its energy for EVERY configuration, so all
    # original competing-sector lower bounds remain available.
    eta = sp.Rational(4, 5)
    signed = {}
    for orientation, favored in ((1, 1), (-1, 1), (1, -1), (-1, -1)):
        preferred_c = -2*favored*orientation
        trial_n = sp.sign(preferred_c)
        other_n = -trial_n
        preferred_mismatch = abs(other_n - preferred_c)
        weaker_mismatch = abs(other_n + preferred_c)
        signed[orientation, favored] = trial_n
        print(f"split s={orientation} favored_h={favored}"
              " (trial n, trial preferred mismatch, opposite n,"
              " opposite preferred/weak mismatches):",
              trial_n, abs(trial_n-preferred_c), other_n,
              preferred_mismatch, weaker_mismatch)
        checks.append((abs(trial_n-preferred_c), preferred_mismatch,
                       weaker_mismatch) == (1, 3, 1))
    checks.extend([signed[-1, 1] == -signed[1, 1],
                   signed[-1, -1] == signed[1, 1],
                   signed[1, -1] == -signed[1, 1]])

    strong_gain = 1 - preferred_mismatch*x_r
    weaker_gain = eta - weaker_mismatch*x_r
    opposite_floor = y_r - weaker_gain**2
    split_margin = sp.factor(opposite_floor - unit_upper)
    print("split eta, opposite-sector positive gain candidates:",
          eta, strong_gain, weaker_gain, 0)
    for label, scaled in (
        ("weak minus strong gain", radius**2*(weaker_gain - strong_gain)),
        ("weak gain", radius**2*weaker_gain),
        ("opposite n gap", radius**4*split_margin),
    ):
        polynomial = sp.cancel(scaled)
        coefficients = interval_coefficients(polynomial, radius, lower, upper)
        print(f"split interval {label} (positive-scaled polynomial, min Bernstein):",
              polynomial, min(coefficients))
        checks.append(min(coefficients) > 0)
    for point in (lower, sp.Integer(1), upper):
        wrong_sign = opposite_floor.subs(radius, point)
        margin = split_margin.subs(radius, point)
        print(f"split r={point} (opposite n lower, trial upper, margin):",
              wrong_sign, unit_upper.subs(radius, point), margin)
        checks.append(margin > 0)
    print("split C_h=0 control (favored n=0 constant, rough floor all n):",
          total_floor.subs({d: 0, n: 0, mass: 1}), neutral_floor)

    if not all(checks):
        raise ArithmeticError("An emitted bound, comparison or symmetry failed.")


if __name__ == "__main__":
    main()
