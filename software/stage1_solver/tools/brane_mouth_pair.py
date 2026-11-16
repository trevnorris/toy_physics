"""Integrate the momentum removed at a flat mouth in a reduced 3-D potential flow.

The positive sink density is the volume rate removed per brane volume. Each
radial profile is supported on its own mouth and is unchanged in the pair.
The removed material carries the local brane velocity. Ordinary pressure is
the only stress in this ideal calculation. The physical mouth normal points
along w, while the momentum component of interest points along x.
"""

import sympy as sp


def main():
    r, s, d, a, b = sp.symbols("r s d a b", positive=True)
    rho, qt, qs = sp.symbols("rho Q_t Q_s", positive=True)
    mu = sp.Symbol("mu", real=True)
    pressure = sp.Symbol("pressure", real=True)
    source_sink = sp.Function("s_source")
    target_sink = sp.Function("s_target")

    # The field point is a distance R from the source centre. Integrate the
    # Green function over a source shell of radius s < R (d > a + b).
    R, distance = sp.symbols("R ell", positive=True)
    kernel = 1 / (4 * sp.pi * distance)
    angular_measure = 2 * sp.pi * distance / (R * s)
    source_shell = sp.integrate(kernel * angular_measure, (distance, R - s, R + s))
    source_rate = sp.Integral(4 * sp.pi * s**2 * source_sink(s), (s, 0, b))
    # Each shell has the same angular mean. Its volume rate is therefore
    # the weight in this average, with Q_s defined by source_rate.
    potential = sp.factor(source_shell * qs / (4 * sp.pi))

    # On a target shell, R runs from d-r to d+r. The x component of
    # grad R in this coordinate frame is (r*mu-d)/R.
    cosine = (d**2 + r**2 - R**2) / (2 * d * r)
    velocity_x = sp.diff(potential, R) * (r * cosine - d) / R
    angle_mean = sp.simplify(sp.integrate(velocity_x * R / (2 * d * r), (R, d-r, d+r)))
    own_radial_velocity = sp.Function("v_self")(r)
    own_mean = sp.integrate(own_radial_velocity * mu / 2, (mu, -1, 1))
    target_rate = sp.Integral(4 * sp.pi * r**2 * target_sink(r), (r, 0, a))
    carried_alone = rho * target_rate * own_mean
    carried_pair = (rho * target_rate * (own_mean + angle_mean)).subs(target_rate, qt)

    axis_x, normal_w = sp.eye(4).col(0), sp.eye(4).col(3)
    pressure_x = (axis_x.T * (pressure * sp.eye(4)) * normal_w)[0]
    delta = sp.factor(carried_pair - carried_alone + pressure_x)
    print("source shell angular integral of 1/(4*pi*distance) =", sp.factor(source_shell))
    print("Q_s =", source_rate)
    print("potential outside finite source =", potential)
    print("target shell mean [alone, companion velocity] =", [own_mean, angle_mean])
    print("carried x momentum per time [alone, pair] =", [carried_alone, carried_pair])
    print("isotropic pressure x projection at mouth =", pressure_x)
    print("paired minus alone mouth load =", delta)
    scaled = sp.factor(4 * sp.pi * delta / (rho * qt * qs))
    print("scaled load at d = 4, 8 =", [scaled.subs(d, x) for x in (4, 8)])
    print("zero companion, doubled separation / original =",
          delta.subs(qs, 0), sp.simplify(delta.subs(d, 2*d) / delta))
    if own_mean != 0 or delta.subs(qs, 0) != 0 or delta.subs(qt, 0) != 0:
        raise ArithmeticError("A load survives in the isolated radial aperture control")


if __name__ == "__main__":
    main()
