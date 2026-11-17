"""Measure the odd-in-speed energy signal on a fixed cut, separately from static force."""
from math import cos, fsum, gamma, hypot, pi, sin


def area(n):
    return 2 * pi ** (n / 2) / gamma(n / 2)


def moving_state(n, rho, q, x, y):
    r = hypot(x, y)
    scale = q / (rho * area(n) * r ** n)
    return scale * x, scale * y


def surface(n, radius, offset, density):
    steps = 1024
    terms = []
    for k in range(steps + 1):
        angle = pi * k / steps
        nx, ny = cos(angle), sin(angle)
        scale = 1 if k in (0, steps) else (4 if k % 2 else 2)
        terms.append(scale * ny ** (n - 2) * density(offset + radius * nx, radius * ny, nx, ny))
    return area(n - 1) * radius ** (n - 1) * pi * fsum(terms) / (3 * steps)


def subtract(observable, q1, q2):
    return fsum((observable(q1, q2), -observable(q1, 0.), -observable(0., q2), observable(0., 0.)))


def power(n, rho, q1, q2, d, radius, s, include_unsteady=True, offset=0.):
    def state(a, b):
        def density(x, y, nx, ny):
            ax, ay = moving_state(n, rho, a, x, y)
            bx, by = moving_state(n, rho, b, x - d, y)
            ux, uy = ax + bx, ay + by
            kinetic = rho * (ux ** 2 + uy ** 2) / 2
            pressure = 0.59 - kinetic
            if include_unsteady:
                pressure -= rho * s * bx  # partial_t phi_2 = s * u_2x at V_2=-s
            return (0.33 + kinetic + pressure) * (ux * nx + uy * ny)
        return surface(n, radius, offset, density)
    return subtract(state, q1, q2)


def static_load(n, rho, q1, q2, d):
    def state(a, b):
        def density(x, y, nx, ny):
            ax, ay = moving_state(n, rho, a, x, y)
            bx, by = moving_state(n, rho, b, x - d, y)
            ux, uy = ax + bx, ay + by
            pressure = 0.59 - rho * (ux * ux + uy * uy) / 2
            return -pressure * nx - rho * ux * (ux * nx + uy * ny)
        return surface(n, 0.87, 0., density)
    return subtract(state, q1, q2)


def main():
    s = 0.2
    passed, lost_unsteady, shape_changes = [], [], []
    print("n radius      J(0)        J(s)         J(-s)        J(s/2)      J(-s/2)      slope_full   slope_half   rest_load    no_unsteady  offset_slope")
    for n, rho, q1, q2 in [(3, 1.4, -1.1, -1.9), (4, 0.85, 1.1, 1.9)]:
        d = 3.25
        load = static_load(n, rho, q1, q2, d)
        for radius in (0.7, 1.25):
            args = (n, rho, q1, q2, d, radius)
            j0, jp, jm, jhp, jhm = (power(*args, v) for v in (0., s, -s, s/2, -s/2))
            full, half = (jp-jm)/(2*s), (jhp-jhm)/s
            omitted = (power(*args, s, False)-power(*args, -s, False))/(2*s)
            # Same fixed, enclosing sphere with its center moved slightly toward mouth 2.
            offset = (power(*args, s, offset=0.25)-power(*args, -s, offset=0.25))/(2*s)
            print(f"{n} {radius:5.2f} " + " ".join(f"{val: .6e}" for val in
                (j0, jp, jm, jhp, jhm, full, half, load, omitted, offset)))
            passed.extend((abs(j0), abs(full-load), abs(half-load), abs(jp+jm), abs(jhp+jhm)))
            lost_unsteady.append(abs(omitted-load))
            shape_changes.append(abs(offset-load))
    print("maximum centered residual = %.3e" % max(passed))
    print("minimum slope error with unsteady pressure removed = %.3e" % min(lost_unsteady))
    print("minimum slope change for a +0.25 center displacement = %.3e" % min(shape_changes))
    if max(passed) > 1.e-10 or min(lost_unsteady) < 1.e-4 or min(shape_changes) < 1.e-4:
        raise ArithmeticError("centered-surface or negative-control check")


if __name__ == "__main__":
    main()
