"""Axisymmetric quadrature of primitive Bernoulli energy and momentum fluxes."""
from math import cos, fsum, gamma, hypot, pi, sin


def area(n):
    return 2 * pi ** (n / 2) / gamma(n / 2)


def field(n, rho, q, dx, y):
    r = hypot(dx, y)
    scale = q / (rho * area(n) * r ** n)
    return scale * dx, scale * y


def integrate(n, rho, rates, d, cut, radius, speed, ref_pressure=0.71, ref_energy=0.29):
    """Four states on one geometric sphere. Normals point away from its center."""
    count = 1024
    integrated = [[] for _ in range(3)]
    for k in range(count + 1):
        theta = pi * k / count
        nx, ny = cos(theta), sin(theta)
        x, y = cut * d + radius * nx, radius * ny
        states = [[] for _ in range(3)]
        for active1, active2, sign in [(1, 1, 1), (1, 0, -1), (0, 1, -1), (0, 0, 1)]:
            ax, ay = field(n, rho, active1 * rates[0], x, y)
            bx, by = field(n, rho, active2 * rates[1], x - d, y)
            ux, uy = ax + bx, ay + by
            normal_speed = ux * nx + uy * ny
            kinetic = rho * (ux * ux + uy * uy) / 2
            energy = kinetic + ref_energy
            pressure = ref_pressure - kinetic + rho * speed * bx
            lab_current = (energy + pressure) * normal_speed
            reynolds = energy * cut * speed * nx
            static_pressure = ref_pressure - kinetic
            static_force = -(static_pressure * nx + rho * ux * normal_speed)
            for values, value in zip(states, (lab_current - reynolds, lab_current, static_force)):
                values.append(sign * value)
        coeff = 1 if k in (0, count) else 2 if k % 2 == 0 else 4
        weight = coeff * ny ** (n - 2)
        for values, terms in zip(integrated, states):
            values.append(weight * fsum(terms))
    measure = area(n - 1) * radius ** (n - 1) * pi / (3 * count)
    return tuple(fsum(values) * measure for values in integrated)


def check():
    d = 3.25
    max_error = 0.
    print("n   r1   r2     V2         J1        J2_lab    J2_relative    F1_static    F2_static     P2_center")
    for n, rho, rates in [(3, 1.4, (-1.1, -1.9)), (4, 0.85, (1.1, 1.9))]:
        amplitude = rates[0] * rates[1] / (rho * area(n) * d ** (n - 1))
        for r1, r2 in [(0.7, 1.15), (1.25, 0.5)]:
            for speed in (-0.2, 0.2):
                j1, _, f1 = integrate(n, rho, rates, d, 0, r1, speed)
                j2, lab2, f2 = integrate(n, rho, rates, d, 1, r2, speed)
                p2 = speed * f2
                print(f"{n} {r1:4.2f} {r2:4.2f} {speed:6.2f} {j1: .9e} {lab2: .9e} {j2: .3e} "
                      f"{f1: .9e} {f2: .9e} {p2: .9e}")
                errors = [j1 + speed * amplitude, j2, f1 - amplitude, f2 + amplitude, j1 - p2]
                max_error = max(max_error, *(abs(value) for value in errors))
    print(f"maximum absolute error across energy and static-force comparisons = {max_error:.3e}")
    if max_error >= 1.e-10:
        raise ArithmeticError("surface integrals")


if __name__ == "__main__":
    check()
