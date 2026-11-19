"""Charge an imposed test distribution using the brane Maxwell zero mode."""

import sympy as sp


def shell_cost():
    r, a, b = sp.symbols("r a b", positive=True)
    # Normalize by e^2/(8 pi epsilon). Gauss flux is the enclosed fraction.
    density = 3 / (4 * sp.pi * (b**3 - a**3))
    fraction = sp.factor(sp.integrate(4 * sp.pi * r**2 * density, (r, a, r)))
    numerator, denominator = sp.fraction(fraction)
    integrand = sp.expand(numerator**2 / r**2)
    inside = sp.integrate(integrand, (r, a, b)) / denominator**2
    tail = sp.integrate(1 / r**2, (r, b, sp.oo))
    cost = sp.factor(inside + tail)
    print("shell enclosed fraction =", fraction)
    print("tail coefficient =", tail)
    print("one-unit shell coefficient =", cost)
    assert sp.simplify(fraction - (r**3 - a**3)/(b**3 - a**3)) == 0
    assert sp.simplify(sp.limit(cost, b, a) - 1/a) == 0
    return a, b, cost


def additions(a, b, cost):
    width = sp.Rational(1, 20)
    radii = [sp.Rational(1, 2), 1, 2]
    windows = []
    for radius in radii:
        coefficient = cost.subs({a: radius, b: radius + width})
        for direction in (-1, 1):
            energies = [coefficient * (direction * n)**2 for n in range(3)]
            steps = [right - left for left, right in zip(energies, energies[1:])]
            print(f"R={radius} step sign={direction:+} first, second =", *steps)
        windows.append((coefficient, 3 * coefficient))
    lower = max(lo for lo, hi in windows)
    upper = min(hi for lo, hi in windows)
    print("constant linear reward across radii needs greater than, less than =", lower, upper)
    print("open intersection exists =", bool(lower < upper))


def reflection_test():
    # At any fixed in-brane wave number, one real coefficient of each depth mode
    # suffices. Ge and Go stand for the positive static Green weights.
    n, ge, go, be, bo, te, to = sp.symbols("n Ge Go Be Bo Te To", real=True)
    def stored(side, units):
        return (ge * (be + units * te)**2 + go * (side * bo + units * side * to)**2) / 2

    for side in (1, -1):
        print(f"side {side:+} field cost with common background =",
              sp.expand(stored(side, n) - stored(side, 0)))
    def odd_difference(extra_reward):
        def energy(s, units):
            return stored(s, units) - s * extra_reward * units
        return sp.expand((energy(1, n) - energy(1, -n) - energy(-1, n) + energy(-1, -n)) / 4)

    reward = sp.symbols("lambda", real=True)
    print("mirror/sign double difference, written field =", odd_difference(0))
    print("mirror/sign double difference, inserted reward =", odd_difference(reward))
    assert odd_difference(0) == 0 and odd_difference(reward) == -reward * n


if __name__ == "__main__":
    additions(*shell_cost())
    reflection_test()
