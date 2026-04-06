from mpmath import exp, pi, quad


def radius_profile(w, a):
    return a * (1.0 + 0.5 * exp(-5.0 * w / a))


def main():
    a = 1.0
    L = 2.0
    rho0 = 1.0
    eps = 0.1

    i3 = quad(lambda w: radius_profile(w, a) ** 3, [0.0, L])
    i5 = quad(lambda w: radius_profile(w, a) ** 5, [0.0, L])

    mass = (4.0 * pi * rho0 / 3.0) * i3
    quadrupole = (4.0 * pi * rho0 * eps / 25.0) * i5
    alpha2 = quadrupole / (mass * a**2)

    print("Rounded funnel parameters:")
    print(f"a = {a}, L = {L}, rho0 = {rho0}, eps = {eps}")
    print(f"M = {float(mass):.12f}")
    print(f"Q = {float(quadrupole):.12f}")
    print(f"alpha_2 = {float(alpha2):.12f}")


if __name__ == "__main__":
    main()
