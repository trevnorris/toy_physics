"""First variation of mass and axial momentum on a fixed 4D exterior cut."""

import sympy as sp


def angular_mean(expression, mu):
    """Average over the ordinary two-sphere of directions within the brane."""
    return sp.expand(sp.integrate(expression, (mu, -1, 1)) / 2)


def main():
    mu, drive, s = sp.symbols("mu U s", real=True)
    j0, j1, vr, a, b, p0, p1 = sp.symbols("j0 j1 vr A B p0 p1")
    j = j0 + drive * j1 * mu
    vx = vr * mu + drive * (a * mu**2 + b * (1 - mu**2))
    pressure = p0 + drive * p1 * mu
    normal_x = s * mu

    def linear_mean(expression):
        return angular_mean(sp.diff(expression, drive).subs(drive, 0), mu)

    mass = linear_mean(j)
    carried = linear_mean(j * vx)
    load = linear_mean(pressure * normal_x)
    print("mass angular mean =", mass)
    print("carried angular mean =", carried)
    print("pressure angular mean =", load)
    assert mass == 0
    assert sp.expand(3 * carried - j1 * vr - j0 * (a + 2 * b)) == 0
    assert 3 * load == s * p1

    # Insert a radial 4D outflow and the gradient of f(R) N_x as a control.
    theta, radius, density = sp.symbols("theta R rho", positive=True)
    q, h, t = sp.symbols("q h t")  # h = f'(R), t = f(R)/R
    sin_theta = sp.sin(theta)
    v0 = q / (2 * sp.pi**2 * radius**3 * density)
    radial = {
        s: sin_theta, j0: density * v0, vr: v0 * sin_theta,
        j1: density * h * sin_theta, a: t + (h - t) * sin_theta**2,
        b: t, p1: -density * v0 * h * sin_theta,
    }
    area = 4 * sp.pi * radius**3 * sin_theta**2
    flux = sp.simplify(sp.integrate(area * (carried + load).subs(radial), (theta, 0, sp.pi)))
    drift, dipole = sp.symbols("V L")
    f = drift * radius + dipole / radius**3
    potential_flux = sp.simplify(flux.subs({h: sp.diff(f, radius), t: f / radius}))
    print("radial Bernoulli control =", flux)
    print("uniform plus decaying dipole =", potential_flux)
    assert sp.expand(flux - q * (h + 3 * t) / 4) == 0
    assert potential_flux == q * drift


if __name__ == "__main__":
    main()
