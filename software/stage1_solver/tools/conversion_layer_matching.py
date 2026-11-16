"""Test one-throat interface balances and the conditional steady order-front integral."""

import sympy as sp


def verify(label, expression):
    result = sp.factor(expression)
    print(f"{label} = {result}", flush=True)
    if result != 0:
        raise AssertionError(label)


def mechanical_match():
    m, n_minus, n_plus = sp.symbols("m n_minus n_plus", positive=True)
    speed, j = sp.symbols("W j", real=True)
    t_minus, t_plus = sp.symbols("T_minus T_plus", real=True)

    def normal_flux(n, t):
        u = speed + j / n
        # Stress and advective momentum minus the swept momentum at the moving cut.
        return t + m * n * u**2 - speed * m * n * u

    reaction = sp.factor(normal_flux(n_minus, t_minus) - normal_flux(n_plus, t_plus))
    print("carrier force on layer per area =", reaction, flush=True)
    verify("normal momentum budget",
           reaction - (t_minus - t_plus + m * j**2 * (1 / n_minus - 1 / n_plus)))
    verify("equal endpoint densities", reaction.subs(n_plus, n_minus) - (t_minus - t_plus))


def order_match():
    z = sp.Symbol("z", real=True)
    q = sp.Function("chi")(z)
    area, density, mobility = [sp.Function(name, positive=True)(z) for name in ("A", "n", "M")]
    potential = sp.Function("f")(q)
    kappa, ell = sp.symbols("kappa ell", positive=True)
    current, extra = sp.symbols("I mu_extra", real=True)
    slope = q.diff(z)
    free_energy = area * (potential + kappa * slope**2 / 2)
    mu_bare = (free_energy.diff(q) - free_energy.diff(slope).diff(z)) / area
    H = area.diff(z) / area
    velocity = current / (area * density)
    law = velocity * slope / mobility + mu_bare + extra
    boundary = potential - kappa * slope**2 / 2
    integral_density = current * slope**2 / (area * density * mobility) + boundary.diff(z) - kappa * H * slope**2 + extra * slope
    verify("order equation times slope minus integral density", law * slope - integral_density)

    # Independently evaluate the nonzero planar weight on the admitted decreasing kink.
    y = sp.Symbol("y", real=True)
    planar_dq_dz = -y * (1 - y) / ell
    weight = sp.integrate(planar_dq_dz, (y, 1, 0))
    print("planar kink integral of slope squared =", weight, flush=True)
    verify("planar kink weight times ell minus 1/6", ell * weight - sp.Rational(1, 6))
    a = sp.Symbol("a", positive=True)
    well = a * y**2 * (1 - y)**2
    verify("planar homogeneous endpoints, energy jump", well.subs(y, 0) - well.subs(y, 1))


if __name__ == "__main__":
    mechanical_match()
    order_match()
