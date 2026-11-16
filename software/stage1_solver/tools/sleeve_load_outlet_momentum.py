"""One throat at 0PN. Vary the confinement energy and evaluate the axial carrier stress."""

import sympy as sp


def check(label, actual, expected):
    residual = sp.factor(actual - expected)
    print(f"{label}: {actual}; residual = {residual}", flush=True)
    if residual != 0:
        raise AssertionError(label)


def confinement():
    # B is the fourth radial moment of the number density at fixed w.
    k, B, R = sp.symbols("k B R", positive=True)
    dB, eta, slope, eps = sp.symbols("dB eta Rp eps", real=True)
    energy = k * B / R**4
    load = -sp.diff(energy, R)
    check("outward confinement load per axial length", load, 4 * k * B / R**5)
    change = load.subs({B: B + eps * dB, R: R + eps * eta})
    check("increment of outward load", sp.factor(sp.diff(change, eps).subs(eps, 0)),
          4 * k * (R * dB - 5 * B * eta) / R**6)
    r, n = sp.symbols("r n", positive=True)
    v_rad = k * (r / R)**4
    body_force = -n * sp.diff(v_rad, R) * slope
    rhs_density = -4 * k * r**4 * n / R**5
    check("axial flare force density plus slope times Path-A source density",
          sp.simplify(body_force + slope * rhs_density), 0)


def outlet():
    w = sp.Symbol("w", real=True)
    m, hb, K = sp.symbols("m hbar K", positive=True)
    j, mu = sp.symbols("j mu", real=True)
    a = sp.Function("a", positive=True)(w)
    V = sp.Function("V")(w)
    n = a**2
    stress_q = -hb**2 * n * sp.diff(sp.log(n), w, 2) / (4 * m)
    flux = m * j**2 / n + K * n**5 + stress_q
    # Stationary, smooth, ungauged, strictly axial carrier. The current is constant.
    amplitude_ode = (2 * m / hb**2) * (
        (V + 5 * K * n**4 / 4 - mu) * a + m * j**2 / (2 * a**3)
    )
    on_shell = sp.expand(flux.subs(a.diff(w, 2), amplitude_ode))
    integral = hb**2 * a.diff(w)**2 / (2 * m) + m * j**2 / (2 * n) + (mu - V) * n - K * n**5 / 4
    check("stationary axial flux", on_shell, integral)
    check("axial flux gradient plus n V'", sp.simplify(
        sp.diff(integral, w).subs(a.diff(w, 2), amplitude_ode) + n * V.diff(w)), 0)

    # Local regular data with common density, mu and scalar pressure. No global throat is assumed.
    settings = {m: 1, hb: 1, K: sp.Rational(4, 5), mu: 1, V: 0, a: 1}
    pb = sp.Rational(4, 5)
    for current, grad, target in [
        (0, 0, 0), (sp.Rational(1, 5), 0, sp.Rational(1, 50)),
        (0, sp.Rational(1, 5), sp.Rational(1, 50))
    ]:
        result = integral.subs({j: current, a.diff(w): grad}).subs(settings) - pb
        check(f"excess per area at (j,a')=({current},{grad})", result, target)

    eps = sp.Symbol("eps")
    n0, dn, du = sp.symbols("n0 dn du", real=True)
    drift = m * (n0 + eps * dn) * (eps * du)**2
    check("convective linear coefficient at rest", drift.expand().coeff(eps, 1), 0)
    check("convective quadratic coefficient at rest", drift.expand().coeff(eps, 2), m * n0 * du**2)


if __name__ == "__main__":
    confinement()
    outlet()
