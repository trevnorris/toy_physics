"""Three tests of orientation selection. The biased terms are detection controls."""

from collections import Counter
from itertools import product

import sympy as sp


def show(label, value):
    print(f"{label} = {value}", flush=True)


def standing_wave():
    t, z = sp.symbols("t z", real=True)
    length, speed, tension, amplitude = sp.symbols("L c T A", positive=True)
    k = sp.pi / (2 * length)
    frequency = speed * k
    for side, spin in product((1, -1), repeat=2):
        # z is the lab depth. The two tangent components use fixed lab axes.
        phase = frequency * t
        field = amplitude * sp.sin(k * side * z) * sp.Matrix(
            [sp.cos(phase), spin * sp.sin(phase)]
        )
        slope = field.diff(z)
        assert field.subs(z, 0) == sp.zeros(2, 1)
        assert slope.subs(z, side * length) == sp.zeros(2, 1)
        assert sp.simplify(field.diff(t, 2) - speed**2 * field.diff(z, 2)) == sp.zeros(2, 1)
        rim = tension * slope.subs(z, 0)
        load_sq = sp.trigsimp(rim.dot(rim))
        flux = -tension * field.diff(t).dot(slope)
        mean_flux = sp.simplify(sp.integrate(flux, (t, 0, 2 * sp.pi / frequency)) * frequency / (2 * sp.pi))
        show(f"side {side:+} polarization {spin:+} (frequency, squared rim traction, mean power flux)",
             (frequency, load_sq, mean_flux))

    # A general reciprocal elastic end may mix the circular components.
    a, b, d = sp.symbols("a b d", real=True)
    stiffness = sp.Matrix([[a, b], [b, d]])
    for spin in (1, -1):
        e = sp.Matrix([1, sp.I * spin]) / sp.sqrt(2)
        show(f"circular {spin:+} diagonal end stiffness", sp.expand((e.H * stiffness * e)[0]))

    # Sensitivity to an added term in the modal action, supported off the node.
    mass, omega0, gamma, omega = sp.symbols("M omega0 gamma omega", positive=True)
    x, y = sp.Function("x")(t), sp.Function("y")(t)
    for side, spin in product((1, -1), repeat=2):
        lag = (mass * (x.diff(t)**2 + y.diff(t)**2 - omega0**2 * (x**2 + y**2))
               + side * gamma * (x * y.diff(t) - y * x.diff(t))) / 2
        orbit = {x: sp.cos(omega * t), y: spin * sp.sin(omega * t)}
        eqs = [sp.diff(lag.diff(u.diff(t)), t) - lag.diff(u) for u in (x, y)]
        factors = [sp.simplify(eq.subs(orbit).doit() / orbit[u]) for eq, u in zip(eqs, (x, y))]
        assert sp.simplify(factors[0] - factors[1]) == 0
        shifted = (sp.sqrt(gamma**2 + 4 * mass**2 * omega0**2) - side * spin * gamma) / (2 * mass)
        assert sp.simplify(factors[0].subs(omega, shifted)) == 0
        show(f"added gyroscopic term, side {side:+} polarization {spin:+}, frequency", shifted)


def material_contact():
    w = sp.symbols("w", real=True)
    c, k, kappa = sp.symbols("c k kappa", positive=True)
    h = sp.symbols("h", real=True)
    potential, state = sp.Function("V"), sp.Function("f")
    def operator(f, v):
        return -c * sp.diff(f, w, 2) + v * f

    mirror_difference = operator(state(-w), potential(-w)) - operator(state(w), potential(w)).subs(w, -w)
    show("reflected differential operator residual", sp.simplify(mirror_difference))
    refl, trans = sp.symbols("r t")
    for side in (1, -1):
        incoming = sp.exp(sp.I * side * k * w) + refl * sp.exp(-sp.I * side * k * w)
        outgoing = trans * sp.exp(sp.I * side * k * w)
        contact = c * side * (outgoing - incoming).diff(w).subs(w, 0) - h * trans
        sol = sp.solve([(incoming - outgoing).subs(w, 0), contact], (refl, trans), dict=True)[0]
        show(f"side {side:+} transmitted current / incident current", sp.factor(abs(sol[trans])**2))
    jump = c * (sp.diff(sp.exp(-kappa*w), w) - sp.diff(sp.exp(kappa*w), w)).subs(w, 0)
    decay = sp.solve(jump - h, kappa)[0]
    show("bound mode (inverse length, energy) when inverse length is positive", (decay, -c * decay**2))
    real_tail = sp.exp(-kappa * w)
    show("bound mode number current (arbitrary normalization)", sp.im(sp.conjugate(real_tail) * real_tail.diff(w)))


def flux_endpoints():
    sectors = Counter()
    for shared, *outer in product((-1, 1), repeat=7):
        # Twice each outward flux is an integer. The shared link has opposite normals.
        q_left = (shared + sum(outer[:3])) // 2
        q_right = (-shared + sum(outer[3:])) // 2
        if q_left + q_right == 0:
            sectors[q_left] += 1
    charges = sorted(sectors)
    show("neutral pair sectors (left charge, number of configurations)", sorted(sectors.items()))
    for coupling in (0, 1, 2):  # units lambda/U
        for side in (1, -1):
            energies = {q: q**2 - 2 * coupling * side * q for q in charges}
            show(f"lambda/U {coupling} left side {side:+} ground charge",
                 [q for q in charges if energies[q] == min(energies.values())])
    ratio = sp.symbols("lambda_over_U", real=True)
    relative_costs = [q**2 - 2 * ratio * q - (1 - 2 * ratio) for q in charges if q != 1]
    show("strict unit-selection interval", sp.reduce_inequalities([cost > 0 for cost in relative_costs], ratio))


if __name__ == "__main__":
    standing_wave()
    material_contact()
    flux_endpoints()
