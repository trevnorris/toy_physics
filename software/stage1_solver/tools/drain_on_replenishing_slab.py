"""Finite spherical shell in a prescribed, flat slab. Solve constituent continuity.

The bulk supplies the two faces uniformly. The free input q is the additional
number per unit time entering the core relative to the smooth reference field.
All radial velocities and fluxes are signed outward unless labelled inward.
"""
import sympy as sp


r, a, h, t = sp.symbols("r a h t", positive=True)
w = sp.Symbol("w", real=True)
n0, b, ell, q = sp.symbols("n0 b ell q", positive=True)
n = n0 + b * t
jw = -ell * w
radial_current = sp.Function("j")

# Integrate the local continuity equation. The regular background extends
# through the reference sphere at a. The throat changes the flux at that cut.
radial_divergence = -sp.diff(n, t) - sp.diff(jw, w)
regular_j = sp.integrate(r**2 * radial_divergence, r) / r**2
core_flux = -4 * sp.pi * a**2 * h * regular_j.subs(r, a) + q
ode = sp.Eq(sp.diff(r**2 * radial_current(r), r), r**2 * radial_divergence)
throat_j = sp.dsolve(ode, ics={radial_current(a): -core_flux / (4 * sp.pi * a**2 * h)}).rhs


def shell(j):
    area = 4 * sp.pi * r**2
    projected_flux = sp.integrate(j, (w, -h/2, h/2))
    phi = area * projected_flux
    core = -phi.subs(r, a)
    faces = sp.integrate(area * (jw.subs(w, h/2) - jw.subs(w, -h/2)), (r, a, r))
    storage = sp.diff(sp.integrate(area * h * n, (r, a, r)), t)
    return (projected_flux / (h * n), phi, core, faces, storage)


def show(label, values):
    print(label, [sp.factor(value) for value in values], flush=True)


background = shell(regular_j)
throat = shell(throat_j)
excess_v = sp.factor(throat[0] - background[0])
background_budget = sum(background[2:])
throat_budget = sum(throat[2:])
show("background [outward radial speed, outward flux, inward core, outward faces, storage]", background)
show("throat [outward radial speed, outward flux, inward core, outward faces, storage]", throat)
show("direct extra inward number flux, extra inward from budgets, extra outward speed",
     [background[1] - throat[1], throat_budget - background_budget, excess_v])
show("no-extra-intake speed difference", [excess_v.subs(q, 0)])
show("raw inward flux, budget if storage omitted, missing term",
     [-throat[1], throat[2] + throat[3], -throat[1] - throat[2] - throat[3]])
params = {a: 1, h: 2, n0: 3, b: sp.Rational(3, 100), ell: sp.Rational(3, 25),
          q: 24 * sp.pi, t: 0}
for radius in (1, 5, 10):
    sample = [background[0], throat[0], excess_v, (background[1] - throat[1]) / sp.pi]
    show(f"R={radius} [background speed, throat speed, difference, extra inward flux/pi]",
         [value.subs(params).subs(r, radius) for value in sample])
for state in (background, throat):
    if sp.simplify(sum(state[1:])) != 0:
        raise AssertionError("Finite-shell material budget did not close")
