"""Solve the spherical G0 annulus for a unit generalized source at two radii."""

import sympy as sp


def one_mouth(a, width, stiffness, contact):
    outer = a + width
    eta = 3 / (4 * sp.pi * (outer**3 - a**3))
    beta = sp.sqrt(contact * eta / stiffness)
    # f = 1 - contact*h. In the annulus t = r*f has t'' = beta^2*t.
    # f is regular and constant inside; outside t = r - tail.
    r = sp.symbols("r", positive=True)
    inner_f = 1 / (sp.cosh(beta * width) + a * beta * sp.sinh(beta * width))
    t = inner_f * (a * sp.cosh(beta * (r-a)) + sp.sinh(beta * (r-a)) / beta)
    residuals = (sp.diff(t, r, 2) - beta**2 * t, t.subs(r, a) - a * inner_f,
                 sp.diff(t, r).subs(r, a) - inner_f, sp.diff(t, r).subs(r, outer) - 1)
    if any(sp.simplify(term) != 0 for term in residuals):
        raise AssertionError("radial equation or matching failed")
    tail = outer - t.subs(r, outer)
    escape = 4 * sp.pi * stiffness * tail / contact
    # Compute the average of h separately from the far tail. An elementary
    # primitive of r*t follows from t'' = beta^2*t.
    primitive = (r * sp.diff(t, r) - t) / beta**2
    avg_f = 4 * sp.pi * eta * (primitive.subs(r, outer) - primitive.subs(r, a))
    if sp.simplify(avg_f - escape) != 0:
        raise AssertionError("far conormal disagrees with integrated shell load")
    return escape.evalf(32), ((1 - avg_f)/contact).evalf(32)


responses = [one_mouth(a, sp.Rational(1, 20), sp.Rational(7, 8), 1) for a in (1, 2)]
for a, (z, response) in enumerate(responses, 1):
    print(f"radius {a}, escape = {z:.12f}, self response = {response:.12f}")
for label, strengths, sign in [
    ("held unit mean bend", [z / response for z, response in responses], 1),
    ("unit integrated conormal, replacement", [z for z, _ in responses], 1),
    ("only the unit committed source", [z for z, _ in responses], -1),
]:
    pairs = [sign * strengths[i] * strengths[j] for i, j in ((0, 0), (0, 1), (1, 1))]
    print(f"{label}, 4 pi kappa d^2 Fodd -> " + ", ".join(f"{b:.12f}" for b in pairs))
z1, self1 = responses[0]
z2, self2 = responses[1]
print(f"required bend(2)/bend(1) for equal positive strength = {(self2/z2)/(self1/z1):.12f}")
print(f"required conormal(2)/conormal(1) in replacement = {z1/z2:.12f}")
