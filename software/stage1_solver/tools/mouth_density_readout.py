"""Check sampling of manufactured number-density fields on a one-throat aperture.

At a fixed radius in the isolated aperture, mu is the cosine to the incident
velocity axis and s is its signed speed. The entrance samples a gradient in w.
The manufactured fields test the readouts only. They obey no assumed core law.
"""
import sympy as sp


def main():
    s, mu, w = sp.symbols("s mu w", real=True)
    n0, grad, tilt, h, even, cubic = sp.symbols("n0 g eta h e k", real=True)
    density = n0 + grad*w + s*h*mu + even*s**2*mu**2 + cubic*s**3*mu**3
    entrance = s*tilt*mu
    sampled = density.subs(w, entrance)
    difference = sampled - sampled.subs(s, 0)
    weight = sp.integrate(mu**2, (mu, -1, 1))
    projection = sp.cancel(sp.integrate(mu*difference, (mu, -1, 1))/(s*weight))
    axial = sp.cancel((difference.subs(mu, 1) - difference.subs(mu, -1))/(2*s))
    print("sampled density difference =", sp.expand(difference))
    print("angular [numerator, normalization] =",
          [sp.integrate(mu*difference, (mu, -1, 1)), s*weight])
    print("[projected, opposing-point] slopes =", [projection, axial])
    print("[projected, opposing-point] zero-drive limits =",
          [sp.limit(reading, s, 0) for reading in (projection, axial)])
    print("half-speed and full-speed projected slopes =", [projection.subs(s, s/2), projection])
    true_response = sp.diff(sampled, s).subs(s, 0)
    fixed_response = sp.diff(density, s).subs({w: 0, s: 0})
    print("linear density [fixed location, moving entrance, difference] =",
          [fixed_response, true_response, sp.simplify(true_response - fixed_response)])

    # Local inward relative speed and signed current at corresponding points.
    speed0 = sp.Symbol("v0", positive=True)
    dv = sp.Symbol("v1", real=True)
    inward_speed = speed0 + s*dv*mu
    current = sampled*inward_speed
    dj = sp.diff(current, s).subs(s, 0)
    inferred = sp.simplify((dj - n0*sp.diff(inward_speed, s))/speed0)
    print("current [baseline, derivative] =", [current.subs(s, 0), sp.expand(dj)])
    print("density derivative from current and relative speed =", inferred)
    if any(sp.simplify(sp.limit(reading, s, 0) - true_response.subs(mu, 1)) != 0
           for reading in (projection, axial)) or sp.simplify(inferred-true_response) != 0:
        raise ValueError("The measured linear slopes disagree.")


if __name__ == "__main__":
    main()
