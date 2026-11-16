"""Steady held pair. All fluxes are projected toward the companion.

Integrate the full momentum current on a finite control sphere in the
reduced 3D two-monopole exterior, and the normal stress on a native graph
face with arbitrary radial profiles. The control sphere is outside the core.
"""
import sympy as s


def exterior():
    b, d, q1, q2, rho, t = s.symbols("b d Q1 Q2 rho t", positive=True)
    # Target at zero, companion at +d. Integrate by distance t to companion,
    # which ranges from d-b to d+b. The signed velocity is the gradient of
    # the positive intake potential. Hold the target's own flux fixed.
    mu = (d**2 + b**2 - t**2)/(2*d*b)
    own = -q1/(4*s.pi*b**2)
    comp = q2/(4*s.pi*t**3)
    ux = own*mu + comp*(d-b*mu)
    un = own + comp*(d*mu-b)
    speed_sq = own**2 + comp**2*t**2 + 2*own*comp*(d*mu-b)
    jac = 2*s.pi*b*t/d

    def inward(current):
        primitive = s.integrate(s.expand(-jac*current), t)
        return s.factor(primitive.subs(t, d+b) - primitive.subs(t, d-b))

    convective = inward(rho*ux*un)
    pressure = inward(-rho*speed_sq*mu/2)
    total = s.factor(convective + pressure)
    isolated = total.subs(q2, 0)
    change = s.factor(total - isolated)
    print("target only, companion only, pair on target contour =", isolated,
          total.subs(q1, 0), total)
    print("bilinear inward convective, pressure =",
          s.factor(s.diff(convective, q1, q2)*q1*q2),
          s.factor(s.diff(pressure, q1, q2)*q1*q2))
    print("pair minus isolated =", change)
    print("sphere-radius derivative =", s.factor(s.diff(change, b)))
    for distance, radius in [(10, 1), (10, 3), (20, 1), (20, 3)]:
        scaled = 4*s.pi*change/(rho*q1*q2)
        print(f"scaled d={distance} b={radius} =", scaled.subs({d: distance, b: radius}))


def upper_face():
    # Profiles and signed outward stress amplitude on the upper graph. These
    # are independent perturbation inputs to the balance until a face law
    # relates them. Include variation of both the traction and the face.
    r = s.Symbol("r", positive=True)
    mu, eps = s.symbols("mu eps", real=True)
    h0, h1, p0, p1 = (s.Function(n)(r) for n in ("h0", "h1", "p0", "p1"))
    h = h0 + eps*h1*mu
    p = p0 + eps*p1*mu
    dx_h = mu*s.diff(h, r) + (1-mu**2)*s.diff(h, mu)/r
    angular = 2*s.pi*s.integrate(-r**2*p*dx_h, (mu, -1, 1))
    print("upper face isolated radial integrand =", angular.subs(eps, 0))
    print("upper face first variation per dr =", s.factor(s.diff(angular, eps).subs(eps, 0)))


if __name__ == "__main__":
    exterior()
    upper_face()
