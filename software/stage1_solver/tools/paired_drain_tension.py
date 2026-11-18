"""Rate-dependent material traction on a 3D sphere around one held point drain."""

import sympy as sp


def main():
    # b < d. Q is the positive inward volume rate in the 3D brane.
    b, d, t, eta, qs, qt = sp.symbols("b d t eta Q_s Q_t", positive=True)
    axial = sp.symbols("x", real=True)
    a = d/2
    # Between two equal sinks about the midpoint, the axial velocities oppose.
    v_between = qs*((a-axial)**-2 - (a+axial)**-2)/(4*sp.pi)
    stretch = sp.factor(sp.diff(v_between, axial).subs(axial, 0))
    print("equal pair midpoint velocity, axial stretch rate =", v_between.subs(axial, 0), stretch)
    assert stretch == 8*qs/(sp.pi*d**3)

    # t is the distance from the integration point to the companion. The
    # isolated target stress cancels pointwise in the paired-minus-alone flux.
    mu = (d**2+b**2-t**2)/(2*b*d)
    # Contract the harmonic companion Hessian with the sphere normal n and x.
    # Offset from source is b*n-d*ex; its axial and normal components follow.
    sx, sn = b*mu-d, b-d*mu
    tx = eta*qs/(2*sp.pi) * (3*sx*sn/t**5 - mu/t**3)
    area = 2*sp.pi*b*t/d
    primitive = sp.integrate(sp.expand(tx*area), t)
    equator = sp.sqrt(d**2+b**2)
    front = sp.factor(primitive.subs(t, equator) - primitive.subs(t, d-b))
    back = sp.factor(primitive.subs(t, d+b) - primitive.subs(t, equator))
    print("source-facing hemisphere =", front)
    print("opposite hemisphere =", back)
    print("complete sphere =", sp.simplify(front+back))
    assert sp.simplify(front-eta*qs*b**2/(d**2+b**2)**sp.Rational(3, 2)) == 0
    assert sp.simplify(front+back) == 0

    # Sensitivity of the existing P6 viscosity term only. Evaluate the
    # isolated radial strain with a possible dipolar change in its coefficient.
    mu, u, beta = sp.symbols("mu U beta", real=True)
    target_velocity = -qt/(4*sp.pi*b**2)
    radial_stress_change = 2*beta*u*mu*sp.diff(target_velocity, b)
    force_change = 2*sp.pi*b**2*sp.integrate(mu*radial_stress_change, (mu, -1, 1))
    print("viscosity dipole contribution per unit incident speed =", sp.factor(force_change/u))
    assert sp.simplify(force_change - 4*qt*beta*u/(3*b)) == 0


if __name__ == "__main__":
    main()
