"""Virtual work of the parent confinement potential for a sideways mouth shift."""

import sympy as sp


def main():
    x, mu, speed = sp.symbols("X mu U", real=True)
    radius, r, ell = sp.symbols("a r ell", positive=True)
    # Intersect a ray from the old origin with the translated sphere.
    rim = x * mu + sp.sqrt(radius**2 - x**2 * (1 - mu**2))
    rim_shift = sp.diff(rim, x).subs(x, 0)
    volume_shift = 2 * sp.pi * sp.integrate(radius**2 * rim_shift, (mu, -1, 1))

    n0, n1, n2, n3 = sp.symbols("n0 n1 n2 n3", real=True)
    density = n0 + speed * (n1 * mu + n2 * sp.legendre(2, mu)) + n3*speed**3*mu
    wall = sp.Function("Vwall")
    moved_potential = wall((r - rim) / ell)
    vr = sp.diff(moved_potential.subs(x, 0), r)
    # Only the geometric coordinate is varied in this partial derivative.
    work = -sp.diff(r**2 * density * moved_potential, x).subs(x, 0)
    reaction = sp.simplify(work - r**2 * density * vr * mu)
    radial_force = sp.simplify(2 * sp.pi * sp.integrate(work, (mu, -1, 1)))
    vr_symbol = sp.Symbol("V_r", real=True)
    per_speed = sp.factor(sp.diff(radial_force, speed).subs(speed, 0).subs(vr, vr_symbol))
    cubic = sp.factor((sp.diff(radial_force, speed, 3)/6).subs(vr, vr_symbol))

    w, width, length = sp.symbols("w h L", positive=True)
    f0, f1 = sp.symbols("f0 f1", real=True)
    # Both trial translations fix the far end and coincide in value and
    # slope at both ends. They differ only inside. Neither is a wall solution.
    s = w / length
    first = (1 - s)**2 * (1 + 2*s)
    second = first + s**2 * (1 - s)**2
    difference = second - first
    traces = [sp.simplify(expr.subs(w, end)) for end in (0, length)
              for expr in (difference, sp.diff(difference, w))]
    extra_work = sp.integrate(f0 * difference, (w, 0, length))

    collar = first.subs(length, width)
    collar_work = sp.integrate((f0 + f1*w) * collar, (w, 0, width))
    contact_weight = collar.subs(w, 0)
    # Elastic energy for the translation eta=q(w)*mu, using densitized T_w.
    slope, stiffness = sp.symbols("q_w T_w", real=True)
    elastic = 2 * sp.pi * sp.integrate(stiffness * (slope * mu)**2 / 2, (mu, -1, 1))
    # At the lower end, minus the variation of this energy gives +dE/dq_w.
    near_force = sp.diff(elastic, slope)
    print("displacement projection, volume change =", rim_shift, volume_shift)
    print("force minus opposite material body force (radial measure included) =", reaction)
    print("undriven radial force integrand =", radial_force.subs(speed, 0))
    print("force per unit w and U, before radial integration =", per_speed)
    print("cubic force coefficient at fixed cross-sections =", cubic)
    print("trial extension trace differences (value,slope at 0,L) =", traces)
    print("difference of distributed works per unit mouth displacement =", extra_work)
    print("regular collar work per unit mouth displacement =", sp.factor(collar_work))
    print("shrinking collar limit, endpoint force weight =", sp.limit(collar_work, width, 0), contact_weight)
    print("axial elastic force on the lower endpoint in physical dipole normalization =", near_force)
    if (rim_shift != mu or volume_shift != 0 or reaction != 0 or radial_force.subs(speed, 0) != 0
            or sp.simplify(per_speed - 4*sp.pi*r**2*n1*vr_symbol/3) != 0
            or sp.simplify(cubic - 4*sp.pi*r**2*n3*vr_symbol/3) != 0
            or traces != [0]*4 or sp.simplify(extra_work - f0*length/30) != 0
            or sp.limit(collar_work, width, 0) != 0 or contact_weight != 1
            or sp.simplify(near_force - 4*sp.pi*stiffness*slope/3) != 0):
        raise ArithmeticError("Inconsistent angular projection or endpoint work")


if __name__ == "__main__":
    main()
