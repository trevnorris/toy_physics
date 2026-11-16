"""Interaction momentum crossing a plane between two steady point intakes."""

import sympy as sp


def main():
    x = sp.symbols("x", real=True)
    s, left, right = sp.symbols("s l b", positive=True)
    q1, q2, rho = sp.symbols("Q_1 Q_2 rho", positive=True)

    def plane_velocity(strength, center):
        potential = strength / (4 * sp.pi * sp.sqrt((x - center)**2 + s**2))
        return sp.Matrix([sp.diff(potential, c).subs(x, 0) for c in (x, s)])

    def current_xx(velocity):
        return rho * (velocity[0]**2 - velocity.dot(velocity) / 2)

    v1 = plane_velocity(q1, -left)
    v2 = plane_velocity(q2, right)
    stress = sp.factor(current_xx(v1 + v2) - current_xx(v1) - current_xx(v2))
    ring_current = 2 * sp.pi * s * stress
    # The cosine between the two source-to-ring directions parametrizes the plane.
    r1, r2 = sp.Matrix([left, s]), sp.Matrix([-right, s])
    cosine = r1.dot(r2) / (r1.norm() * r2.norm())
    per_cosine = sp.factor(ring_current / sp.diff(cosine, s))
    if per_cosine.has(s):
        raise ValueError("The flux per unit cosine is not constant on this separator")
    through_plane = sp.factor(per_cosine * (sp.limit(cosine, s, sp.oo) - cosine.subs(s, 0)))
    midplane = sp.integrate(ring_current.subs(right, left), (s, 0, sp.oo))
    arrival1, arrival2 = -through_plane, through_plane

    # k_i is a sample mouth law: signed outward export as a multiple of local arrival.
    # All components here are in the common x direction, from mouth 1 to mouth 2.
    k1, k2 = sp.symbols("k_1 k_2", real=True)
    export = sp.factor(k1 * arrival1 + k2 * arrival2)
    forces = (arrival1 * (1-k1), arrival2 * (1-k2))
    pair_force = sp.factor(sum(forces))
    witness = {k1: 0, k2: 1}
    witness_ratios = [sp.cancel(z.subs(witness) / arrival1) for z in (*forces, export)]

    print("interaction xx momentum current on the separator =", stress, flush=True)
    print("outward x flux from the left across the whole plane =", through_plane, flush=True)
    print("bisector direct integral =", midplane, flush=True)
    print("signed exterior arrivals at mouths 1 and 2 =", (arrival1, arrival2), flush=True)
    print("sum of drives for the sample export law =", pair_force, flush=True)
    print("summed signed export =", export, flush=True)
    print("example (force 1, force 2, net export) / positive arrival =", witness_ratios, flush=True)

    assert sp.simplify(through_plane.subs(right, left) - midplane) == 0
    assert sp.simplify(through_plane + rho*q1*q2 / (4*sp.pi*(left+right)**2)) == 0
    assert sp.simplify(pair_force + export) == 0 and witness_ratios == [1, 0, -1]


if __name__ == "__main__":
    main()
