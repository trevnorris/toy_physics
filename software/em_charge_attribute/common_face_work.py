"""Virtual normal translation of the two slab faces, wall ends and mouth source."""

import sympy as sp


def check(label, expr):
    value = sp.simplify(expr)
    print(label, "=", value, flush=True)
    if value != 0:
        raise ArithmeticError(label)


def endpoints():
    # z, y are *laboratory* displacements of the mouth and far attachment.
    z, y, beta, orient, load = sp.symbols("z y beta s P", real=True)
    wall = orient*load*(y-z)
    force_along = -sp.diff(wall.subs(y, beta*z), z)
    print("wall work per unit near displacement =", force_along)
    check("rigid wall work", force_along.subs(beta, 1))
    print("fixed far end wall work =", force_along.subs(beta, 0))

    # Identifiability control only. This extra elastic support and bulk force
    # are diagnostic inputs. An unknown prestress can preserve the same
    # initial geometry even when the bulk is exerting a different force.
    k = sp.symbols("k", positive=True)
    bulk, rest = sp.symbols("B d", real=True)
    potential = wall + k*(y-z-rest)**2/2 - bulk*y
    equilibrated = sp.solve(potential.diff(y), y)[0]
    force = sp.factor(-potential.diff(z).subs(y, equilibrated))
    neutral_rest = sp.solve(equilibrated-z, rest)[0]
    print("far end following factor =", equilibrated.diff(z))
    print("slab force after far force balance =", force)
    print("support reference for same near/far displacement =", neutral_rest)
    check("same endpoint positions on either loaded preparation",
          equilibrated.subs(rest, neutral_rest)-z)
    check("total load vs eliminated energy",
          force+sp.diff(potential.subs(y, equilibrated), z))


def cap_projection():
    # Graph of a cap over its 3-ball. For an axisymmetric graph the axial
    # projected area is independent of slope. Positive p acts lumen to bulk.
    r, a, slope = sp.symbols("r a slope", positive=True)
    sigma, bending = sp.symbols("sigma k_b", positive=True)
    da = 4*sp.pi*r**2*sp.sqrt(1+slope**2)
    area = sp.integrate(da/sp.sqrt(1+slope**2), (r, 0, a))
    print("cap projected area =", area)
    axial_pull = 4*sp.pi*sigma*a**2 + 16*sp.pi*bending/9
    required = sp.expand(axial_pull/area)
    print("pressure-only cap support for straight wall =", required)
    check("support reconstructs wall pull", area*required-axial_pull)
    print("pressure doubling-radius ratios [area term, bending term] =",
          [sp.factor(t.subs(a,2*a)/t) for t in
           (required.subs(bending,0),required.subs(sigma,0))])


def mouth_source():
    w, z, ell, orient = sp.symbols("w z ell s", real=True, nonzero=True)
    t = w-z
    kernel, excess = sp.Function("o"), sp.Function("D")
    product = kernel(t)*excess(t)
    check("translated orientation integrand: shift derivative plus spatial derivative",
          sp.diff(product,z)+sp.diff(product,w))
    h, delta, overlap, k, g = sp.symbols("h delta c k_m g", real=True)
    scalar = h + overlap*delta/ell
    energy = k*scalar**2/2-g*orient*scalar
    force = -energy.diff(delta).subs(delta,0)
    print("mouth force for uniform embedding overlap =", sp.factor(force))
    check("transported scalar value has zero mouth work", force.subs(overlap,0))


if __name__ == "__main__":
    endpoints()
    cap_projection()
    mouth_source()
