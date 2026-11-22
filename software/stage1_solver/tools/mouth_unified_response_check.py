"""Independent four-dimensional elastic comparison and relaxed two-face strip tangent."""
import json

import numpy as np
from scipy.linalg import solve_banded

from mouth_ordered_strip_checks.signed_moments import Strip


def material_test(bulk, shear, released):
    """Unit axial strain in an isotropic material with specified free transverse axes."""
    volumetric = np.full((4, 4), .25)
    stiffness = 4. * bulk * volumetric + 2. * shear * (np.eye(4) - volumetric)
    strain = np.array([1., 0., 0., 0.])
    free = np.arange(1, 1 + released)
    if released:
        strain[free] = np.linalg.solve(stiffness[np.ix_(free, free)], -stiffness[free, 0])
    stress = stiffness @ strain
    return dict(free_laterals=released, strain=strain.tolist(),
                stress=stress.tolist(), axial_modulus=float(stress[0]))


def relaxed(strip, curvature, axial):
    """Minimize the scalar energy plus axial * chi * strain^2/2 with 0 <= chi <= 1."""
    c = strip.seed.copy()
    strain = curvature * strip.s
    force = .5 * axial * strain**2 * strip.quad
    elastic_gradient = np.zeros_like(c)
    elastic_gradient[:-1] += np.sum(strip.left * force, axis=1)
    elastic_gradient[1:] += np.sum(strip.right * force, axis=1)
    for _ in range(80):
        energy, moment, grad, bands = strip.values(c, curvature, 0.)
        grad += elastic_gradient
        active = ~strip.free | ((c == 0.) & (grad >= 0.)) | ((c == 1.) & (grad <= 0.))
        grad[active] = 0.
        if np.max(np.abs(grad)) < 2.e-12:
            mid = c[:-1, None] * strip.left + c[1:, None] * strip.right
            weighted_inertia = np.sum(strip.quad * strip.s**2 * mid)
            return (float(energy + np.sum(mid * force)),
                    float(moment + axial * curvature * weighted_inertia),
                    float(weighted_inertia), c, float(max(abs(grad))))
        pinned_link = active[1:] | active[:-1]
        bands[0, 1:][pinned_link] = 0.
        bands[2, :-1][pinned_link] = 0.
        bands[1, active] = 1.
        c = np.clip(c - solve_banded((1, 1), bands, grad), 0., 1.)
    raise ArithmeticError("Bounded strip did not equilibrate")


def strip_test(n):
    strip = Strip(n, 18.)
    _, _, i2, flat, _ = relaxed(strip, 0., 0.)
    flat_order = flat[:-1, None] * strip.left + flat[1:, None] * strip.right
    i0 = float(np.sum(strip.quad * flat_order))
    _, _, _, bands = strip.values(flat, 0., 0.)
    # The native energy is affine in curvature for fixed order. Schur reduction
    # gives its full relaxed tangent from the mixed force on the flat solution.
    h = .01
    gplus = strip.values(flat, h, 0.)[2]
    gminus = strip.values(flat, -h, 0.)[2]
    mixed = (gplus - gminus) / (2.*h)
    neutral = float(-mixed @ solve_banded((1, 1), bands, mixed))
    rows = []
    for bulk in (.001, .25, 1., 2.):
        for released in (0, 1, 3):
            bench = material_test(bulk, 1., released)
            axial = bench["axial_modulus"]
            tangent = neutral + axial * i2
            pairs = []
            for k in (.02, .01, .0005, .00025):
                a, b = [relaxed(strip, sgn*k, axial) for sgn in (-1., 1.)]
                pairs.append(dict(curvature=k, energy_minus=a[0], energy_plus=b[0],
                                  moment_minus=a[1], moment_plus=b[1],
                                  slope=(b[1]-a[1])/(2.*k),
                                  max_stationarity=max(a[4], b[4])))
            # The finite-curvature source is constant in chi. The bounded
            # outer order tails can switch contact sets over the wider pair.
            direct = (4.*pairs[-1]["slope"]-pairs[-2]["slope"])/3.
            rows.append(dict(bulk=bulk, **bench, tangent=tangent, finite_pairs=pairs,
                             extrapolated_signed_tangent=direct, difference=direct-tangent))
    return dict(elements_per_width=n, neutral_tangent=neutral, order_weighted_I0=i0, order_weighted_I2=i2,
                critical_axial_modulus=-neutral/i2, material_rows=rows)


def main():
    grids = [strip_test(n) for n in (10, 20, 40)]
    fine, medium = grids[-1], grids[-2]
    estimates = {key: (4.*fine[key] - medium[key])/3.
                 for key in ("neutral_tangent", "order_weighted_I0", "order_weighted_I2", "critical_axial_modulus")}
    output = dict(description="Unit ordered shear modulus and transverse inertia in four spatial dimensions. "
                              "The lattice bulk modulus and shear modulus both carry one factor of local "
                              "order. These are wave-off tangent elastic laws for registered material.",
                  delta=.04, thickness=1., bath_extent_in_delta=18.,
                  grids=grids, extrapolation=estimates)
    print(json.dumps(output, allow_nan=False), flush=True)
    if max(abs(r["difference"]) for grid in grids for r in grid["material_rows"]) > 2.e-8:
        raise ArithmeticError("Finite signed strip tangent differs from flat Hessian reduction")


if __name__ == "__main__":
    main()
