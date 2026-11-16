"""Slow companion force, to bilinear order in the two intake rates, in solver units."""

import json
from pathlib import Path

import numpy as np
from scipy.special import eval_jacobi, roots_legendre, spherical_jn

from lateral_companion_boundary_term import axial_kernel
from lateral_mouth_radius import Target
from lateral_reservoir_field import LateralSlab
from lateral_source_rate import verified_packet
from mouth_feed_weak_check import dx4


ARCHIVE = Path("docs/_measurements")
TARGET_RADIUS = 1.5


def surface_moment(a, order, cutoff, radial_nodes):
    """Isolated l=1 Robin moment for a unit current on the unit-depth slab."""
    x, weights = roots_legendre(radial_nodes)
    r = (x + 1)*a/2
    radial_weight = weights*a/2 * r**2
    t = (r/a)**2
    g = (1-t)**3
    normalization = 4*np.pi*np.dot(radial_weight, g)
    mean_square = 4*np.pi*np.dot(radial_weight, r*r*g)/normalization
    if not np.isclose(mean_square, 3*a*a/11, rtol=1e-12):
        raise ArithmeticError("Intake second moment fails its beta-integral check")

    # Test the surface equation w+B_2(g w)=w0 with g times each dipole basis
    # function. All radial integrals omit the common angular factor 4π/3.
    basis = np.array([r/a*eval_jacobi(n, 3, 1.5, 2*t-1) for n in range(order)])
    mass = (basis*(radial_weight*g)) @ basis.T
    nodes, legendre_weights = roots_legendre(8)
    k = (np.arange(cutoff)[:, None]+(nodes+1)/2).ravel()
    dk = np.tile(legendre_weights/2, cutoff)
    transform = spherical_jn(1, np.outer(k, r))
    modes = transform @ ((radial_weight*g)*basis).T
    source = transform @ (radial_weight * 6*r/a**2 * (1-t)**2/normalization)
    b0 = 1/(k*np.tanh(k))
    energy = np.sqrt(k*k+4)
    b2 = 1/(energy*np.tanh(energy))
    measure = (2/np.pi)*dk*k*k
    matrix = mass + (modes.T * (measure*b2)) @ modes
    rhs = modes.T @ (measure * (b0-b2)*source/4)
    coeff = np.linalg.solve(matrix, rhs)
    return (4*np.pi/3)*np.dot(radial_weight*r*g, coeff @ basis)


def weak_grid_response(slab, target, reference, distance):
    """Per product of rates, retaining finite distance, both mouths, depth and grid."""
    rate = np.sum(slab.dv*target.speed)
    empty = np.zeros(slab.eigenvalues.shape)
    p = slab.elliptic(empty.copy(), target.dspeed/rate, 0)
    g = np.maximum(1-(slab.x[:, None]+distance)**2-slab.r**2, 0)**3
    weight = g/np.roll(reference, -round(distance/slab.dx), axis=0)
    weight /= np.sum(slab.dv*weight)
    restoring = g+target.g
    face = np.zeros_like(g)
    for iteration in range(1, 100):
        u = slab.elliptic(-p, -restoring*face, 2)
        defect = np.max(abs(u[-1]-face))
        if defect < 1e-14:
            break
        face = u[-1]
    else:
        raise ArithmeticError("Real-amplitude boundary equation did not converge")
    h = slab.elliptic(u.copy(), np.zeros_like(g), 0)
    response = np.sum(slab.dv*weight*dx4(h[-1], slab.dx))
    # On these compact footprints the antisymmetry of the axial stencil
    # gives an independent integrated-by-parts force readout.
    by_parts = -np.sum(slab.dv*dx4(weight, slab.dx)*h[-1])
    if not np.isclose(response, by_parts, rtol=1e-11, atol=1e-16):
        raise ArithmeticError("Discrete receiver projection mismatch")
    return response, float(defect), iteration


def main():
    settings = [(6, 70, 120), (10, 130, 210), (14, 200, 270)]
    moments = [[surface_moment(a, *s) for a in (1, TARGET_RADIUS)] for s in settings]
    if not np.allclose(moments[-2], moments[-1], rtol=1e-6, atol=1e-10):
        raise ArithmeticError(f"Isolated-mouth quadrature disagreement: {moments}")
    # Radii, slab depth, sound speed and healing scale in the solver's units.
    terms = dict(density_gradient=0.25, aperture=(1+TARGET_RADIUS**2)/22,
                 companion_robin=moments[-1][0], target_robin=moments[-1][1])
    coefficient = -sum(terms.values())/(2*np.pi)
    print("derivation " + json.dumps(dict(terms=terms, per_product=coefficient,
          settings=settings, robin_moments=moments), allow_nan=False), flush=True)

    endpoints = json.loads((ARCHIVE/"lateral_companion_separation.json").read_bytes())
    midpoint = json.loads((ARCHIVE/"lateral_companion_midpoint.json").read_bytes())
    for run in midpoint["runs"]:
        extent = run["L"]
        mid = run["measurement"]
        pair, = (p for p in endpoints["calculation"]["pairs"]
                 if p["grid"] == f"e{extent}_c4" and p["target_radius"] == TARGET_RADIUS)
        if (run["guard_exit"] or run["worker_exit"] or not run["limits_verified"]
                or endpoints["guard_exit"] or endpoints["worker_exit"] or not endpoints["limits_verified"]):
            raise ValueError("Invalid comparison run")
        name, = run["input_fields"]
        slab = LateralSlab(extent, 4, 16)
        with verified_packet(Path(name), {Path(name): run["input_fields"][name]}) as packet:
            if (int(packet["d"]) or float(packet["dz"]) != slab.dz or
                    not np.array_equal(slab.x, packet["x"]) or not np.array_equal(slab.r, packet["r"])):
                raise ValueError("Incorrect steady calibration geometry")
            ref = abs(packet["field"][-1])**2
            mean = np.sum(slab.dv*packet["currents"][0]/ref)/(4*np.pi/3)
        target = Target(slab, TARGET_RADIUS, ref, mean)
        for entry in (pair["near"], mid, pair["far"]):
            d = entry["d"]
            if not np.isclose(mean, entry["mean_speed"], rtol=1e-11, atol=0):
                raise ValueError("Intake preparation mismatch")
            product = entry["companion_rate"]*entry["target_rate"]
            prepared = 0.25*target.volume*mean
            result, defect, steps = weak_grid_response(slab, target, ref, d)
            boundary = -product*axial_kernel(d/extent)/extent
            if not np.isclose(boundary, entry["boundary"], rtol=1e-11, atol=0):
                raise ArithmeticError("Reservoir subtraction differs from the archived comparison")
            measured = entry["observed"]-boundary
            direct = entry["direct"]-boundary
            predicted = product*coefficient/d**3
            grid_prediction = product*result-boundary
            row = dict(L=extent, d=d, target_rate=entry["target_rate"], companion_rate=entry["companion_rate"],
                       product=product, prepared_product=prepared, coefficient=coefficient*product,
                       prepared_coefficient=coefficient*prepared, measured_d3=d**3*measured,
                       prediction=predicted, measured=measured, direct=direct, weak_grid=grid_prediction,
                       asymptotic_relative_error=predicted/measured-1,
                       weak_grid_relative_error=grid_prediction/measured-1,
                       weak_grid_raw=product*result, boundary=boundary, max_update=defect, steps=steps)
            print("comparison " + json.dumps(row, allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
