"""Two-sided cylindrical bending of a resting ordered strip, with an optional material strain law."""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_bvp
from scipy.special import expit

from mouth_curved_order_patch import well

DELTA = .04
H = 1.
WELL = 75.
SIGMA = WELL * DELTA / 3.
CALIBRATION = .08
HELD_OUT = .04


def solve_strip(curvature, modulus, extension=0., *, resolution=7, extent=18, profile=False):
    """Hold mid-area (except in flat extension controls), half-order heights and bath pressures."""
    faces = np.array([-H / 2, H / 2]) / DELTA
    edges = [faces[0] - extent, *faces, faces[1] + extent]
    pins = [0., .5, .5, 0.]
    abscissas, weights = leggauss(8)
    readings = dict(energy=0., order_moment=0., crystal_moment=0., order_tension=0.,
                    crystal_tension=0., registered_second_moment=0., registered_thickness=0.)
    samples = []
    max_residual = 0.
    for idx, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        t = np.linspace(lo, hi, int(np.ceil((hi - lo) * resolution)) + 1)
        left, right = expit(t - faces[0]), expit(faces[1] - t)
        guess = np.array([left * right, left * right * (right - left)])

        def ode(t, field):
            chi, slope = field
            strain = extension + curvature * DELTA * t
            area = 1. + strain
            if np.any(area <= 0.):
                raise ValueError("Parallel surface has crossed a focal point")
            order_force = chi * (1. - chi) * (1. - 2. * chi)
            elastic_force = modulus * chi * strain**2 / (2. * WELL * area)
            return np.stack((slope, order_force + elastic_force - curvature * DELTA * slope / area))

        def boundaries(a, b):
            return np.array([a[0] - pins[idx], b[0] - pins[idx + 1]])

        sol = solve_bvp(ode, boundaries, t, guess, tol=3e-9, max_nodes=20000)
        if not sol.success:
            raise ArithmeticError(sol.message)
        max_residual = max(max_residual, max(sol.rms_residuals))

        def components(t):
            chi, slope = sol.sol(t)
            if np.min(chi) < -1.e-8 or np.max(chi) > 1. + 1.e-8:
                raise ValueError("Order left its declared interval")
            z = DELTA * t
            strain = extension + curvature * z
            order_stress = WELL * (well(chi) + slope**2)
            crystal_stress = modulus * chi**2 * strain
            return z, chi, strain, order_stress, crystal_stress

        half = np.diff(sol.x)[:, None] / 2.
        points = (sol.x[:-1, None] + half * (1. + abscissas)).ravel()
        measure = (DELTA * half * weights).ravel()
        z, chi, eps, base, added = components(points)
        readings["energy"] += measure @ ((1. + eps) * base + .5 * eps * added)
        readings["order_moment"] += measure @ (z * base)
        readings["crystal_moment"] += measure @ (z * added)
        readings["order_tension"] += measure @ base
        readings["crystal_tension"] += measure @ added
        readings["registered_second_moment"] += measure @ (z**2 * chi**2)
        readings["registered_thickness"] += measure @ (chi**2)
        if lo <= 0. <= hi:
            readings["center_order"] = float(sol.sol(0.)[0])
        if profile:
            # Retain both limiting readings at a registered (pinned) half-order face.
            zs = np.linspace(lo, hi, int(round(hi - lo)) * 4 + 1)
            z, chi, eps, base, added = components(zs)
            samples.append(dict(region=idx, depth=z.tolist(), fraction=chi.tolist(),
                                order_axial_stress=base.tolist(), crystal_axial_stress=added.tolist()))
    readings.update(curvature=curvature, extension=extension, modulus=modulus,
                    holding_moment=readings["order_moment"] + readings["crystal_moment"],
                    end_tension=readings["order_tension"] + readings["crystal_tension"],
                    max_residual=max_residual)
    if profile:
        readings["stress_profiles"] = samples
    return readings


def experiment(modulus, mesh, extent, include_profiles):
    def solve(k, stretch=0., profiles=False):
        return solve_strip(k, modulus, stretch, resolution=mesh, extent=extent, profile=profiles)

    flat = solve(0., profiles=include_profiles)
    train = [solve(k, profiles=include_profiles) for k in (-CALIBRATION, CALIBRATION)]
    stiffness = (train[1]["holding_moment"] - train[0]["holding_moment"]) / (2. * CALIBRATION)
    intercept = (train[0]["holding_moment"] + train[1]["holding_moment"]) / 2.
    predictions = [intercept + stiffness * k for k in (-HELD_OUT, HELD_OUT)]
    check = [solve(k, profiles=include_profiles) for k in (-HELD_OUT, HELD_OUT)]
    step = 1.e-4
    work = []
    for point in train + check:
        k = point["curvature"]
        lo, hi = solve(k-step), solve(k+step)
        derivative = (hi["energy"] - lo["energy"]) / (2. * step)
        work.append(dict(curvature=k, energy_lo=lo["energy"], energy_hi=hi["energy"],
                         derivative=derivative, stress_moment=point["holding_moment"],
                         difference=derivative-point["holding_moment"]))
    # The small-curvature tangent and its field distribution are independent of the finite-pair fit.
    slopes = []
    for k in (.01, .005):
        minus, plus = solve(-k), solve(k)
        slopes.append((plus["holding_moment"] - minus["holding_moment"]) / (2. * k))
    tangent = (4. * slopes[1] - slopes[0]) / 3.
    # A matched flat-extension experiment gives the separate macroscopic axial stiffness.
    e = .001
    flat_pair = [solve(0., -e), solve(0., e)]
    extension_stiffness = (flat_pair[1]["end_tension"] - flat_pair[0]["end_tension"]) / (2. * e)
    return dict(flat=flat, training=train, fit_stiffness=stiffness, fit_intercept=intercept,
                predictions=predictions, held_out=check, moment_work=work, work_step=step,
                small_pair_stiffness=slopes, tangent=tangent,
                flat_pair=flat_pair, flat_incremental_axial_stiffness=extension_stiffness)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolution", required=True, type=int)
    parser.add_argument("--extent", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--profiles", action="store_true")
    args = parser.parse_args()
    native, registered = [experiment(y, args.resolution, args.extent, args.profiles)
                          for y in (0., SIGMA / H)]
    i2 = native["flat"]["registered_second_moment"]
    yc = -native["tangent"] / i2
    comparison = dict(flat_I2=i2, measured_added_tangent=registered["tangent"]-native["tangent"],
                      predicted_added_tangent=SIGMA/H*i2,
                      critical_YH_over_sigma=yc*H/SIGMA)
    result = dict(thickness=H, delta=DELTA, well_scale=WELL, tension_one_face=SIGMA,
                  curvature_train=CALIBRATION, curvature_test=HELD_OUT,
                  resolution=args.resolution, exterior_extent_in_delta=args.extent,
                  native=native, registered_trial=registered, comparison=comparison)
    args.output.write_text(json.dumps(result, allow_nan=False, separators=(",", ":"))+"\n")
    for name, data in (("active", native), ("registered_trial", registered)):
        for k, pred, actual in zip((-HELD_OUT, HELD_OUT), data["predictions"], data["held_out"]):
            value = actual["holding_moment"]
            print(json.dumps(dict(law=name, curvature=k, prediction=pred, direct=value,
                                  relative_error=abs(pred/value-1.), calibration=data["fit_stiffness"],
                                  tangent=data["tangent"], max_work_error=max(
                                      abs(row["difference"]) for row in data["moment_work"]),
                                  flat_incremental_axial_stiffness=data["flat_incremental_axial_stiffness"])))
    print(json.dumps(comparison), flush=True)
    if any(abs(row["difference"]) > 2.e-8
           for data in (native, registered) for row in data["moment_work"]):
        raise ArithmeticError("Through-thickness stress and independently differenced energy disagree")
    if abs(comparison["measured_added_tangent"] - comparison["predicted_added_tangent"]) > 2.e-8:
        raise ArithmeticError("Flat strip and curvature tangent disagree")


if __name__ == "__main__":
    main()
