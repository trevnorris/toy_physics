"""Test the separation-6 companion load against the separation-4 and -8 predictions."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import jn_zeros

from lateral_companion_boundary_term import axial_kernel, radial_kernel
from lateral_companion_load import CompanionTrace
from lateral_mouth_radius import Target, held_reading, steady
from lateral_paired_acceleration import PairedTranslation
from lateral_reservoir_field import LateralSlab
from lateral_reservoir_readout import mouth_load
from lateral_source_rate import verified_packet


MEASUREMENTS = Path("docs/_measurements")
DISTANCE = 6
RADIUS = 1.5


def predictions(extent):
    archive = json.loads((MEASUREMENTS/"lateral_companion_separation.json").read_bytes())
    grid = f"e{extent}_c4"
    pair, = (p for p in archive["calculation"]["pairs"]
             if p["grid"] == grid and p["target_radius"] == RADIUS)
    result = []
    for endpoint in (pair["near"], pair["far"]):
        ratio = endpoint["d"] / DISTANCE
        result.append({"from_d": endpoint["d"], "static": endpoint["static_load"]*ratio**2,
                       "corrected": endpoint["residual"]*ratio**3,
                       "direct_corrected": endpoint["direct_residual"]*ratio**3})
    print("predictions " + json.dumps(result, allow_nan=False), flush=True)
    return pair, result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("extent", type=int, choices=(24, 48))
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    pair, forecasts = predictions(args.extent)
    preparations = json.loads((MEASUREMENTS/"lateral_paired_acceleration.json").read_bytes())
    prep, = (r for r in preparations["runs"] if r["grid"] == pair["grid"] and r["radius"] == RADIUS)
    if prep["guard_exit"] or prep["worker_exit"] or not prep["limits_verified"]:
        raise ValueError("Unverified reference preparation")
    name, = prep["reference_fields"]
    slab = LateralSlab(args.extent, 4, 16)
    with verified_packet(Path(name), {Path(name): prep["reference_fields"][name]}) as data:
        if (int(data["d"]) != 0 or float(data["dz"]) != slab.dz or
                not np.array_equal(slab.x, data["x"]) or not np.array_equal(slab.r, data["r"])):
            raise ValueError("Calibration geometry mismatch")
        ref = abs(data["field"][-1])**2
        mean = np.sum(slab.dv*data["currents"][0]/ref)/(4*np.pi/3)
    for source in (pair["near"], pair["far"]):
        if (source["dx"] != slab.dx or source["dz"] != slab.dz or
                not np.isclose(mean, source["mean_speed"], rtol=1e-11, atol=0)):
            raise ValueError("Calibration differs from the endpoint runs")
    args.output.mkdir(exist_ok=False)
    target = Target(slab, RADIUS, ref, mean)
    base, currents = steady(slab, target, DISTANCE, ref)
    held_reading(slab, target, base, currents, DISTANCE)  # Checks field residual and bath flux.
    motion = PairedTranslation(target, base, currents, DISTANCE, ref)
    trace = CompanionTrace(slab, base[-4:], currents[0], motion.source_g, motion.source_shape)
    static = trace.nonlinear_load(trace.face)
    # mouth_load uses a centered unit mouth as the last current and counts -x as inward.
    static_parts = -mouth_load(base[-4:], currents[::-1], slab.x+DISTANCE, slab.r, slab.dz, 2)
    if not np.isclose(static, sum(static_parts[:2]), rtol=1e-12, atol=1e-15):
        raise ArithmeticError("Static companion readouts disagree on the footprint")
    p, rp = motion.static_tangent()
    bp = motion.boundary_check(p, 1)
    ix = np.flatnonzero(abs(slab.x) < 10)
    crop = slice(ix[0], ix[-1]+1)
    radial = slice(0, np.count_nonzero(slab.r < 4))
    with (args.output/"base.npz").open("xb") as stream:
        np.savez_compressed(stream, x=slab.x, r=slab.r, dz=slab.dz, d=DISTANCE, target_radius=RADIUS,
                            field=base, currents=currents, local_x=slab.x[crop], local_r=slab.r[radial],
                            p=p[:, crop, radial])
    h, rh, steps = motion.dynamic_tangents(p, 0)
    bh = motion.boundary_check(h, 0)
    carried, stress, direct = trace.load(h[-4:])
    observed = carried + stress
    step = 1e-4/max(1, np.max(abs(h[-1])))
    fd = [trace.acceleration_difference(h[-1], eps) for eps in (step, step/2)]
    if (not np.all(np.isfinite([bp, bh, *static_parts, observed, direct, *fd]))
            or max(bp, bh) > 3e-8 or not np.allclose(fd, observed, rtol=2e-6, atol=2e-11)):
        raise ArithmeticError("Translation boundary or companion force derivative check failed")
    with (args.output/"h_0.npz").open("xb") as stream:
        np.savez_compressed(stream, h=h[:, crop, radial], omega=0.)
    companion_rate, target_rate = np.sum(currents*slab.dv, axis=(1, 2))
    fraction = DISTANCE / args.extent
    kernel = axial_kernel(fraction)
    if abs(kernel-radial_kernel(fraction, jn_zeros(0, 512), 10)) > 2e-12:
        raise ArithmeticError("The two boundary kernel calculations disagree")
    boundary = -target_rate * companion_rate * kernel / args.extent  # Unit slab depth.
    corrected = observed - boundary
    direct_corrected = direct - boundary
    readings = {"static": static, "corrected": corrected, "direct_corrected": direct_corrected}
    comparison = [{"from_d": pred["from_d"], **{key: readings[key]/pred[key]-1 for key in readings}}
                  for pred in forecasts]
    result = dict(L=args.extent, d=DISTANCE, target_radius=RADIUS, companion_radius=1, omega=0,
                  dx=slab.dx, dz=slab.dz, mean_speed=mean, target_rate=target_rate,
                  companion_rate=companion_rate, static_direct=sum(static_parts[2:]), **readings,
                  carried=carried, stress=stress, observed=observed, direct=direct, kernel=kernel,
                  boundary=boundary, fd=fd, tangent_residuals=[rp, rh], boundary_errors=[bp, bh],
                  steps=steps, relative_errors=comparison)
    print("measurement " + json.dumps(result, allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
