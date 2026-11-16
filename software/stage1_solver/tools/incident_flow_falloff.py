"""Measure the one-drain incident flow on the fixed footprint of an absent target."""

import argparse
import json
from pathlib import Path

import numpy as np

from lateral_reservoir_field import LateralSlab
from lateral_source_rate import verified_packet
from lateral_target_rate import checked_solve
from mouth_feed_weak_check import dx4


def checked_grid(slab, packet):
    if (not np.array_equal(packet["x"], slab.x) or
            not np.array_equal(packet["r"], slab.r) or float(packet["dz"]) != slab.dz):
        raise ValueError("Field grid differs from the requested grid")


def check_current(slab, center, currents):
    footprint = np.maximum(1-(slab.x[:, None]-center)**2-slab.r**2, 0)**3
    if currents.shape != (1, *footprint.shape):
        raise ValueError("Expected exactly one opening")
    current = currents[0]
    total = np.sum(slab.dv*current)
    if (not np.isfinite(total) or abs(total-0.25) > 3e-14 or
            np.any(current[footprint == 0] != 0) or np.any(current[footprint > 0] <= 0)):
        raise ValueError("Incorrect mouth position, radius, or rate")
    return current


def calibrate(args):
    slab = LateralSlab(*args.grid)
    field, currents, steps = slab.solve(0)
    check_current(slab, 0, currents)
    np.savez_compressed(args.output, layers=field[-4:], currents=currents,
                        x=slab.x, r=slab.r, dz=slab.dz, d=0)
    print(f"centered_calibration_steps={steps}")


def measure(args):
    slab = LateralSlab(*args.grid)
    path, digest = args.baseline
    with verified_packet(Path(path), {Path(path): digest}) as p:
        checked_grid(slab, p)
        if int(p["d"]) != 0:
            raise ValueError("Expected a centered calibration")
        target_current = check_current(slab, 0, p["currents"])
        density = abs(p["layers" if "layers" in p else "field"][-1])**2
    d = args.distance
    if args.archive:
        path, digest = args.archive
        with verified_packet(Path(path), {Path(path): digest}) as p:
            checked_grid(slab, p)
            if not np.array_equal(p["centers"], [-d]) or not np.array_equal(p["rates"], [0.25]):
                raise ValueError("Archived source preparation mismatch")
            layers, currents = p["layers"], p["currents"]
    else:
        layers, currents = checked_solve(slab, density, [-d], [0.25], args.output)
    check_current(slab, -d, currents)
    face = layers[-1]
    speed = -(face.conj()*dx4(face, slab.dx)).imag/abs(face)**2
    weights = slab.dv*target_current
    momentum_rate = float(np.sum(weights*speed))
    mean_speed = momentum_rate/float(weights.sum())
    if not np.isfinite(mean_speed):
        raise ArithmeticError("Nonfinite incident speed")
    print(json.dumps(dict(grid=args.grid, d=d, mean_speed=mean_speed,
                          incident_reference=momentum_rate)), flush=True)


def compare(args):
    values = {}
    for path in args.readouts:
        row = json.loads(path.read_text().splitlines()[-1])
        values_by_d = values.setdefault(tuple(row["grid"]), {})
        if row["d"] in values_by_d:
            raise ValueError("Duplicate readout")
        values_by_d[row["d"]] = row["mean_speed"]
    grids = [(48, 4, 16), (96, 4, 16), (96, 6, 24)]
    if set(values) != set(grids):
        raise ValueError("Missing domain or refinement control")
    ratios = []
    for grid in grids:
        v = values[grid]
        if set(v) != {4, 8, 16}:
            raise ValueError("Missing separation control")
        ratios.append(np.array([v[8]/v[4], v[16]/v[8]]))
        first, second = ratios[-1]
        print(f"grid={grid} U4={v[4]:.12g} U8={v[8]:.12g} U16={v[16]:.12g} "
              f"U8/U4={first:.12g} U16/U8={second:.12g} "
              f"longer_ratio_closer={abs(second-0.25) < abs(first-0.25)}")
    for label, delta in zip(("move_bath", "refine_mesh"), np.diff(ratios, axis=0)):
        print(f"{label} short_shift_pp={100*delta[0]:.9g} long_shift_pp={100*delta[1]:.9g}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(required=True)
    cal = sub.add_parser("calibrate")
    source = sub.add_parser("source")
    for cmd in (cal, source):
        cmd.add_argument("--grid", nargs=3, type=int, required=True, metavar=("L", "C", "Z"))
    cal.add_argument("--output", required=True, type=Path)
    cal.set_defaults(run=calibrate)
    source.add_argument("--baseline", nargs=2, required=True, metavar=("FILE", "SHA256"))
    source.add_argument("--distance", type=int, choices=(4, 8, 16), required=True)
    mode = source.add_mutually_exclusive_group(required=True)
    mode.add_argument("--archive", nargs=2, metavar=("FILE", "SHA256"))
    mode.add_argument("--output", type=Path)
    source.set_defaults(run=measure)
    summary = sub.add_parser("compare")
    summary.add_argument("readouts", nargs="+", type=Path)
    summary.set_defaults(run=compare)
    args = parser.parse_args()
    args.run(args)
