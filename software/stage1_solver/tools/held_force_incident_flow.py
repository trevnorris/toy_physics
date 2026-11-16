"""Measure how much of a single drain's incident momentum rate appears in a held target load."""

import argparse
import json
from pathlib import Path

import numpy as np

from incident_flow_falloff import check_current, checked_grid
from lateral_reservoir_field import LateralSlab
from lateral_source_rate import verified_packet
from lateral_target_rate import checked_solve, force
from mouth_feed_weak_check import dx4


RATE = 0.25
GRIDS = ((48, 4, 16), (96, 4, 16), (96, 6, 24))
ARCHIVE = Path("docs/_measurements/incident_flow_falloff.json")


def packet(reference):
    path, digest = reference
    return verified_packet(Path(path), {Path(path): digest})


def measure(args):
    slab = LateralSlab(*args.grid)
    if tuple(args.grid) not in GRIDS:
        raise ValueError("No matched incident control on this grid")
    d = args.distance
    with packet(args.baseline) as p:
        checked_grid(slab, p)
        if int(p["d"]) != 0:
            raise ValueError("Expected a centered isolated target")
        target = check_current(slab, 0, p["currents"])
        layers = p["layers"] if "layers" in p else p["field"][-4:].copy()
    density = abs(layers[-1])**2
    alone = force(slab, layers, target, RATE)
    print(f"alone_carried={alone[0]:.12g} alone_stress={alone[1]:.12g} "
          f"alone_direct={alone[2]:.12g}", flush=True)
    if not np.isfinite(alone).all() or max(abs(alone)) > 1e-11:
        raise ArithmeticError("Stationary isolated target failed reflection check")
    with packet(args.source) as p:
        checked_grid(slab, p)
        if not np.array_equal(p["centers"], [-d]) or not np.array_equal(p["rates"], [RATE]):
            raise ValueError("Expected only the companion in the incident control")
        check_current(slab, -d, p["currents"])
        face = p["layers"][-1]
        incident_speed = -(face.conj()*dx4(face, slab.dx)).imag/abs(face)**2
        available = np.sum((slab.dv*target)*incident_speed)
        mean_speed = available/np.sum(slab.dv*target)
    if args.pair:
        with packet(args.pair) as p:
            checked_grid(slab, p)
            if int(p["d"]) != d:
                raise ValueError("Archived pair has the wrong separation")
            layers, currents = p["layers"], p["currents"]
    else:
        if args.output.exists():
            raise FileExistsError(args.output)
        layers, currents = checked_solve(slab, density, [-d, 0], [RATE, RATE], args.output)
    if len(currents) != 2:
        raise ValueError("Expected two active mouths")
    check_current(slab, -d, currents[:1])
    paired_target = check_current(slab, 0, currents[1:])
    carried, stress, direct = force(slab, layers, paired_target, RATE)-alone
    total = carried+stress
    if not np.isfinite([total, direct, available, mean_speed]).all() or available == 0:
        raise ArithmeticError("Unusable load or incident speed")
    print(json.dumps(dict(grid=args.grid, d=d, incident_speed=mean_speed,
                          incident_momentum_rate=available, carried=carried, stress=stress,
                          force=total, coefficient=total/available, direct=direct,
                          direct_coefficient=direct/available)), flush=True)


def compare(args):
    archive = {r["name"]: r for r in json.loads(ARCHIVE.read_text())["runs"]}
    rows = {}
    for path in args.readouts:
        row = json.loads(path.read_text().splitlines()[-1])
        extent, cells, _ = row["grid"]
        d = row["d"]
        old = archive[f"L{extent}C{cells}_d{d}"]
        ref = json.loads(old["stdout"][-1])
        if (old["guard_exit"] or old["worker_exit"] or not old["limits_verified"] or
                ref["d"] != d or ref["grid"] != row["grid"] or
                not np.isclose(row["incident_momentum_rate"], ref["incident_reference"], rtol=1e-12, atol=0) or
                not np.isclose(row["incident_speed"], ref["mean_speed"], rtol=1e-12, atol=0)):
            raise ValueError("Reconstructed incident reference disagrees with the measured control")
        key = tuple(row["grid"]), d
        if key in rows:
            raise ValueError("Duplicate held-mouth readout")
        rows[key] = row
    if set(rows) != {(grid, d) for grid in GRIDS for d in (8, 16)}:
        raise ValueError("Incomplete matched controls")
    for grid in GRIDS:
        low, high = (rows[grid, d] for d in (8, 16))
        def ratio(key):
            return high[key]/low[key]
        print(f"L={grid[0]} c={grid[1]} force16_over_8={ratio('force'):.12g} "
              f"speed16_over_8={ratio('incident_speed'):.12g} K8={low['coefficient']:.12g} "
              f"K16={high['coefficient']:.12g} K_relative_change_ppm={1e6*(ratio('coefficient')-1):.9g} "
              f"direct16_over_8={ratio('direct'):.12g} "
              f"direct_K_relative_change_ppm={1e6*(ratio('direct_coefficient')-1):.9g} "
              f"direct_minus_primary_at8_percent={100*(low['direct']/low['force']-1):.9g} "
              f"direct_minus_primary_at16_percent={100*(high['direct']/high['force']-1):.9g}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(required=True)
    run = sub.add_parser("measure")
    run.add_argument("--grid", nargs=3, type=int, required=True)
    run.add_argument("--distance", type=int, choices=(8, 16), required=True)
    for flag in ("baseline", "source"):
        run.add_argument("--"+flag, nargs=2, metavar=("FILE","SHA256"), required=True)
    mode = run.add_mutually_exclusive_group(required=True)
    mode.add_argument("--output", type=Path)
    mode.add_argument("--pair", nargs=2, metavar=("FILE","SHA256"))
    run.set_defaults(run=measure)
    comp = sub.add_parser("compare")
    comp.add_argument("readouts", nargs="+", type=Path)
    comp.set_defaults(run=compare)
    args = parser.parse_args()
    args.run(args)
