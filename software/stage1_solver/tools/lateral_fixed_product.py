"""Measure the load on a held target for three pairs with the same rate product."""

import argparse
import json
from pathlib import Path

import numpy as np

from lateral_reservoir_field import LateralSlab
from lateral_source_rate import verified_packet
from lateral_target_rate import checked_solve, force


def archived_load(packet, slab, distance, rates):
    if (not np.array_equal(packet["x"], slab.x) or not np.array_equal(packet["r"], slab.r)
            or float(packet["dz"]) != slab.dz):
        raise ValueError("Archived field on a different grid")
    centers = [-distance, 0] if distance else [0]
    if "centers" in packet:
        if not np.array_equal(packet["centers"], centers) or not np.array_equal(packet["rates"], rates):
            raise ValueError("Archived mouth preparation mismatch")
    elif int(packet["d"]) != distance:
        raise ValueError("Archived separation mismatch")
    currents = packet["currents"]
    if len(currents) != len(rates):
        raise ValueError("Archived mouth count mismatch")
    for center, current, rate in zip(centers, currents, rates, strict=True):
        outside = (slab.x[:, None] - center)**2 + slab.r**2 >= 1
        if np.any(current[outside] != 0) or not np.isclose(
                np.sum(current*slab.dv), rate, rtol=0, atol=2e-13):
            raise ValueError("Archived intake location or rate mismatch")
    layers = packet["layers"] if "layers" in packet else packet["field"][-4:]
    return force(slab, layers, currents[-1], rates[-1])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path, help="Original equal-rate fields")
    parser.add_argument("isolated", type=Path, help="Isolated controls from the target-rate sweep")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    archive = json.loads(Path("docs/_measurements/lateral_target_rate.json").read_text())
    hashes = {Path(path): digest for run in archive["runs"] for key in ("reference_fields", "fields")
              for path, digest in run[key].items()}
    with verified_packet(args.reference / "d0.npz", hashes) as packet:
        r = packet["r"]
        slab = LateralSlab(round(r[-1]+r[0]), round(1/(2*r[0])), round(1/float(packet["dz"])))
        alone = {0.25: archived_load(packet, slab, 0, [0.25])}
        calibration = abs(packet["field"][-1])**2
    for target in (0.125, 0.5):
        with verified_packet(args.isolated / f"alone_J{target:g}.npz", hashes) as packet:
            alone[target] = archived_load(packet, slab, 0, [target])
    args.output.mkdir(exist_ok=False, parents=True)
    print(f"extent={r[-1]+r[0]:g} dx={slab.dx:.12g} dz={slab.dz:.12g}", flush=True)
    for target, parts in alone.items():
        print(f"  alone_Jt={target:g} carried={parts[0]:.12g} stress={parts[1]:.12g} "
              f"direct={parts[2]:.12g}", flush=True)
        if not np.isfinite(parts).all() or max(abs(parts)) > 1e-11:
            raise ArithmeticError("Invalid isolated target control")

    preparations = ((0.25, 0.25), (0.5, 0.125), (0.125, 0.5))
    loads = {}
    for distance in (4, 8):
        for source, target in preparations:
            if source == target:
                with verified_packet(args.reference / f"d{distance}.npz", hashes) as packet:
                    paired = archived_load(packet, slab, distance, [source, target])
            else:
                path = args.output / f"d{distance}_Js{source:g}_Jt{target:g}.npz"
                layers, currents = checked_solve(slab, calibration, [-distance, 0], [source, target], path)
                paired = force(slab, layers, currents[-1], target)
            carried, stress, direct = paired - alone[target]
            total = carried + stress
            if not np.isfinite([carried, stress, total, direct]).all():
                raise ArithmeticError("Nonfinite mouth readout")
            loads[distance, source, target] = (total, direct)
            print(f"  d={distance} Js={source:g} Jt={target:g} carried={carried:.12g} "
                  f"stress={stress:.12g} total={total:.12g} direct={direct:.12g}", flush=True)
        reference, reference_direct = loads[distance, 0.25, 0.25]
        if reference != 0 and reference_direct != 0:
            for source, target in preparations[1:]:
                total, direct = loads[distance, source, target]
                print(f"  d={distance} Js={source:g} Jt={target:g} "
                      f"change_pct={100*(total/reference-1):.9g} "
                      f"direct_change_pct={100*(direct/reference_direct-1):.9g}", flush=True)
    for source, target in preparations:
        near, _ = loads[4, source, target]
        far, _ = loads[8, source, target]
        if near != 0:
            print(f"  Js={source:g} Jt={target:g} F8_over_F4={far/near:.12g}", flush=True)


if __name__ == "__main__":
    main()
