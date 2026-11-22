"""Signed fixed-action shape work of both shear laws on the same frozen broad-waist fields."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_exchange_shear import DEGREE, solve
from mouth_live_waist_wave import OccupiedWave
from mouth_sealed_shape_relax import Meridian
from mouth_unified_transverse import ARCHIVE, PORT, cell_difference


def measure(center, archive):
    case = center["case"]
    prep = json.loads(archive["runs"][case]["literal_stdout"].splitlines()[0])
    source = np.load(center["snapshot"])
    n, mesh, radius = prep["mesh"], prep["field_mesh"], prep["radius"]
    coords = np.array(prep["coordinates"])
    if radius not in (1, 2) or prep["point"] != "broad" or prep["J"] != 1:
        raise ValueError("Unexpected source preparation")
    families = ("curl", "symmetric")
    result = {}
    for axis, index, scale in (("waist", 2, 1/radius), ("log_depth", 0, 1.)):
        pairs = {key: [] for key in families}
        prior = {}
        for delta in (.002, .001):
            branches = {key: {} for key in families}
            for sign in (-1, 1):
                change = sign*delta
                x = coords.copy()
                x[index] += scale*change
                shape = Meridian(radius, x)
                wave = OccupiedWave(shape, n, 2*n, mesh, PORT)
                weights = wave.project @ source["order"].T.ravel()
                inside = wave.matrix(weights[:, None, None]*wave.mc)
                mass = wave.cut_m + inside
                curl = wave.cut_k + wave.matrix(weights[:, None, None]*wave.kc)
                diff = wave.matrix(weights[:, None, None]*cell_difference(wave, shape, mesh, DEGREE))
                for fraction, family in enumerate(families):
                    endpoint = center["endpoints"][fraction]
                    omega = complex(endpoint["frequency"], -endpoint["amplitude_loss"])
                    seed = prior.get((family, sign), source["vector"])
                    pencil = (wave.uniform[0], curl-fraction*diff, mass, inside)
                    out, vec, _ = solve(shape, mesh, PORT, 0., omega, incoming_vector=seed,
                                        supplied_system=pencil, target_frequency=omega.real, shape_work=False)
                    overlap = abs(np.vdot(seed, mass@vec))**2/np.vdot(seed, mass@seed).real
                    if overlap < .95 or abs(out["frequency"]/omega.real-1) > .005:
                        raise ArithmeticError(f"Lost signed branch on {case} {axis} {change} {family}")
                    branches[family][sign] = (out["frequency"], shape.length, overlap)
                    prior[family, sign] = vec
            for fraction, family in enumerate(families):
                lo, hi = branches[family][-1], branches[family][1]
                denominator = 2*delta if axis == "waist" else hi[1]-lo[1]
                center_frequency = center["endpoints"][fraction]["frequency"]
                pairs[family].append(dict(offset=delta, low_frequency=lo[0], high_frequency=hi[0],
                                          center_frequency=center_frequency, low_depth=lo[1], high_depth=hi[1],
                                          force=-(hi[0]-lo[0])/denominator,
                                          squared_overlap_min=min(lo[2], hi[2])))
        for rows in pairs.values():
            f1, f2 = [row["force"] for row in rows]
            if abs(f2-f1) > .01*max(abs(f2), 1e-8):
                raise ArithmeticError(f"Signed force changed by over 1% on halving step for {case} {axis}")
        result[axis] = pairs
    print(json.dumps(dict(case=case, frozen_material_order_sha256=center["sha256"],
                          action=prep["J"], signed_wave_forces=result), allow_nan=False), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("centers", type=Path, help="JSON lines produced by mouth_unified_transverse.py")
    args = parser.parse_args()
    archive = json.loads(ARCHIVE.read_bytes())
    for line in args.centers.read_text().splitlines():
        measure(json.loads(line), archive)
