"""Resolve the fixed-bond duct work from the tightly settled zero-displacement fields."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_closed_rim_surface_review import Reader
from mouth_last_ring_solid import interpolation
from mouth_bonded_front_adjoint import angular_weights
from mouth_rim_bank_work import Bank
from mouth_rim_disk_work_review import Disk
from mouth_rim_pinned_sources import Pinned, amplitude
from mouth_rim_pinned_work import pin_order_work
from mouth_supercore_joined_front import file_sha, state_from


def measure(state, wall, original, job, validated):
    meta = json.loads((job / "inputs.json").read_text())
    zero = job / "zero.npz"
    if (file_sha(zero) != validated["zero_state_sha256"] or
            file_sha(original) != meta["field_sha256"] or meta != validated["inputs"]):
        raise ValueError("Fields and geometry differ from the validated rematches")
    disk = Disk(state, wall, original)
    model = disk.model
    reference = state_from(state)[1]
    with np.load(zero, allow_pickle=False) as data:
        driven, resting = data["wound"], data["rest"]
    disk.reader = Reader(model, reference, driven, resting)
    edges = tuple(meta["duct_bounds"])
    knots = meta["terminal_knots"]
    half = meta["full_step"] / 2
    pins = Pinned(model, state, wall, reference, np.zeros(model.size))
    order = pin_order_work(disk, pins, 0, edges, knots, (-half, half))
    carrier = [Bank(disk, 0, knots, edge).work(2)["native_including_vortex_bank"] for edge in edges]
    carrier_difference = np.subtract(*carrier)
    r, z, weight = pins.material_sites()
    radius, slope = disk.sh.radius(z)
    x = r / radius
    order_difference = (driven[2] - resting[2])[model.ids[1, 0]]
    mean = np.einsum("raj,a->jr", order_difference, angular_weights(model.a))
    _, dx, dz = interpolation(model.z[1], model.rr[0], mean, z, x)
    crystal = weight @ (amplitude(0, edges, knots, x, z) * (dz - x*slope*dx/radius))
    predicted = carrier_difference[-1] + crystal + order.sum(axis=1).mean()
    measured = validated["four_minus_eight"]["duct_relative"]["half"]["centered"]
    return dict(grid=meta["grid"], zero_sha256=file_sha(zero), carrier=carrier_difference.tolist(),
                order_static_and_phase=order.tolist(), crystal=crystal, predicted=predicted,
                settled=measured, settled_minus_predicted=measured-predicted)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("state", "wall", "original", "job", "validated"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    cases = json.loads(a.validated.read_text())["cases"]
    case = next(c for c in cases if c["inputs"]["grid"] == a.job.name)
    print(json.dumps(measure(a.state, a.wall, a.original, a.job, case), allow_nan=False))
