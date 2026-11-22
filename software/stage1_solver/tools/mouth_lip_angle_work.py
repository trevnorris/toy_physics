"""Test restoring work invisible to any stationary rim sensing only aperture and tilt."""
import argparse
import hashlib
import json
from math import fsum, log1p
from pathlib import Path

import numpy as np

from mouth_sealed_shape_relax import Meridian
from mouth_shallow_puncture_fields import fields, occupied, read_case


def shallow_center(path, radius, n, original_hash):
    outcome = json.loads((path.parent / "outcome.json").read_text())
    checksum = hashlib.sha256(path.read_bytes()).hexdigest()
    if not outcome["limitsVerified"] or outcome["exitCode"] or outcome["childOutcome"]["exitCode"]:
        raise ValueError("Shallow state lacks a successful guard outcome")
    if checksum not in (path.parent / "stdout").read_text():
        raise ValueError("Missing shallow checkpoint fingerprint in validated output")
    with np.load(path) as data:
        if (str(data["original_snapshot_sha256"].item()) != original_hash
                or data["order"].shape != (2*n+1, n+1)):
            raise ValueError("Wrong material preparation")
        x = data["coordinates"]
        if abs(Meridian(radius, x).length - 1) > 1e-12:
            raise ValueError("Checkpoint does not have unit depth")
        return x, data["order"], data["vector"], complex(data["pole"].item()), checksum


def run(radius, n, shallow_path):
    prep, chi, seed, ref, _, exterior_mesh, original_hash = read_case(radius, n)
    x = np.array(prep["coordinates"])
    source = original_hash
    pole = complex(ref["frequency"], -ref["amplitude_loss"])
    if shallow_path is not None:
        x, chi, seed, pole, source = shallow_center(shallow_path, radius, n, original_hash)
    mesh = prep["field_mesh"]
    center = Meridian(radius, x)
    slope = float(center.profile(0.)[1])
    depths = (-.002, -.001, .001, .002)
    samples = {}
    for step in depths:
        # b'_new(0) = (L*b'_old(0)+4*dq)/(L+dL).
        shifted = x.copy()
        shifted[0] += log1p(step/center.length)
        shifted[2] += slope*step/(4*radius)
        trial = Meridian(radius, shifted)
        rim, derivative, _ = trial.profile(0.)
        if abs(rim-radius) > 1e-12 or abs(derivative-slope) > 1e-11:
            raise ArithmeticError("The virtual displacement moved or rotated the rim")
        _, _, wave = occupied(trial, n, mesh, chi, pole, seed)
        energies, metrics = fields(trial, n, mesh, exterior_mesh, chi, wave)
        if abs(metrics["current"]-metrics["terminal"]) > 1e-8:
            raise ArithmeticError("Intake and exit disagree on the trial")
        samples[step] = energies
        print(json.dumps(dict(radius=radius, mesh=n, center_depth=center.length,
                              depth_step=step, radial_step=slope*step/4, slope=derivative,
                              base_input_sha256=source, energies=energies, **metrics),
                         allow_nan=False), flush=True)
    for low, high in ((-.002, -.001), (.001, .002), (-.001, .001), (-.002, .002)):
        a, b = samples[low], samples[high]
        forces = {k:(a[k]-b[k])/(high-low) for k in a}
        print(json.dumps(dict(radius=radius, mesh=n, center_depth=center.length,
                              interval=(low, high), components=forces,
                              force_at_fixed_lip=fsum(forces.values())),
                         allow_nan=False), flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--radius", type=int, choices=(1, 2), required=True)
    p.add_argument("--mesh", type=int, choices=(48, 72), required=True)
    p.add_argument("--shallow-center", type=Path)
    args = p.parse_args()
    run(args.radius, args.mesh, args.shallow_center)
