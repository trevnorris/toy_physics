"""Check whether the saved brane-side states supply a free translation experiment.

The translation control changes the input geometry by one grid cell. It
tests the centroid and the stationary boundary residual, with no physical
time advancement. Length, force and current use the held experiment's units.
"""

from hashlib import sha256
import json
from pathlib import Path

import numpy as np


def manifest_fields(path):
    data = json.loads(path.read_text())
    return {Path(name): digest for run in data["runs"]
            for name, digest in run.get("fields", {}).items()}


def centroid(weights, x):
    return float(np.sum(weights*x) / weights.sum())


def main():
    docs = Path("docs/_measurements")
    baseline_hashes = manifest_fields(docs/"mouth_feed_weak_run.json")
    driven_hashes = manifest_fields(docs/"mouth_intake_redistribution.json")
    baseline_dir = Path("_scratch/s11c/mouth_feed_depth01/fields")
    driven_dir = Path("_scratch/s11c/mouth_redistribution_depth01/fields")
    for law in ("flat_density_slope", "restoring_slope"):
        for speed in (0, 0.01, 0.02):
            folder, hashes = ((baseline_dir, baseline_hashes) if speed == 0 else
                              (driven_dir, driven_hashes))
            path = folder/f"n32_{law}_u{speed:g}.npz"
            if sha256(path.read_bytes()).hexdigest() != hashes[path]:
                raise ValueError(f"Invalid source packet {path}")
            with np.load(path) as pkt:
                dx, dz = float(pkt["dx"]), float(pkt["dz"])
                field, mouth, current = pkt["field"], pkt["mouth"], pkt["current"]
                if float(pkt["speed"]) != speed:
                    raise ValueError("Drive does not match the source packet")
            coord = (np.arange(mouth.shape[0]) + 0.5 - mouth.shape[0]/2)*dx
            x = coord[:, None, None]
            rate = dx**3 * current.sum()
            feed_rate = dx**3 * np.imag(field[0].conj()*field[1]).sum()/dz
            center, intake_center = centroid(mouth, x), centroid(current, x)
            if speed == 0:
                initial_mouth = mouth.copy()
                initial_intake_center = intake_center
            changed = np.count_nonzero(mouth != initial_mouth)
            print(f"{law} U={speed:g} cells_changed={changed} mouth_x={center:.12g} "
                  f"intake_x_shift={intake_center-initial_intake_center:.12g} "
                  f"mouth_rate={rate:.12g} feed_rate={feed_rate:.12g}", flush=True)
            if changed or abs(rate-0.25) > 1e-12 or abs(feed_rate-rate) > 2e-7:
                raise ArithmeticError("Saved geometry or drainage disagrees with the preparation")
            if speed != 0:
                continue
            # Discrete translation of the whole stationary solution on a
            # periodic brane is a relabelling. The old and translated local
            # normal-current conditions coincide only if the sink moves too.
            shifted_current = np.roll(current, 1, axis=0)
            shifted_mask = np.roll(mouth, 1, axis=0)
            dx_test = centroid(shifted_mask, x) - center
            # Reconstruct the face current from the carrier on the feed
            # side, without importing the imposed Neumann derivative.
            phase_flux = np.imag(field[-1].conj() *
                                 (11*field[-1]-18*field[-2]+9*field[-3]-2*field[-4]))/(6*dz)
            translated_flux = np.roll(phase_flux, 1, axis=0)
            native_error = np.max(abs(phase_flux-current))
            fixed_error = np.max(abs(translated_flux-current))
            shifted_error = np.max(abs(translated_flux-shifted_current))
            print(f"{law} prescribed translation={dx:g} centroid_change={dx_test:.12g} "
                  f"normal_flux_max_mismatch(original,translated_field_only,translated_together)="
                  f"{native_error:.9g},{fixed_error:.9g},{shifted_error:.9g}", flush=True)
            if abs(dx_test-dx) > 1e-13 or shifted_error != native_error or fixed_error <= native_error:
                raise ArithmeticError("The geometric translation control failed")


if __name__ == "__main__":
    main()
