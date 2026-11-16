"""Compare fixed local flux and fixed normal-speed shape at the same total intake.

Reuse each saved zero-drive solution to calibrate the normal speed pointwise.
Only the two driven, redistributing states are new solves. The unaffected
three serialized states are read directly for baseline and fixed-flux data.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_feed_weak_flow import FeedPatch


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--reference-manifest", type=Path, required=True)
    parser.add_argument("--level", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.reference_manifest.read_text())
    hashes = {Path(path).resolve(): digest for run in manifest["runs"]
              for path, digest in run.get("fields", {}).items()}
    args.output.mkdir(parents=True, exist_ok=False)
    for law in ("flat_density_slope", "restoring_slope"):
        stem = f"n{args.level}_{law}_u"
        baseline_path = args.reference / (stem + "0.npz")
        if hashlib.sha256(baseline_path.read_bytes()).hexdigest() != hashes[baseline_path.resolve()]:
            raise ValueError("Baseline does not match the recorded run")
        with np.load(baseline_path) as packet:
            field = packet["field"]
            n0 = abs(field[-1])**2
            patch = FeedPatch(args.level, field.shape[0], law)
            if (field.shape != patch.shape or float(packet["dz"]) != patch.dz
                    or float(packet["dx"]) != patch.dx or float(packet["speed"]) != 0
                    or not np.array_equal(packet["current"], patch.current)
                    or not np.array_equal(packet["mouth"], patch.mouth) or np.min(n0) <= 0):
                raise ValueError("Incompatible zero-drive reference")
        del field
        for speed in (0.01, 0.02):
            field, steps, error = patch.solve(speed, reference_density=n0)
            result = patch.read(field, speed)
            # Re-evaluate the converged law from the saved face, avoiding
            # even the last fixed-point iteration lag in the output current.
            raw = patch.weight * abs(field[-1])**2 / n0
            current = 0.25 * raw / (patch.dx**3 * raw.sum())
            if (np.max(abs(current-patch.current)) > 2e-10 or result["min_density"] <= 0
                    or np.any(current[patch.mouth] <= 0) or np.any(current[~patch.mouth] != 0)
                    or abs(result["feed_rate"]-0.25) > 2e-7):
                raise ArithmeticError("Inward current, convergence or feed balance failed")
            np.savez_compressed(args.output/(stem + f"{speed:g}.npz"), field=field, current=current,
                                dx=patch.dx, dz=patch.dz, mouth=patch.mouth, speed=speed)
            print(f"n={args.level} depth={field.shape[0]-1} law={law} U={speed:g} "
                  f"steps={steps} update={error:.3g} feed={result['feed_rate']:.12g} "
                  f"mouth={patch.dx**3 * current.sum():.12g} min_density={result['min_density']:.9g}",
                  flush=True)


if __name__ == "__main__":
    main()
