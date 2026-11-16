"""Compare lateral and normal refinement at each fixed value of the other spacing."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np


def read_run(folder):
    outcome = json.loads((folder/"outcome.json").read_text())
    if (outcome["exitCode"] != 0 or outcome["childOutcome"]["exitCode"] != 0
            or not outcome["limitsVerified"]):
        raise ValueError(f"Unsuccessful guarded run in {folder}")
    log = (folder/"stdout").read_bytes()
    fields = dict(token.split("=") for token in log.decode().split())
    cells, normal = int(fields["cells"]), int(fields["normal_cells"])
    with np.load(folder/"mouth_layers.npz", allow_pickle=False) as packet:
        # Check the actual coordinates as well as the worker's labels.
        nr = 24*cells
        if (int(fields["extent"]) != 24 or float(packet["radius"]) != 1
                or float(packet["dz"]) != 1/normal
                or not np.array_equal(packet["x"], (np.arange(2*nr)+0.5-nr)*(1/cells))
                or not np.array_equal(packet["r"], (np.arange(nr)+0.5)*(1/cells))):
            raise ValueError(f"Different physical preparation or grid in {folder}")
    if abs(float(fields["rate"])-0.25) > 3e-14 or abs(float(fields["bath_inward"])-0.25) > 2e-8:
        raise ValueError(f"Different mouth or reservoir rate in {folder}")
    keys = ("lab", "whole_layer", "layer_cut", "layer_halo", "carry", "stress")
    values = {k: float(fields[k]) for k in keys}
    values["whole_minus_boundary"] = values["whole_layer"]-values["lab"]
    return (cells, normal), values, sha256(log).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("diagonal_manifest", type=Path)
    parser.add_argument("crossed_runs", type=Path, nargs=2)
    args = parser.parse_args()
    prior = json.loads(args.diagonal_manifest.read_text())
    grid = {}
    for record in prior["runs"]:
        key, values, digest = read_run(Path(record["stdout_path"]).parent)
        if digest != record["stdout_sha256"]:
            raise ValueError("Diagonal run differs from archived measurement")
        grid[key] = values
    for folder in args.crossed_runs:
        key, values, _ = read_run(folder)
        if key in grid:
            raise ValueError(f"Repeated grid {key}")
        grid[key] = values
    if set(grid) != {(8,32), (8,48), (12,32), (12,48)}:
        raise ValueError("Incomplete 2-by-2 grid")
    for (cells, normal), values in sorted(grid.items()):
        print(f"cells={cells} normal_cells={normal} " +
              " ".join(f"{k}={v:.12g}" for k,v in values.items()))
    for key in grid[8,32]:
        x32 = grid[12,32][key]-grid[8,32][key]
        x48 = grid[12,48][key]-grid[8,48][key]
        z8 = grid[8,48][key]-grid[8,32][key]
        z12 = grid[12,48][key]-grid[12,32][key]
        print(f"quantity={key} lateral_at_32={x32:.12g} lateral_at_48={x48:.12g} "
              f"normal_at_8={z8:.12g} normal_at_12={z12:.12g} interaction={x48-x32:.12g}")


if __name__ == "__main__":
    main()
