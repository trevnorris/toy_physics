"""Subtract the slow cylinder response from paired companion loads at d=4 and 8."""

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import jn_zeros

from lateral_companion_boundary_term import axial_kernel, radial_kernel, recorded_values


MEASUREMENTS = Path("docs/_measurements")
GRIDS = ("e24_c4", "e48_c4", "e24_c6")
RADII = (0.75, 1.0, 1.5)


def archived_runs(name):
    manifest = json.loads((MEASUREMENTS / f"lateral_{name}.json").read_bytes())
    for run in manifest["runs"]:
        log = ("\n".join(run["stdout"]) + "\n").encode()
        if (run["guard_exit"] or run["worker_exit"] or not run["limits_verified"]
                or run["swap_bytes"] != 0 or hashlib.sha256(log).hexdigest() != run["stdout_sha256"]):
            raise ValueError(f"Unverified archived run: {run['stdout_path']}")
    return manifest["runs"]


def corrected_load(paired, load, distance, extent, roots):
    radius = paired["radius"]
    rate = recorded_values(paired["stdout"], f"state_d={distance} ")
    setup = recorded_values(load["stdout"], f"radius={radius:g} d={distance} ")
    slow = recorded_values(load["stdout"], f"d={distance} radius={radius:g} omega=0 ")
    geometry = recorded_values(paired["stdout"], f"radius={radius:g} extent=")
    if (float(geometry["extent"]) != extent or geometry["dx"] != setup["dx"]
            or geometry["dz"] != setup["dz"]):
        raise ValueError("Readout and preparation geometry differ")
    shape = axial_kernel(distance/extent)
    independent = radial_kernel(distance/extent, roots, 10)
    if not np.isclose(shape, independent, rtol=0, atol=2e-12):
        raise ArithmeticError("Cylinder mode sums disagree")
    target_rate = float(rate["target_rate"])
    companion_rate = float(setup["companion_rate"])
    # Unit slab depth. Use the steady rate of this particular paired field.
    boundary = -target_rate * companion_rate * shape / extent
    observed = float(slow["cross_per_A"])
    direct = float(slow["direct"])
    if (companion_rate != 0.25 or target_rate <= 0 or
            not np.isclose(observed, float(slow["carry_per_A"]) + float(slow["stress_per_A"]),
                           rtol=1e-10, atol=1e-15)):
        raise ValueError("Invalid mouth rate or signed force components")
    return {
        "d": distance, "dx": float(setup["dx"]), "dz": float(setup["dz"]),
        "target_rate": target_rate, "mean_speed": float(rate["mean_speed"]),
        "companion_rate": companion_rate, "static_load": float(setup["static_companion_load"]),
        "kernel": shape, "mode_sum_difference": abs(shape-independent),
        "boundary": boundary, "observed": observed, "residual": observed-boundary,
        "direct": direct, "direct_residual": direct-boundary,
    }


def main():
    pair_runs = archived_runs("paired_acceleration")
    load_runs = archived_runs("companion_load")
    roots = jn_zeros(0, 512)
    pairs = []
    for grid in GRIDS:
        load, = (run for run in load_runs if run["grid"] == grid)
        for radius in RADII:
            paired, = (run for run in pair_runs if run["grid"] == grid and run["radius"] == radius)
            if any(load["reference_fields"].get(path) != sha
                   for path, sha in {**paired["reference_fields"], **paired["fields"]}.items()):
                raise ValueError("Rate and force refer to different paired fields")
            extent = int(grid.split("_")[0][1:])
            near, far = [corrected_load(paired, load, d, extent, roots) for d in (4, 8)]
            if near["mean_speed"] != far["mean_speed"]:
                raise ValueError("Target intake preparation changed between separations")
            pairs.append({
                "grid": grid, "L": extent, "target_radius": radius, "near": near, "far": far,
                "magnitude_ratio_8_over_4": abs(far["residual"]/near["residual"]),
                "direct_magnitude_ratio_8_over_4": abs(far["direct_residual"]/near["direct_residual"]),
                "static_ratio_8_over_4": far["static_load"]/near["static_load"],
            })

    changes = []
    for radius in RADII:
        small, large = [next(pair["far"] for pair in pairs
                             if pair["grid"] == grid and pair["target_radius"] == radius)
                        for grid in ("e24_c4", "e48_c4")]
        predicted = small["boundary"] - large["boundary"]
        observed = small["observed"] - large["observed"]
        changes.append({
            "target_radius": radius, "d": 8, "L": [24, 48], "predicted_change": predicted,
            "measured_change": observed, "direct_change": small["direct"]-large["direct"],
            "relative_difference_predicted_vs_measured": (predicted-observed)/observed,
        })
    print(json.dumps({"slab_depth": 1, "omega": 0, "companion_radius": 1,
                      "pairs": pairs, "boundary_changes": changes}, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
