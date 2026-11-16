"""Slow exterior boundary contribution for a centered translating point drain.

The reservoir cylinder has unit radius and endcaps at X=±1 in scaled
coordinates. K(s) is the regular axial mixed derivative of (-Δ_D)^(-2)
at X=-s, Y=0. The singular free-space kernel has zero axial mixed derivative
at distinct collinear points. K(0) denotes the regular limit along the axis.
"""

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import i0e, j1, jn_zeros


ARCHIVE = Path("docs/_measurements/lateral_companion_domain.json")


def axial_kernel(s, count=24):
    # Images of the screened free-space squared inverse between the two caps
    # give -1/(8π) after taking the zero-screening limit. Add the sidewall
    # using the remaining even Dirichlet sine modes on [-1,1].
    k = np.pi * np.arange(1, count + 1)
    inverse_i0 = np.exp(-k) / i0e(k)
    wall = -np.sum(np.cos(k*s) * inverse_i0**2) / (4*np.pi)
    return -1/(8*np.pi) + wall


def radial_kernel(s, roots, images):
    # Independent order of separation. Each transverse J0 mode produces a
    # one-dimensional squared inverse (1+j|X-Y|) exp(-j|X-Y|)/(4j³).
    # Differentiating in X and Y gives both image families the same sign.
    if not 0 < s < 1:
        raise ValueError("The radial sum requires distinct interior points")
    weights = 1 / (4*np.pi * roots * j1(roots)**2)
    distances = np.abs(s - 2*np.arange(-images, images+1))
    arguments = roots[:, None] * distances
    return np.sum(weights[:, None] * (1-arguments) * np.exp(-arguments))


def recorded_values(lines, prefix):
    line, = (line for line in lines if line.startswith(prefix))
    return dict(token.split("=") for token in line.split())


def main():
    archive_bytes = ARCHIVE.read_bytes()
    archive = json.loads(archive_bytes)
    large = max(archive["runs"], key=lambda run: run["L"])
    inputs = recorded_values(large["stdout"], "state_d=")
    setup = recorded_values(large["stdout"], "L=")
    if large["guard_exit"] or large["worker_exit"] or not large["limits_verified"]:
        raise ValueError("Invalid archived calibration run")
    distance = float(setup["d"])
    if distance != float(inputs["state_d"]) or float(setup["target_radius"]) != 1.5:
        raise ValueError("Archived mouth geometry mismatch")
    target_rate = float(inputs["target_rate"])
    companion_rate = float(setup["companion_rate"])
    depth = 1.0
    load_factor = -target_rate * companion_rate / depth
    central = axial_kernel(0.0)
    slope = load_factor * central
    result = {
        "domain_input": str(ARCHIVE),
        "domain_input_sha256": hashlib.sha256(archive_bytes).hexdigest(),
        "calibration_L": large["L"], "target_rate": target_rate,
        "companion_rate": companion_rate, "separation": distance, "depth": depth,
        "regular_axial_kernel_at_center": central,
        "endcaps_kernel": -1/(8*np.pi),
        "sidewall_kernel_at_center": central + 1/(8*np.pi),
        "predicted_inverse_L_coefficient": slope,
    }

    roots = jn_zeros(0, 1024)
    rows = []
    for row in sorted(archive["comparisons"], key=lambda row: row["L"]):
        if row["omega"] != 0:
            continue
        extent = row["L"]
        s = distance / extent
        kernel = axial_kernel(s)
        radial = radial_kernel(s, roots, 12)
        errors = [abs(kernel-radial), abs(kernel-axial_kernel(s, 12)),
                  abs(kernel-radial_kernel(s, roots[:512], 8))]
        if not np.all(np.isfinite(errors)) or max(errors) > 2e-12:
            raise ArithmeticError(f"Cylinder kernel checks failed for {extent=}: {errors}")
        boundary = load_factor * kernel / extent
        rows.append({
            "L": extent, "kernel": kernel, "radial_crosscheck_abs_error": errors[0],
            "max_crosscheck_abs_error": max(errors),
            "predicted_boundary": boundary, "measured": row["cross"], "direct": row["direct"],
            "measured_minus_boundary": row["cross"]-boundary,
            "direct_minus_boundary": row["direct"]-boundary,
        })
    if [row["L"] for row in rows] != [24, 48, 72, 96, 128]:
        raise ValueError("Incomplete or repeated slow-limit comparison")
    steps = []
    for small, large in zip(rows, rows[1:]):
        inverse_step = 1/small["L"] - 1/large["L"]
        predicted = small["predicted_boundary"] - large["predicted_boundary"]
        observed = small["measured"] - large["measured"]
        direct = small["direct"] - large["direct"]
        steps.append({
            "sizes": [small["L"], large["L"]],
            "predicted_change_small_minus_large": predicted,
            "measured_change_small_minus_large": observed,
            "direct_change_small_minus_large": direct,
            "measured_inverse_L_slope": observed/inverse_step,
            "direct_inverse_L_slope": direct/inverse_step,
            "leading_slope_relative_difference": slope/(observed/inverse_step)-1,
            "finite_offset_relative_difference": predicted/observed-1,
            "direct_finite_offset_relative_difference": predicted/direct-1,
        })
    result.update(rows=rows, adjacent_sizes=steps)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
