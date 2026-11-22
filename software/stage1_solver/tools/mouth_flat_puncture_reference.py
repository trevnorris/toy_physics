"""Finite-strain geometric tests on stationary axisymmetric material registrations."""
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, quad
from scipy.optimize import brentq, minimize_scalar

from mouth_sealed_shape_relax import Meridian


ROOT = Path(__file__).resolve().parents[3]
RECORD = ROOT / "research/pde_ledger_v3/_measurements/mouth_live_waist_order.json"
THICKNESS = 1.  # Reference ordered slab, same flat preparation as the stretch test.


def read_preparation(manifest, radius):
    key = f"r{radius}_broad_72"
    archived = manifest["runs"][key]
    original = archived["literal_stdout"]
    if (archived["worker_exit"] or archived["guard_exit"] or not archived["limits_verified"]
            or hashlib.sha256(original.encode()).hexdigest() != archived["stdout_sha256"]):
        raise ValueError(f"Unverified held geometry {key}")
    prep, = (row for row in map(json.loads, original.splitlines()) if row["kind"] == "preparation")
    if (prep["radius"], prep["point"], prep["J"]) != (radius, "broad", 1):
        raise ValueError("Unexpected radius, state or occupied wave")
    return Meridian(radius, prep["coordinates"]), archived["stdout_sha256"]


def measure(shape):
    end = shape.length

    def geometry(w):
        return shape.profile(w / end)

    def radius(w):
        return float(geometry(w)[0])

    def slope(w):
        return float(geometry(w)[1])

    def arc(lo, hi):
        return quad(lambda w: np.hypot(1., slope(w)), lo, hi,
                    epsabs=1.e-10, epsrel=1.e-10)[0]

    # The analytic meridian on the two held profiles has one interior minimum.
    grid = np.linspace(0., end, 1025)
    radii, slopes, _ = geometry(grid)
    changes = np.flatnonzero(slopes[:-1] * slopes[1:] < 0)
    if len(changes) != 1 or slopes[0] >= 0 or slopes[-1] <= 0:
        raise ValueError("Unexpected extrema in the held profile")
    i = changes.item()
    waist = brentq(slope, grid[i], grid[i+1])
    narrow = radius(waist)
    exit_radius = radius(end)
    # Flat face: material radius S decreases in the direction going in from
    # the exterior. lambda_h=r/S in both angular directions; lambda_s=ds/|dS|.
    # If all lie in [1/B,B], then length(u,v) <= B(S_u-S_v)
    # <= B^2*r_u - r_v. For any given u the numerator length(u,v)+r_v
    # grows with v because |dr/ds|<=1, so the terminal v is strongest.
    lengths = cumulative_trapezoid(np.hypot(1., slopes), grid, initial=0.)
    ratios = (lengths[-1] - lengths + exit_radius) / radii
    j = int(np.argmax(ratios))

    def root_required(w):
        return np.sqrt((arc(w, end) + exit_radius) / radius(w))

    extremes = [0., end]
    if 0 < j < len(grid) - 1:
        refined = minimize_scalar(lambda w: -root_required(w),
                                  bounds=(grid[j-1], grid[j+1]), method="bounded",
                                  options={"xatol": 1.e-10})
        if not refined.success:
            raise ArithmeticError("Material stretch bound did not converge")
        extremes.append(float(refined.x))
    upstream = max(extremes, key=root_required)
    bound = root_required(upstream)
    arc_total = arc(0., end)
    # Different preparation: the whole traced meridian was the cylindrical
    # bore through a slab of thickness H, with one fixed material hole radius.
    # Even optimal initial radius cannot remove the angular stretch spread.
    max_radius = max(radius(0.), exit_radius)
    bore_hoop_bound = np.sqrt(max_radius / narrow)
    bore_meridian_bound = max(arc_total/THICKNESS, THICKNESS/arc_total)
    return dict(depth=end, arc_length=arc_total, thickness=THICKNESS,
                waist_depth=waist, waist_radius=narrow, exit_radius=exit_radius,
                waist_angular_sectional_curvature=1./(narrow*narrow*(1+slope(waist)**2)),
                flat_face=dict(witness_upstream_depth=upstream, witness_upstream_radius=radius(upstream),
                               remaining_arc=arc(upstream, end), endpoint_radius=exit_radius,
                               reciprocal_stretch_lower_bound=bound, log_strain_lower_bound=np.log(bound),
                               independent_grid_bound=np.sqrt(np.max(ratios)),
                               waist_to_exit_bound=root_required(waist)),
                straight_bore=dict(angular_bound=bore_hoop_bound, meridional_bound=bore_meridian_bound,
                                   largest_reciprocal_stretch_lower_bound=max(bore_hoop_bound,bore_meridian_bound)))


def main():
    manifest = json.loads(RECORD.read_text())
    for radius in (1, 2):
        shape, source = read_preparation(manifest, radius)
        record = measure(shape)
        print(json.dumps(dict(entrance_radius=radius, source_stdout_sha256=source, **record),
                         allow_nan=False), flush=True)
        face = record["flat_face"]
        if (abs(face["independent_grid_bound"] - face["reciprocal_stretch_lower_bound"]) > 2.e-5
                or face["reciprocal_stretch_lower_bound"] < face["waist_to_exit_bound"] - 1.e-8):
            raise ArithmeticError("Integrated meridian and sampled geometry disagree")


if __name__ == "__main__":
    main()
