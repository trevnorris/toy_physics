"""First-bend wall displacement on the released patch, in the lateral wall tangent."""
import argparse
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator

from mouth_carrier_lateral_inertia import taper
from mouth_connected_passage import MIDDLE, WAIST
from mouth_empty_bend_stiffness import INTERVALS, shapes
from mouth_middle_tag_grip import cubic_trace
from mouth_middepth_rearrangement import RECORD, wall_table
from mouth_order_wave_pulse import digest, emit
from mouth_segment_contact_hold import AMPLITUDE, PITCH
from neutral_patch_sites import COUNT, WITHIN
from neutral_patch_turning_neighbor import ENTRY

PREVIOUS = Path("research/pde_ledger_v3/_measurements/mouth_middepth_rearrangement.json")


def run(args):
    prior = json.loads(PREVIOUS.read_text())
    for file, expected in prior["input_sha256"].items():
        if digest(Path(file)) != expected:
            raise ValueError(f"Changed patch input {file}")
    reference = [json.loads(line) for line in prior["guarded_runs"]["morse"]["literal_stdout"]]
    setup = next(row for row in reference if row["test"] == "local_setup")
    sections = setup["reference_gap"]
    period = setup["time_per_section"]
    rows = [json.loads(line) for line in json.loads(RECORD.read_text())["guarded_runs"]["resolution"]["literal_stdout"]]
    tangents = next(row["neutral_tangents"] for row in rows if row["test"] == "setup")
    trace, intake = cubic_trace(34)
    np.testing.assert_allclose(intake, .25, rtol=1e-10)

    # Map the native bare section spacing into arclength along the actual cubic wall.
    # Include the taper end and gauge center as interpolation knots.
    low, high = INTERVALS[0]
    depth = np.unique(np.r_[np.linspace(WAIST, 2., 8193), low, high, MIDDLE])
    arc = cumulative_trapezoid(np.hypot(1., trace.geometry.radius(depth)[1]), depth, initial=0)
    registry = (arc - arc[np.searchsorted(depth, MIDDLE)]) / (sections * PITCH)
    to_depth = PchipInterpolator(registry, depth, extrapolate=False)
    table_size = 2 * args.sub
    # Half-row stagger and half-time integration steps share the same spatial grid.
    coord = np.arange((args.length + 2) * table_size + 1) / table_size - args.length / 2
    depths = to_depth(coord)
    def wall_shape(w):
        w = np.atleast_1d(w)
        radius = trace.geometry.radius(w)[0]
        first = shapes(radius, w, "moving-box", trace.geometry.far, trace.geometry.roof)[0][:, 0]
        np.testing.assert_allclose(first, 1 - taper(w, low, high)[0], atol=1e-14)
        if args.shape == "uniform":
            return np.ones_like(first)
        return 1 - first if args.shape == "entrance-share" else first

    values = wall_shape(depths)
    if not np.all(np.isfinite(values)) or not np.all((values >= -1e-12) & (values <= 1. + 1e-12)):
        raise ValueError("Patch extends outside the mapped wall")
    end_values = {"first-bend": [0., 1.], "uniform": [1., 1.], "entrance-share": [1., 0.]}
    np.testing.assert_array_equal(wall_shape([low, MIDDLE]), end_values[args.shape])

    def shape(x):
        return float(wall_shape(to_depth(x))[0])

    emit(test="geometry", shape=args.shape, shape_interval=[low, high], radius=1., intake=intake,
         middle=MIDDLE, gauge_depths=to_depth([-6., 6.]).tolist(),
         gauge_shape=[shape(-6.), shape(0.), shape(6.)], end_in_section_spacings=float(registry[np.searchsorted(depth, high)]),
         left_depth=float(depths[0]), left_shape=float(values[0]), pitch=PITCH,
         amplitude_in_mouth_radii=args.amplitude, amplitude_in_bond_lengths=args.amplitude/PITCH,
         feed=args.feed, length=args.length, lanes=args.lanes, sub=args.sub, seed=args.seed,
         period=period, prior_sha256=digest(PREVIOUS), binary_sha256=digest(args.binary))
    potential, alpha, rest = WITHIN
    params = [args.length, args.lanes, args.sub, args.seed, "all", sections, setup["hoop"], period,
              setup["order"], setup["contact_per_section"], args.amplitude/PITCH, rest, alpha,
              potential, ENTRY, COUNT, setup["cutoff"], 8, 8, 1]
    with tempfile.TemporaryDirectory() as folder:
        contacts, profile = [Path(folder)/name for name in ("contacts", "profile")]
        wall_table(contacts, tangents)
        np.savetxt(profile, values, fmt="%.17g")
        subprocess.run([str(args.binary.resolve()), *map(str, params), str(contacts), str(profile),
                        args.feed, str(registry[np.searchsorted(depth, high)])], check=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--binary", required=True, type=Path)
    p.add_argument("--length", type=int, default=32, choices=[32,48,64,96])
    p.add_argument("--lanes", type=int, default=4, choices=[4,6,8])
    p.add_argument("--sub", type=int, default=128, choices=[128,256])
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--amplitude", type=float, default=AMPLITUDE/2)
    p.add_argument("--shape", choices=["first-bend","uniform","entrance-share"], default="first-bend")
    p.add_argument("--feed", choices=["centered","following"], default="centered")
    run(p.parse_args())
