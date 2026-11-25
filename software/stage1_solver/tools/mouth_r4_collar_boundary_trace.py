"""Check one-sided natural trace tangents of the registered radius-four collar."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_r4_aligned_derivatives import load
from mouth_r4_depth_cell_geometry import digest


def measure(line, coordinate):
    # Trace within each registry rectangle is exactly affine on the edge in
    # (eta,t). The source adds 4 to the natural radial argument.
    slopes = np.diff(line, axis=0) / np.diff(coordinate)[:, None]
    speed = np.sqrt(np.sum(slopes * slopes, axis=1))
    unit = slopes / speed[:, None]
    jumps = np.linalg.norm(np.diff(slopes, axis=0), axis=1)
    unit_jumps = np.linalg.norm(np.diff(unit, axis=0), axis=1)
    position = np.argmax(jumps) + 1
    return dict(knots=coordinate.tolist(), natural_positions=line.tolist(),
                piecewise_tangent_vectors=slopes.tolist(),
                maximum_full_tangent_jump=float(np.max(jumps)),
                maximum_unit_tangent_jump=float(np.max(unit_jumps)),
                worst_knot=float(coordinate[position]),
                tangent_from_left=slopes[position - 1].tolist(), tangent_from_right=slopes[position].tolist(),
                minimum_uniform_tangent_error_for_any_C1_extension=float(np.max(jumps) / 2),
                strictly_positive_trace_speed=bool(np.all(speed > 0)))


def run(study, out):
    out.mkdir(exist_ok=False)
    _, t, _, _ = load(study, 'primary', out)
    buffers = t.b.buffers
    e, z = buffers['reg_e'], buffers['reg_t']
    lab = buffers['labels'].reshape(len(z), len(e), 2).copy()
    lab[:, :, 0] += 4
    cut = buffers['x'][(len(buffers['x']) - 1) * 2 // 3]
    if not np.isclose(cut, np.exp(-t.b.data.pressure * .5), atol=1e-15, rtol=0):
        raise ValueError('Unexpected weld label')
    ix = np.flatnonzero(e >= 0)
    if e[ix[0]] != 0 or e[ix[-1]] != .5:
        raise ValueError('Unexpected active collar registry')
    result = dict(reviewer_scope='Native parameterized boundary traces and unit tangent directions. No change to any coordinate or field.',
                  script_sha256=digest(Path(__file__)), weld_chart=float(cut),
                  weld_eta_half=measure(lab[:, -1], z),
                  outer_wall_eta_zero=measure(lab[:, ix[0]], z),
                  entrance_t_zero=measure(lab[0, ix], e[ix]),
                  terminal_t_one=measure(lab[-1, ix], e[ix]))
    (out / 'result.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    run(args.study.resolve(), args.output.resolve())
