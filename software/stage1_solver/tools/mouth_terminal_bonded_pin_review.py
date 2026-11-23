"""Check compatibility of a physical crystal pin and the imposed R4 wall motion."""

import argparse
import json
from pathlib import Path

import numpy as np

from mouth_last_ring_solid import source
from mouth_supercore_joined_front import file_sha


def geometry(state, wall, reference):
    shape, displacement, _, _ = source(state, wall)
    measured = json.loads(reference.read_text())['data']
    records = []
    for grid in ('double', 'local'):
        data = measured[grid]
        if data['source']['wall'] != file_sha(wall):
            raise ValueError('Signed wall work uses a different physical geometry')
        rows = [r for r in data['results'] if r['strip'] == 'near' and r['quadrature'] == 2]
        knots = np.array(rows[0]['depth_knots'])
        if len(rows) != 4 or not all(np.array_equal(r['depth_knots'], knots) for r in rows):
            raise ValueError('Different physical bank motion inside a control pair')
        records.append(knots)
    if not np.array_equal(*records):
        raise ValueError('Doubled and locally refined grids use different walls')
    lo, mid, hi = records[0]

    def shift(z):
        return np.maximum(0., np.minimum((z-lo)/(mid-lo), (hi-z)/(hi-mid)))

    # The crystal outer row has both degrees of freedom held to the wall.
    # Read physical marker positions without interpolating a carrier field.
    outer = displacement[:, -1]
    if np.any(outer != 0.):
        raise ValueError('Physical side displacement is already nonzero')
    markers = np.linspace(0., shape.length, len(outer))
    g = shift(markers)
    ids = np.flatnonzero(g > 0)
    if not len(ids):
        raise ValueError('No solid marker in this strip')

    r, slope = shape.radius(np.array([mid]))
    normal_speed = -slope[0] / np.hypot(1., slope[0])
    # Coordinates of the displaced geometric bank in the OLD physical crystal.
    # An axial step eps keeps its radius and changes z by eps at the peak.
    normal_mismatch = []
    for eps in (-1e-5, 1e-5, -5e-6, 5e-6):
        rr, _ = shape.radius(np.array([mid+eps]))
        signed_offset = (r[0]-rr[0]) / np.hypot(1., slope[0])
        normal_mismatch.append([eps, float(signed_offset / eps)])
    return dict(
        state_sha256=file_sha(state), wall_sha256=file_sha(wall),
        bank_work_sha256=file_sha(reference), depth_knots=records[0].tolist(),
        fixed_ring_ids=ids.tolist(), imposed_marker_shift_per_unit=g[ids].tolist(),
        old_marker_motion_max=float(np.abs(outer).max()),
        both_wall_normals=[float(normal_speed), float(-normal_speed)],
        signed_normal_offset_per_step=normal_mismatch,
        inside_enclosure_tent_at_tip=float(shift(shape.length)),
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('state', 'wall', 'reference'):
        parser.add_argument(key, type=Path)
    args = parser.parse_args()
    print(json.dumps(geometry(args.state,args.wall,args.reference), allow_nan=False))
