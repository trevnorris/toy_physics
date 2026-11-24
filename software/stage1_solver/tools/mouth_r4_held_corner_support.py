"""Determine which nodal labels control the limiting one-sided R4 corner volume."""
import argparse
from collections import Counter
from itertools import product
import json
from pathlib import Path

import mpmath as mp
import numpy as np

from mouth_cp10_label_grid import get_cp10
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
from mouth_r4_eulerian_corner_audit import point
from mouth_r4_refinement_review import _collar_reference
from mouth_r4_stage2_outcome import verify_guard


def audit(scratch):
    run = scratch / 'r4_barrier_held01'
    case = run / 'primary'
    verify_guard(run / 'primary_guard', 1)
    accepted = read(case / 'last_active.json')
    stop = read(case / 'held_stop.json')
    state = case / stop['file']
    if not digest(state) == stop['sha256'] == accepted['sha256'] == digest(case / accepted['file']):
        raise ValueError('Last accepted material differs from the retained stop')
    probes = stop['stop']['probes']
    floor = [p for p in probes if p['alpha'] == 1e-4]
    if len(floor) != 4 or {p['predictor'] for p in floor} != {'fixed_melt', 'harmonic_increment'}:
        raise ValueError('Missing restarted or retained-direction trials at the requested floor')
    cells = {tuple(p['witness']['cell']) for p in floor}
    local = {tuple(p['witness']['chart_fractions']) for p in floor}
    if len(cells) != 1 or len(local) != 1 or any(p['witness']['chart_det'] >= 0 for p in floor):
        raise ValueError('The four trial determinants do not identify one blocking corner')
    cell = np.array(cells.pop())
    bits = np.array(local.pop())
    if not ((bits == 0) | (bits == 1)).all():
        raise ValueError('The failed trial witness is not at a vertex')
    body, cp, axes = get_cp10(scratch)
    with np.load(state) as a:
        if str(a['parent_sha256']) != CHECKPOINT_SHA or int(a['radius']) != 4 or not np.array_equal(a['x'], cp['x']):
            raise ValueError('Not the admitted fixed R4 geometry')
        labels = a['labels']
        chi = a['field'][2, body.ids].transpose(0, 2, 1)
    vertex = cell + bits.astype(int)
    def adjacent(v):
        for q in product((0, 1), repeat=3):
            base = v - q
            if np.all(base >= 0) and np.all(base + 1 < chi.shape):
                yield base, q
    def maximum(c):
        return float(chi[tuple(slice(i, i+2) for i in c)].max())
    determining = []
    D = np.empty((3, 3))
    vertices = [vertex]
    for axis in range(3):
        neighbor = vertex.copy()
        neighbor[axis] += 1 - int(2*bits[axis])
        width = axes[axis][neighbor[axis]] - axes[axis][vertex[axis]]
        D[:, axis] = (labels[tuple(neighbor)] - labels[tuple(vertex)]) / width
        vertices.append(neighbor)
    for v in vertices:
        incident = [maximum(c) for c, _ in adjacent(v)]
        determining.append(dict(node=v.tolist(), chi=float(chi[tuple(v)]),
                                incident_cells=len(incident), positive_order_cells=sum(m > 0 for m in incident)))
    with mp.workdps(65):
        rows = [point(body, axes, labels, chi, c, q, _collar_reference) for c,q in adjacent(vertex)]
    corner = next(row for row in rows if row['cell'] == cell.tolist())
    if abs(corner['chart_det'] - np.linalg.det(D)) > 2e-14:
        raise ArithmeticError('Incident-edge differences and vertex-interpolant derivatives disagree')
    return dict(accepted_sha256=digest(state), cell=cell.tolist(), vertex=vertex.tolist(),
                determinant_from_four_nodes=float(np.linalg.det(D)), determining_nodes=determining,
                accepted_corner=corner, accepted_incident_cells=rows,
                immutable_under_melt_extension=all(v['positive_order_cells']>0 for v in determining),
                trials_at_floor=floor, trial_reason_counts=dict(Counter(p['reason'] for p in probes)),
                source_sha256=digest(Path(__file__)),
                interpretation='A trilinear cell determinant at a corner depends on that vertex and its three adjacent edge vertices. A melt-only extension cannot change any vertex incident to an element with positive order.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('scratch', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    result = audit(args.scratch)
    with args.output.open('x') as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write('\n')
    print(result['cell'], result['vertex'], result['determining_nodes'], result['accepted_corner'], flush=True)
