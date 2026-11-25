"""Independent one-sided energy and material-triangle readback at the first axial crossing."""
import argparse
import ctypes as ct
import json
from pathlib import Path

import numpy as np

from mouth_r4_aligned_derivatives import load
from mouth_r4_depth_cell_geometry import NAMES, compile_trace, digest
from mouth_r4_full_label_precision import P
from mouth_r4_smooth_label_derivatives import CrystalBanks


def inspect(area, output):
    output.mkdir(exist_ok=False)
    cross = json.loads((area / 'crossings/result.json').read_text())
    raw = area / 'directions/directions.npz'
    if digest(raw) != cross['directions_sha256']:
        raise ValueError('Direction digest changed')
    first = cross['first']
    cells = np.asarray(first['cells'], dtype=np.int32)
    if first['knot'] != 'endpoint' or len(cells) != 2:
        raise ValueError('Expected a shared endpoint in the saved production bank')
    with np.load(raw, allow_pickle=False) as saved:
        labels, direction, order = [saved[key].astype(np.longdouble) for key in ('labels', 'pj', 'order')]
    _, trial, _, _ = load(area.parent / 'r4_rematch39_01', 'primary', output)
    bank = CrystalBanks(trial.b, order, output)
    sample = compile_trace(output)
    alpha = np.longdouble(first['alpha'])
    fractions = [np.longdouble(i) * np.longdouble(s) for s in ('0.0001', '0.00001') for i in (-2, -1, 1, 2)]
    physical_work = []
    result = []
    for cell in cells:
        values = []
        branches = []
        for f in fractions:
            state = np.ascontiguousarray(labels + (1 + f) * alpha * direction)
            rows = np.empty((8192, len(NAMES)), dtype=np.longdouble)
            center = np.empty(3, dtype=np.longdouble)
            n = sample(ct.byref(bank.data), state.ctypes.data_as(P), bank.order.ctypes.data_as(P),
                       bank.rule, bank.weights.ctypes.data_as(P), cell.ctypes.data,
                       rows.ctypes.data_as(P), len(rows), center.ctypes.data_as(P))
            if n <= 0:
                raise ArithmeticError(('Production physical quadrature failed', cell.tolist(), n))
            rows = rows[:n]
            ray = rows[(rows[:, 0] == first['ray'][0]) & (rows[:, 1] == first['ray'][1])]
            # The full, finite quadrature weight is already included in column 15.
            values.append(ray[:, 15].sum(dtype=np.longdouble))
            branches.append(ray[::3, [20, 21, 22]].astype(int).tolist())
        slopes = []
        for i in (0, 4):
            h = alpha * (fractions[i + 1] - fractions[i])
            slopes.append(dict(halfwidth=float(h), before=float((values[i + 1] - values[i]) / h),
                               after=float((values[i + 3] - values[i + 2]) / h)))
        result.append(dict(cell=cell.tolist(), bank=slopes, reference_triangles_before=branches[5],
                           reference_triangles_after=branches[6]))
        physical_work.append(np.asarray(values))
    summed = np.sum(physical_work, axis=0, dtype=np.longdouble)
    pair = []
    for i in (0, 4):
        h = alpha * (fractions[i + 1] - fractions[i])
        left, right = (summed[i + 1] - summed[i]) / h, (summed[i + 3] - summed[i + 2]) / h
        pair.append(dict(halfwidth=float(h), before=float(left), after=float(right), jump=float(right - left)))
    data = dict(state_sha256=cross['state_sha256'], crossing_sha256=digest(area / 'crossings/result.json'),
                script_sha256=digest(Path(__file__)), direction_sha256=digest(raw), ray=first['ray'],
                alpha=float(alpha), relative_offsets=list(map(float, fractions)), triangle_keys=['eta', 't', 'upper'],
                measured='Physical energy of the two incident quadrature rays only, sampled through depth_cell_trace.',
                cells=result, pair=pair)
    (output / 'result.json').write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('area', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    inspect(args.area.resolve(), args.output.resolve())
