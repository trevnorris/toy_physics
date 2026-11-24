"""Compare the eight trilinear order coefficients in the R4 blocking cell across one rematch."""
import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np

from mouth_cp10_label_grid import get_cp10
from mouth_exact_checkpoint_control import digest, read


def compare(scratch, corner_file):
    held = scratch / 'r4_barrier_held01'
    stop = read(held / 'primary/held_stop.json')
    corner = read(corner_file)
    raw = read(scratch / 'r4_barrier_work01/PRECISION_GATE.json')
    prior = scratch / 'r4_eulerian_stage2_01/primary/label_stall.npz'
    final = held / 'primary' / stop['file']
    if (digest(prior) != raw['cases']['primary']['source_sha256']
            or digest(final) != stop['sha256'] or corner['accepted_sha256'] != stop['sha256']):
        raise ValueError('Order profiles must address the same admitted states')
    body, checkpoint, _ = get_cp10(scratch)
    cell = np.array(corner['cell'])
    nodes = []
    with np.load(prior) as old, np.load(final) as new:
        if not np.array_equal(old['x'], new['x']) or not np.array_equal(new['x'], checkpoint['x']):
            raise ValueError('Physical order grid changed')
        before, after = [s['field'][2, body.ids].transpose(0,2,1) for s in (old, new)]
        for offset in product((0,1), repeat=3):
            v = cell + offset
            nodes.append(dict(local=offset, node=v.tolist(), before=float(before[tuple(v)]),
                              after=float(after[tuple(v)])))
    return dict(held_before_sha256=digest(prior), accepted_sha256=digest(final),
                corner_source_sha256=digest(corner_file), cell=cell.tolist(), vertices=nodes,
                before_positive_vertices=sum(v['before']>0 for v in nodes),
                after_positive_vertices=sum(v['after']>0 for v in nodes))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('scratch', type=Path)
    p.add_argument('corner', type=Path)
    p.add_argument('output', type=Path)
    args = p.parse_args()
    result = compare(args.scratch, args.corner)
    with args.output.open('x') as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write('\n')
    print(result)
