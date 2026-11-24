"""One-sided order variations at the held R4 corner with the production melt mask."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

from mouth_r4_fixed_label_order import digest, read


def run(area, state, order_result, result):
    pin = area/'pin/software/stage1_solver/tools'
    manifest = read(area/'sources.json')
    if any(digest(pin/name) != sha for name, sha in manifest.items()):
        raise ValueError('Incompatible physical source')
    if digest(state) != read(order_result)['saved_sha256']:
        raise ValueError('Incompatible post-order state')
    sys.path.insert(0, str(pin))
    from mouth_continuous_shape_relax import ShapeTrial
    from mouth_continuous_depth_mesh import axes_for_depth
    from mouth_r4_eulerian_release import Driver, block_axes, melted_cells, order_grid
    from mouth_r4_eulerian_solid import SpatialCrystal
    from mouth_r4_barrier_shape import unit_work

    with np.load(state) as src:
        q, labels, fields = (src[key] for key in ('x', 'labels', 'field'))
        epsilon = float(src['epsilon'])
    driver = Driver(area.parent, state.parent, epsilon)
    shape = driver.shape(q)
    trial = ShapeTrial(driver.base, driver.prep, q, 8,24,6,3,2,precondition=False,
                       axes=axes_for_depth(shape.length,driver.reference_depth))
    order = order_grid(trial, fields)
    mask = melted_cells(order)
    axes = block_axes(trial.b)
    _, _, source, _ = SpatialCrystal(trial.b).evaluate(labels,order,False)
    inner = np.bincount(trial.b.ids.ravel(),weights=source.transpose(0,2,1).ravel(),
                        minlength=trial.model.size)
    trial.model.material_load = trial.ext + inner
    # Smooth physical order work follows the exact production variational derivative.
    points = {}
    for name, node in [('corner',(16,10,1)),('opposite',(15,9,2))]:
        ix = trial.number[trial.b.ids[node[0],node[2],node[1]]]
        dofs = np.flatnonzero(trial.number==ix)
        if not trial.free[ix]:
            raise ValueError('Selected order coordinate is externally fixed')
        samples = []
        _, g = trial.model.evaluate(fields)
        g0 = g[2,dofs].sum()
        for h in [1.e-6,1.e-9,1.e-12,1.e-15]:
            probe = fields.copy()
            probe[2,dofs] += h
            _, g1 = trial.model.evaluate(probe)
            probe[2,dofs] = fields[2,dofs] + h/2
            _, gm = trial.model.evaluate(probe)
            # Nonzero h/2 and h have the same incident-cell partition.
            nextmask = melted_cells(order_grid(trial,probe))
            lost = mask & ~nextmask
            physical = float(h*(g0+4*gm[2,dofs].sum()+g1[2,dofs].sum())/6)
            if np.any(nextmask & ~mask):
                raise ArithmeticError('Positive order increment added fully melted cells')
            delta_gauge = -epsilon*float(unit_work(shape,trial.polar_shift,driver.reference_depth,
                                     axes,labels,lost,np.zeros(6))) if np.any(lost) else 0.
            samples.append(dict(h=h,physical_work=physical,gauge_mask_work=delta_gauge,
                                total_work=physical+delta_gauge,lost_cells=int(lost.sum())))
        points[name]=dict(node=node,order=float(order[node]),physical_slope=float(g0),
                          samples=samples)
    with result.open('x') as f:
        json.dump(dict(state_sha256=digest(state),source_sha256=digest(Path(__file__)),
                       epsilon=epsilon,sites=points),f,indent=2,allow_nan=False)
        f.write('\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('area','state','order_result','result'):
        p.add_argument(name,type=Path)
    args=p.parse_args()
    run(args.area.resolve(),args.state.resolve(),args.order_result.resolve(),args.result.resolve())
