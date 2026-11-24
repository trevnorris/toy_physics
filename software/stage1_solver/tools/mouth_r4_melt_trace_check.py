"""Independent work and orientation on the same fixed trace with two melted-zone gauges."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_cp10_label_grid import get_cp10
from mouth_exact_checkpoint_control import digest, read
from mouth_r4_eulerian_relax_labels import LabelRelaxation, absmax
from mouth_r4_eulerian_solid import SpatialCrystal
from mouth_r4_stage2_outcome import verify_guard


def check(scratch, case):
    folder = scratch/'r4_slaved_melt_barrier01'
    verify_guard(folder/(case+'_guard'), 0)
    r = read(folder/case/'readout.json')
    start = scratch/'r4_eulerian_stage2_01'/case
    frozen = start/'label_stall.npz'
    if digest(frozen) != r['sha256'] or digest(frozen) != read(start/'label_stall.json')['sha256']:
        raise ValueError('Endpoint differs from the registered trace')
    with np.load(frozen) as f:
        eps=float(f['epsilon'])
        if eps!=r['epsilon'] or int(f['radius'])!=4 or int(f['winding'])!=1:
            raise ValueError('Wrong gauge or preparation')
        labels, fields, x = (f[k].copy() for k in ('labels','field','x'))
    body,cp,_=get_cp10(scratch)
    if not np.array_equal(cp['x'],x):
        raise ValueError('Not the prescribed fixed shape')
    order=fields[2,body.ids].transpose(0,2,1).copy(order='C')
    opt=LabelRelaxation(body,labels,order,eps)
    state=folder/case/'fixed_active.npz'
    with np.load(state) as f:
        extended=f['labels']
    solid=SpatialCrystal(body)
    E0,g0,o0,_=solid.evaluate(labels,order)
    E1,g1,o1,_=solid.evaluate(extended,order)
    stats=r['static']
    if (not np.array_equal(labels[opt.stencil], extended[opt.stencil])
            or E1!=E0 or not np.array_equal(g0[opt.stencil], g1[opt.stencil])
            or abs(E1-stats['energy_before'])>1.e-12):
        raise ValueError('Physical trace, energy or displacement forces changed')
    changed=o1-o0
    site=tuple(int(i) for i in np.unravel_index(np.argmax(np.abs(changed)),order.shape))
    comparisons={}
    for test_eps in (1.e-8,1.e-10):
        # At fixed trace, multiplying an isolated strictly-melted objective
        # by either positive constant preserves its stationary point.
        gauge=test_eps*stats['unit_gauge']
        comparisons[str(test_eps)]=dict(gauge_energy=gauge,ratio=gauge/E1,
            max_active_gauge_gradient=test_eps*stats['gauge_force_on_free_active']/eps)
    return dict(source_sha256=r['sha256'], extension_sha256=digest(state),
                fixed_trace_label_change=absmax((extended-labels)[opt.stencil]),
                fixed_trace_elastic_energy_change=E1-E0,
                fixed_trace_elastic_force_max_change=absmax((g1-g0)[opt.stencil]),
                maximum_order_source_change=dict(node=list(site),order=float(order[site]),
                    before=float(o0[site]),after=float(o1[site]),change=float(changed[site])),
                order_source_change_on_positive_labels=absmax(changed[order>0]),
                same_trace_gauge_rescalings=comparisons)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('scratch',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    data={case:check(args.scratch.resolve(),case) for case in ('primary','secondary')}
    with args.output.open('x') as f:
        json.dump(data,f,indent=2,allow_nan=False)
        f.write('\n')
    for name,r in data.items():
        print(name,r['maximum_order_source_change'],flush=True)
