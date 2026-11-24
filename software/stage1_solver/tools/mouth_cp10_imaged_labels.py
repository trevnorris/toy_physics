"""Read the checkpoint-10 common label and compare old and projected crystal energies
at identical physical locations. Re-evaluate the eight fixed-field current readings.
"""
import argparse
import ctypes as ct
import json
import subprocess
from pathlib import Path

import numpy as np

from mouth_exact_checkpoint_control import CHECKPOINT_SHA, PORTS, acquired, digest


def compare(root, labels, folder):
    folder.mkdir(exist_ok=False, parents=True)
    build, state, q, u, _, saved = acquired(root/'r4_melted_survival01', root/'r4_exact_support01/patch')
    trial = build(state)
    body = trial.b
    # Use the pinned body's ABI, obtained after acquiring checkpoint 10.
    from mouth_free_shape_solid import Data, Results, Pointer, pointer
    library = folder/'imaged_labels.so'
    subprocess.run(['g++', '-O3', '-std=c++17', '-shared', '-fPIC', str(Path(__file__).with_suffix('.cpp')), '-o', str(library)], check=True)
    kernel = ct.CDLL(str(library)).compare_images
    kernel.argtypes = [ct.POINTER(Data), ct.POINTER(Results)]+[Pointer]*4
    kernel.restype = ct.c_int
    with np.load(labels, allow_pickle=False) as f:
        if str(f['source_sha256']) != CHECKPOINT_SHA or not np.array_equal(f['order'],u[2,body.ids].transpose(0,2,1)):
            raise ValueError('The physical order or label source differs from checkpoint 10')
        xi = np.ascontiguousarray(f['labels'])
    cells = np.zeros((len(body.t)-1,len(body.x)-1,len(body.a)-1,2,2))
    groups = np.zeros((2,3))
    load = np.zeros(body.ids.shape)
    body.work = lambda data, results: kernel(data, results, pointer(xi), pointer(cells), pointer(groups), pointer(load))
    native = body.evaluate(q, u[2], load=False, gradient=False)[0]
    projected_load = np.bincount(body.ids.ravel(), weights=load.ravel(), minlength=trial.model.size)
    flow = trial.measure(u, projected_load+trial.ext)
    ports = {k:float(flow[k]) for k in PORTS}
    delta = {k:ports[k]-saved['current'][k] for k in PORTS}
    sums = np.sum(cells.reshape(-1,2),axis=0,dtype=np.longdouble)
    diff = cells[...,1]-cells[...,0]
    largest = np.unravel_index(abs(diff).argmax(),diff.shape)
    old_held = saved['energy']
    if abs(sums[0]-native)>1e-12 or groups[:,0].sum() == 0:
        raise ArithmeticError('Native quadrature control failed')
    out = dict(checkpoint=CHECKPOINT_SHA, label_hash=digest(labels), native=float(sums[0]), projected=float(sums[1]),
               energy_difference=float(sums[1]-sums[0]), total_at_same_points=float(old_held+sums[1]-sums[0]),
               maximum_element_delta=float(np.max(abs(diff))), most_changed_element=list(map(int,largest)),
               most_changed_energies=cells[largest].tolist(), port_values=ports, port_deltas=delta,
               categories=['same reference side', 'reassigned by reconstructed bond'],
               columns=['active_points','native_energy','projected_energy'], by_side_correspondence=groups.tolist())
    (folder/'result.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('root','labels','output'):p.add_argument(key,type=Path)
    args=p.parse_args()
    compare(args.root.resolve(),args.labels.resolve(),args.output.resolve())
