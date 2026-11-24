"""Isolate work from shared order-knot motion at native material quadrature points."""
import argparse
import ctypes as ct
import json
from pathlib import Path
import subprocess

import numpy as np

from mouth_cp10_nested_energy import acquire_labels
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, acquired, digest


def calculate(root, degree, output):
    build, shape, q, fields, _, _ = acquired(root/'r4_melted_survival01', root/'r4_exact_support01/patch')
    from mouth_free_shape_solid import Data, Grid, Pointer, Results, pointer
    solid = build(shape).b
    if np.any(solid.buffers['base']):
        raise ValueError('Expected absolute spatial order')
    output.mkdir(parents=True, exist_ok=False)
    library = output/'depth_work.so'
    subprocess.run(['g++','-O3','-std=c++17','-shared','-fPIC',str(Path(__file__).with_suffix('.cpp')),
                    '-o',str(library)], check=True)
    work = ct.CDLL(str(library)).order_work
    byte_pointer = ct.POINTER(ct.c_ubyte)
    work.argtypes = [ct.POINTER(Data),ct.POINTER(Results),ct.POINTER(Grid),ct.POINTER(Pointer),
                     Pointer,byte_pointer,Pointer,ct.POINTER(ct.c_longlong)]
    work.restype = ct.c_int
    paths = [root/'r4_cp10_labels01/labels.npz',
             *(root/f'r4_nested_labels01/{level}/labels.npz' for level in (2,4))]
    tables = [acquire_labels(path,solid,level) for path,level in zip(paths,(1,2,4))]
    axes = (Grid*9)(*(Grid(len(a),pointer(a)) for grid,_ in tables for a in grid))
    labels = (Pointer*3)(*(pointer(xi) for _,xi in tables))
    t = solid.buffers['field_t']
    length = solid.data.length
    tdot = np.where(t*length<length-.12,(t*length/(length-.12))**25,1.)-t
    chi = fields[2][solid.ids]
    last = np.where(np.any(chi>0,axis=(0,1)))[0][-1]
    # Every parent prism maps into the convex hull of its vertices. Exclude
    # only prisms entirely beyond the support of the last positive order knot.
    corner_t = solid.t[:,None,None]+solid.unpack(q)[...,1]
    min_t = np.minimum.reduce([corner_t[j:j+len(solid.t)-1,i:i+len(solid.x)-1,k:k+len(solid.a)-1]
                               for j in (0,1) for i in (0,1) for k in (0,1)])
    active = np.ascontiguousarray(min_t <= t[min(last+1,len(t)-1)],dtype=np.uint8)
    nodes, weights = np.polynomial.legendre.leggauss(degree)
    nodes, weights = (nodes+1.)/2., weights/2.
    solid.data.gauss = Grid(degree,pointer(nodes))
    solid.data.weights = pointer(weights)
    totals = np.zeros((4,3))
    count = ct.c_longlong()
    solid.work = lambda d,r: work(d,r,axes,labels,pointer(tdot),active.ctypes.data_as(byte_pointer),
                                 pointer(totals),ct.byref(count))
    solid.evaluate(q,fields[2],load=False,gradient=False)
    if count.value<=0 or not np.isfinite(totals).all():
        raise ArithmeticError('Invalid depth-work readout')
    record = dict(checkpoint_sha256=CHECKPOINT_SHA,label_sha256=[digest(p) for p in paths],
                  degree=degree,ordered_points=count.value,
                  rows=['native','1x','2x','4x'],
                  columns=['geometric_stress','order_knot_work','sum_absolute_pointwise_order_work_error'],
                  sums=totals.tolist(),projected_minus_native=(totals[1:,:2]-totals[0,:2]).tolist())
    with (output/'result.json').open('x') as stream:
        json.dump(record,stream,indent=2,allow_nan=False)
    print(json.dumps(record,allow_nan=False),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path)
    parser.add_argument('degree',type=int,choices=(6,8))
    parser.add_argument('output',type=Path)
    opts=parser.parse_args()
    calculate(opts.root.resolve(),opts.degree,opts.output.resolve())
