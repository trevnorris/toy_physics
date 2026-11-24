"""Additional fixed-normalized-chart depth virtual work at the original CP10 material points."""
import argparse
import ctypes as ct
import json
from pathlib import Path
import subprocess

import numpy as np

from mouth_cp10_nested_energy import acquire_labels
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, acquired, digest


def calculate(root,rows_path,level,out):
    build,x,q,u,_,_ = acquired(root/'r4_melted_survival01',root/'r4_exact_support01/patch')
    body=build(x).b
    from mouth_free_shape_solid import Data, Results, Pointer, pointer
    axes,xi=acquire_labels(rows_path,body,level)
    t=body.buffers['field_t'];L=body.data.length
    rate=np.where(t*L<L-.12,(t*L/(L-.12))**25,1.)-t
    out.mkdir(exist_ok=False,parents=True)
    lib=out/'depth.so'
    subprocess.run(['g++','-O3','-std=c++17','-shared','-fPIC',str(Path(__file__).with_suffix('.cpp')),'-o',str(lib)],check=True)
    k=ct.CDLL(str(lib)).nested_depth_work
    k.argtypes=[ct.POINTER(Data),ct.POINTER(Results)]+[Pointer,ct.c_int]*3+[Pointer]*3
    k.restype=ct.c_int
    report=np.zeros(7)
    knot_args=[v for z in axes for v in (pointer(z),len(z))]
    body.work=lambda data,results: k(data,results,*knot_args,pointer(xi),pointer(rate),pointer(report))
    e=body.evaluate(q,u[2],load=False,gradient=False)[0]
    if int(report[6])!=5291652 or abs(e-.20652001564580955)>1e-13:raise ValueError('Wrong native quadrature')
    def work(row):return dict(stress=float(row[0]),order_knot_motion=float(row[1]),total=float(sum(row)))
    native=work(report[:2]); projected=work(report[2:4]);gradient=work(report[4:6])
    legacy=json.loads((Path(__file__).resolve().parents[3]/'research/pde_ledger_v3/_measurements/mouth_r4_two_metric_labels.json').read_bytes())
    bath=legacy['data']['total_shape_forces']['field_analytic'][0]
    target=legacy['data']['total_shape_forces']['native_saved'][0]
    data=dict(checkpoint=CHECKPOINT_SHA,label_sha256=digest(rows_path),factor=level,active_points=int(report[6]),
              frozen_inner_crystal=True,kinematics='Hold the label function and its normalized chart gradient at each native image; transport order knots and physical depth.',
              native=native,projected=projected,native_label_spatial_gradient=gradient,
              projected_minus_native={k:projected[k]-native[k] for k in native},analytical_bath=bath,
              saved_full_depth_force=target,native_plus_bath_minus_saved=native['total']+bath-target,
              projected_plus_bath_minus_saved=projected['total']+bath-target)
    (out/'result.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    print(json.dumps(data,allow_nan=False),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path);p.add_argument('labels',type=Path);p.add_argument('factor',type=int,choices=(1,2,4));p.add_argument('output',type=Path)
    a=p.parse_args();calculate(a.root.resolve(),a.labels.resolve(),a.factor,a.output.resolve())
