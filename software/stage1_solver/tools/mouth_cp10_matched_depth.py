"""Common moving-knot CP10 depth-force comparison on the saved two-metric crystal."""
import argparse
import ctypes as ct
import json
from pathlib import Path
import subprocess

import numpy as np

from mouth_cp10_nested_forces import subdivide
from mouth_cp10_spatial_crystal import load_context, signed_geometry
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest
from mouth_free_shape_solid import Data, Pointer, pointer


def run(area, out, level, degree):
    out.mkdir(parents=True, exist_ok=False)
    src = Path(__file__).with_suffix('.cpp')
    lib = out/'kernel.so'
    subprocess.run(['g++', '-O3', '-fPIC', '-shared', '-std=c++17', str(src), '-o', str(lib)], check=True)
    native = area/'r4_cp10_labels01/labels.npz'
    lab = native if level == 1 else area/f'r4_nested_labels01/{level}/labels.npz'
    ctx = load_context(area, native)
    b = ctx[2]
    with np.load(area/'r4_melted_survival01/coarse/state_00010.npz') as f:
        solid = np.ascontiguousarray(f['solid'])
    if solid.shape != b.chart_shape:
        raise ValueError('Crystal grid disagreement')
    b.data.offset = pointer(solid)
    with np.load(lab) as f:
        if str(f['source_sha256']) != CHECKPOINT_SHA:
            raise ValueError('Checkpoint identity')
        ax = [np.ascontiguousarray(f[k]) for k in 'xta']
        xi = np.ascontiguousarray(f['labels'])
        chi = np.ascontiguousarray(f['order'])
        mask = f['active_stencil']
        for c, a in zip('xta', ax):
            if np.max(abs(subdivide(b.buffers['field_'+c], level)-a))>5e-16:
                raise ValueError('Wrong knot motion')
        if level>1:
            minus,plus = [np.ascontiguousarray(f['dxi_dt_'+s]) for s in ('minus','plus')]
        else:
            with np.load(area/'r4_cp10_labels01/label_transport.npz') as transport:
                if str(transport['label_sha256']) != digest(lab):
                    raise ValueError('Signed tangent provenance')
                minus,plus = [np.ascontiguousarray(transport['dxi_dt'][...,i,:]) for i in (0,1)]
    fun = ct.CDLL(str(lib)).matched_depth
    fun.restype=ct.c_int
    fun.argtypes=([ct.POINTER(Data),Pointer,ct.c_int,Pointer,ct.c_int,ct.c_int]+[Pointer]*5+
                  [ct.c_double]*3+[Pointer]*2+[ct.c_int]*3+[Pointer])
    gp,gw=np.polynomial.legendre.leggauss(degree)
    gp=np.ascontiguousarray((gp+1)/2);gw=np.ascontiguousarray(gw/2)
    record=dict(checkpoint_sha256=CHECKPOINT_SHA,labels_sha256=digest(lab),cxx_sha256=digest(src),factor=level,
                degree=degree,coordinate_order='native polar angle / normalized physical depth / horizontal affine triangle section',steps=[])
    def save():
        target=out/'partial.json'
        target.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
        target.replace(out/'readout.json')
    for h in (4.e-7,2.e-7):
        ts=[];ys=[];ls=[]
        for sign in (-1,1):
            body=signed_geometry(ctx,0,sign*h)
            st=body.buffers['field_t']
            ts.append(np.ascontiguousarray(subdivide(st,level)))
            delta=subdivide(st-b.buffers['field_t'],level)[None,:,None,None]
            slope=np.where(delta>=0,plus,minus)
            inc=np.multiply(delta,slope,out=np.zeros_like(xi),where=delta!=0)
            mapped=np.ascontiguousarray(xi+inc)
            if not np.all(np.isfinite(mapped[mask])):raise ValueError('Nonfinite ordered nodal map')
            ys.append(mapped)
            ls.append(body.data.length)
        vals=np.zeros(10)
        rc=fun(ct.byref(b.data),pointer(ax[0]),len(ax[0]),pointer(ax[2]),len(ax[2]),len(ax[1]),
               pointer(ts[0]),pointer(ts[1]),pointer(ys[0]),pointer(ys[1]),pointer(chi),*ls,h,
               pointer(gp),pointer(gw),degree,0,len(b.a)-1,pointer(vals))
        if rc or not np.all(np.isfinite(vals)) or vals[9]<=0 or vals[7]<=0:
            raise ArithmeticError(f'Invalid matched material support: {rc}, {vals.tolist()}')
        em,ep=vals[:2],vals[2:4];forces=vals[4:6]; energy_fd=(em-ep)/(2*h)
        if np.max(abs(energy_fd-forces))>2e-8:
            raise ArithmeticError(f'Energy/paired-force accumulation disagreement {energy_fd-forces}')
        record['steps'].append(dict(h=h,energy_minus=em.tolist(),energy_plus=ep.tolist(),
            forces_native_projected=forces.tolist(),projected_minus_native=float(forces[1]-forces[0]),
            summed_energy_fd_minus_paired_force=(energy_fd-forces).tolist(),
            material_chart_volume_through_global_order_bound=vals[6],shared_points=int(vals[7]),
            signed_positive_order_evaluations=int(vals[8]),min_Jinv_over_positive_samples=vals[9]))
        save()
    record['completed']=True
    save()


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('area',type=Path);p.add_argument('out',type=Path)
    p.add_argument('level',type=int,choices=(1,2,4));p.add_argument('degree',type=int,choices=(6,8))
    a=p.parse_args();run(a.area.resolve(),a.out.resolve(),a.level,a.degree)
