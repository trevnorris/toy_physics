"""Frozen R4 smooth-map control for energy and the five modes with stationary order knots."""
import argparse
import ctypes as ct
import json
from pathlib import Path
import subprocess

import numpy as np

from mouth_cp10_spatial_crystal import MODES, load_context, signed_geometry
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
from mouth_free_shape_solid import Data, Grid, Pointer, pointer


def run(root, out, degree):
    out.mkdir(parents=True, exist_ok=False)
    src=Path(__file__).with_suffix('.cpp')
    lib=out/'forces.so'
    subprocess.run(['g++','-std=c++17','-O3','-shared','-fPIC',str(src),'-o',str(lib)],check=True)
    fn=ct.CDLL(str(lib)).clearance_five
    fn.argtypes=[ct.POINTER(Data),ct.c_int,Pointer,ct.POINTER(ct.c_ubyte),Grid,Pointer,
                 ct.POINTER(ct.c_longdouble),Pointer,ct.POINTER(ct.c_longlong)]
    fn.restype=ct.c_int
    ctx=load_context(root,root/'r4_cp10_labels01/labels.npz')
    base,chi=ctx[2],ctx[-1]
    mask_dir=root/'r4_clearance_regions01'
    with np.load(mask_dir/'masks.npz') as m:
        mask=np.ascontiguousarray(m['bucket'],dtype=np.uint8)
    old_mask=root/'r4_bucket_geometry01/masks.npz'
    with np.load(old_mask) as m:
        previous=m['bucket']
    if mask.shape!=tuple(s-1 for s in chi.shape) or set(np.unique(mask))!={0,1,2}:
        raise ValueError('Unexpected region mask')
    hs=(1e-4,5e-5)
    signed=[signed_geometry(ctx,m,sign*h) for h in hs for m in range(1,6) for sign in (-1,1)]
    for body in signed:
        if not np.array_equal(body.buffers['field_t'],base.buffers['field_t']):
            raise ValueError('Unexpected depth grid motion at fixed axial length')
    bank=(Data*21)(base.data,*(b.data for b in signed))
    p,w=np.polynomial.legendre.leggauss(degree)
    p=np.ascontiguousarray((p+1)/2);w=np.ascontiguousarray(w/2)
    common_path=src.parents[3]/'research/pde_ledger_v3/_measurements/mouth_r4_two_metric_labels.json'
    field=read(common_path)['data']['held_fields']
    if field['checkpoint']!=CHECKPOINT_SHA:raise ValueError('Wrong bath state')
    bath=np.array([rec['analytic'] for rec in field['modes'][1:]])
    bath_fd=np.array([[rec['central'][str(h)]['force'] for rec in field['modes'][1:]] for h in hs])
    record=dict(checkpoint=CHECKPOINT_SHA,common_field_hash=digest(common_path),source=digest(src),
                common_quadrature=digest(src.with_name('mouth_r4_smooth_points.hpp')),
                mask_sha256=digest(mask_dir/'masks.npz'),old_mask_sha256=digest(old_mask),
                quadrature_degree=degree,force_modes=MODES[1:],signed_steps=hs,
                columns=['stored_crystal_energy','order_work',*MODES[1:]],results={})
    for level in (1,2,4):
        values=np.zeros((*mask.shape,2,27),dtype=np.longdouble)
        J=np.zeros(2); count=ct.c_longlong()
        rc=fn(bank,level,pointer(chi),mask.ctypes.data_as(ct.POINTER(ct.c_ubyte)),Grid(degree,pointer(p)),
              pointer(w),values.ctypes.data_as(ct.POINTER(ct.c_longdouble)),pointer(J),ct.byref(count))
        if rc or not np.isfinite(values).all() or J.min()<=0 or count.value<=0:
            raise ArithmeticError('Invalid two-metric energy or force')
        sums=values.sum(axis=(0,1,2),dtype=np.longdouble)
        energy=sums[:,7:].reshape(2,2,5,2)
        signed_fd=(energy[:,:,:,0]-energy[:,:,:,1])/(2*np.asarray(hs)[None,:,None])
        analytic=sums[:,2:7]
        total=analytic+bath
        total_fd=signed_fd+bath_fd
        rich=(4*total_fd[:,1]-total_fd[:,0])/3
        tails=values[...,1,:7]-values[...,0,:7]
        def group(tags):
            return np.stack([tails[tags==k].sum(axis=0,dtype=np.longdouble) for k in range(3)])
        by_region=group(mask)
        delta=sums[1,:7]-sums[0,:7]
        close=float(np.max(abs(by_region.sum(axis=0)-delta)))
        if close>3e-14:raise ArithmeticError('Region accounting did not close')
        def rows(x):return np.asarray(x,dtype=float).tolist()
        record['results'][level]=dict(material_spatial=rows(sums[:,:7]),delta=rows(delta),
            native_region_delta=rows(group(previous)),synthetic_region_delta=rows(by_region),regional_sum_delta=close,
            analytic_total=rows(total),central_total=rows(total_fd),richardson_total=rows(rich),
            total_analytic_minus_richardson=rows(total-rich),
            crystal_analytic_minus_central=rows(analytic[:,None,:]-signed_fd),
            signed_crystal_energies=rows(energy),sampled_min_physical_J=rows(J),points=count.value)
        temp=out/'writing.json'
        temp.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
        temp.replace(out/'result.json')
        print(json.dumps(dict(level=level,delta=rows(delta),closure=float(np.max(abs(total-rich))))),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('degree',choices=(6,8),type=int)
    args=parser.parse_args()
    run(args.root.resolve(),args.output.resolve(),args.degree)
