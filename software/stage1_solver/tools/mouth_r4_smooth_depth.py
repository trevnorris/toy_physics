"""Smooth analytic inverse-chart depth work, its signed energy check and spatial bins."""
import argparse
import ctypes as ct
import json
from pathlib import Path
import subprocess

import numpy as np

from mouth_cp10_spatial_crystal import load_context,signed_geometry
from mouth_exact_checkpoint_control import CHECKPOINT_SHA,digest,read
from mouth_free_shape_solid import Data,Grid,Pointer,pointer


def run(root,output,degree):
    output.mkdir(parents=True,exist_ok=False)
    cp=root/'r4_melted_survival01/coarse/state_00010.npz'
    if digest(cp)!=CHECKPOINT_SHA: raise ValueError('Wrong checkpoint')
    context=load_context(root,root/'r4_cp10_labels01/labels.npz')
    body=context[2]
    mask_path=root/'r4_clearance_regions01/masks.npz'
    if digest(mask_path)!='9accf8221d8839d4aa81aa8c818903fe7bb1fc590246538cc22223f9b6b8a254':
        raise ValueError('Unexpected smooth weld and front partition')
    with np.load(mask_path) as saved:
        tags=np.ascontiguousarray(saved['bucket'],dtype=np.uint8)
        for ax in 'xta':
            if not np.array_equal(saved[ax],body.buffers['field_'+ax]):
                raise ValueError('Wrong chart for masks')
    order=np.ascontiguousarray(context[-1])
    if tags.shape!=tuple(n-1 for n in order.shape) or np.any(tags>2): raise ValueError('Wrong mask size')
    source=Path(__file__).with_suffix('.cpp');shared=source.with_name('mouth_r4_smooth_points.hpp');lib=output/'kernel.so'
    subprocess.run(['g++','-O3','-std=c++17','-shared','-fPIC',str(source),'-o',str(lib)],check=True)
    kernel=ct.CDLL(str(lib));analytic=kernel.smooth_depth_jet;fd=kernel.smooth_depth_fd
    byte=ct.POINTER(ct.c_ubyte); countp=ct.POINTER(ct.c_longlong)
    analytic.argtypes=[ct.POINTER(Data),ct.c_int,Pointer,byte,Grid,Pointer,ct.c_int,Pointer,countp];analytic.restype=ct.c_int
    fd.argtypes=[ct.POINTER(Data),Grid,ct.c_int,Pointer,byte,Grid,Pointer,Pointer,countp];fd.restype=ct.c_int
    gp,w=np.polynomial.legendre.leggauss(degree); gp,w=np.ascontiguousarray((1+gp)/2),np.ascontiguousarray(w/2)
    rule=Grid(len(gp),pointer(gp)); mask=tags.ctypes.data_as(byte)
    common=read(source.parents[3]/'research/pde_ledger_v3/_measurements/mouth_r4_matched_knot.json')['common_bath_depth_force']
    half_steps=(4e-7,2e-7)
    signed_shapes={s*h:signed_geometry(context,0,s*h) for h in half_steps for s in (-1,1)}
    record=dict(checkpoint_sha256=CHECKPOINT_SHA,mask_sha256=digest(mask_path),cpp_sha256=digest(source),
                partition_sha256=digest(shared),degree=degree,buckets=['A','B','C'],banks=['analytic source','spatial Q1'],
                columns=['stored crystal energy','full log-depth crystal force','order work','remaining depth work'],
                common_bath_depth_force=common,rows={})
    for level in (1,2,4):
        rows=[];counts=[]
        for s in (-1,1):
            arr=np.zeros((3,2,3)); count=ct.c_longlong()
            code=analytic(ct.byref(body.data),level,pointer(order),mask,rule,pointer(w),s,pointer(arr),ct.byref(count))
            if code or count.value<=0 or not np.isfinite(arr).all():raise ArithmeticError(f'Failed analytic map ({code})')
            rows.append(arr);counts.append(count.value)
        static=np.mean(rows,axis=0)
        with_rem=np.concatenate((static,(static[:,:,1]-static[:,:,2])[:,:,None]),axis=-1)
        entry=dict(analytic_by_region=with_rem.tolist(),analytic_counts=counts,
                   directional_spread=np.abs(rows[0]-rows[1]).tolist(),
                   error_by_region=(with_rem[:,1]-with_rem[:,0]).tolist(),
                   full_force=(common+np.sum(static[:,:,1],axis=0)).tolist(),fd=[])
        for h in half_steps:
            energy=[]
            for sign in (-1,1):
                ct_body=signed_shapes[sign*h]
                arr=np.zeros((3,2));count=ct.c_longlong()
                code=fd(ct.byref(ct_body.data),body.data.field_t,level,pointer(order),mask,rule,pointer(w),pointer(arr),ct.byref(count))
                if code or count.value<=0 or not np.isfinite(arr).all():raise ArithmeticError(f'Failed signed map ({code})')
                energy.append(arr)
            force=(energy[0]-energy[1])/(2*h)
            entry['fd'].append(dict(h=h,signed_energies=np.asarray(energy).tolist(),force=force.tolist(),
                analytic_minus_fd_by_region=(static[:,:,1]-force).tolist(),
                analytic_minus_fd_total=(static[:,:,1]-force).sum(axis=0).tolist()))
        record['rows'][str(level)]=entry
        with (output/f'level_{level}.json').open('x') as stream:json.dump(entry,stream,indent=2,allow_nan=False)
    with (output/'result.json').open('x') as stream:
        json.dump(record,stream,indent=2,allow_nan=False);stream.write('\n')
    for level,value in record['rows'].items():
        print(level,'error',np.sum(value['error_by_region'],axis=0).tolist(),'analytic-FD',value['fd'][-1]['analytic_minus_fd_total'])

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path);p.add_argument('output',type=Path);p.add_argument('degree',type=int,choices=[6,8]);a=p.parse_args()
    run(a.root.resolve(),a.output.resolve(),a.degree)
