"""Measure local CP10 matched depth work on the original and fixed-bin quadratures.

Both methods share all smooth integration intervals. Only the very narrow signed
ribbons that cross a frozen native depth knot need a second local integration.
The source augmentation preserves the exact quadrature that produced eda4ffed.
"""
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


def add_bins(source):
    text=source.read_text()
    if digest(source)!='1bdeba71c1c584cb43046ce222434a503f756f3c38e963949eeeb3318d91d6a0':
        raise ValueError('Signed material kernel no longer matches eda4ffed')
    def once(old,new):
        nonlocal text
        if text.count(old)!=1:
            raise ValueError(f'Expected one instrumentation site: {old}')
        text=text.replace(old,new)
    once('namespace matched {','#include "mouth_cp10_bucket_transport.hpp"\nnamespace matched {')
    once('double h,weld;Rows rows;', 'double h,weld;Rows rows; CellWork &bins;')
    once('List parts{0.,1.};','const auto parent=bins.locate(PM); const bool straddles=bins.straddles(PM);\n            List parts{0.,1.};')
    once('}rows.count++;','}rows.count++; bins.deposit(parent,straddles,w,energies,h);')
    once('double z=zz[q-1]+dz*gauss.v[T];Line end[2];','''
                const double lower=zz[q-1], upper=zz[q];
                auto knot=std::upper_bound(d.field_t.v,d.field_t.v+d.field_t.n,lower);
                bool crossing=knot<d.field_t.v+d.field_t.n && *knot<upper;
                if(crossing && knot+1<d.field_t.v+d.field_t.n && knot[1]<upper)
                    throw std::runtime_error("two fixed depth boundaries in one signed slab");
                for(int part=0;part<(crossing?3:1);part++){
                const Rows before=rows;
                bins.partition=!crossing?CellWork::common:part==0?CellWork::original:CellWork::split;
                const double start=part==2?*knot:lower;
                const double width=part==0?dz:part==1?*knot-lower:upper-*knot;
                double z=start+width*gauss.v[T];Line end[2];''')
    once('end[0],end[1],aw*dz*weights[T]);',
         'end[0],end[1],aw*width*weights[T]);\n                if(part>0)rows=before;\n                }')
    once('int first,int past,double *out){','int first,int past,double *out,double *cells,const double *t0,const unsigned char *weld_mask){')
    once('Audit job{','CellWork bins(*d,{nx,x},{nt,t0},{na,a},weld_mask);\n    Audit job{')
    once('yp},chi,lp}},h,cut,{}};','yp},chi,lp}},h,cut,{},bins};')
    once('auto &r=job.rows;', 'bins.export_cells(cells);\n    auto &r=job.rows;')
    return text


def run(area,output,level,degree):
    output.mkdir(parents=True,exist_ok=False)
    tools=Path(__file__).resolve().parent
    old=tools/'mouth_cp10_matched_depth.cpp'
    cpp=output/'kernel.cpp'; cpp.write_text(add_bins(old))
    lib=output/'kernel.so'
    subprocess.run(['g++','-O3','-fPIC','-shared','-std=c++17','-I',str(tools),str(cpp),'-o',str(lib)],check=True)
    base=area/'r4_cp10_labels01/labels.npz'
    ctx=load_context(area,base); body=ctx[2]
    with np.load(area/'r4_melted_survival01/coarse/state_00010.npz') as cp:
        solid=np.ascontiguousarray(cp['solid']);body.data.offset=pointer(solid)
    path=base if level==1 else area/f'r4_nested_labels01/{level}/labels.npz'
    with np.load(path) as fd:
        assert str(fd['source_sha256'])==CHECKPOINT_SHA
        x,t,a=(np.ascontiguousarray(fd[c]) for c in 'xta')
        xi,chi,mask=(np.ascontiguousarray(fd[c]) for c in ('labels','order','active_stencil'))
        if level>1:
            minus,plus=(np.ascontiguousarray(fd['dxi_dt_'+s]) for s in ('minus','plus'))
        else:
            with np.load(base.parent/'label_transport.npz') as tr:
                assert str(tr['label_sha256'])==digest(base)
                minus,plus=(np.ascontiguousarray(tr['dxi_dt'][...,i,:]) for i in (0,1))
    oldrow=json.loads((tools.parents[2]/'research/pde_ledger_v3/_measurements/mouth_r4_matched_knot.json').read_text())['matched'][f'{level}_{degree}']
    assert oldrow['checkpoint_sha256']==CHECKPOINT_SHA and oldrow['labels_sha256']==digest(path)
    prior=oldrow['steps'][-1]; h=prior['h']; assert h==2e-7
    signed_t=[];signed_y=[];lengths=[]
    for sign in (-1,1):
        shifted=signed_geometry(ctx,0,sign*h); st=shifted.buffers['field_t']
        signed_t.append(np.ascontiguousarray(subdivide(st,level)))
        dt=subdivide(st-body.buffers['field_t'],level)[None,:,None,None]
        rate=np.where(dt>=0,plus,minus)
        inc=np.multiply(dt,rate,out=np.zeros_like(xi),where=dt!=0)
        state=np.ascontiguousarray(xi+inc)
        if not np.isfinite(state[mask]).all(): raise ValueError('Invalid active label')
        signed_y.append(state);lengths.append(shifted.data.length)
    rule,weights=np.polynomial.legendre.leggauss(degree)
    rule=np.ascontiguousarray((1+rule)/2);weights=np.ascontiguousarray(weights/2)
    fn=ct.CDLL(str(lib)).matched_depth
    fn.restype=ct.c_int
    fn.argtypes=([ct.POINTER(Data),Pointer,ct.c_int,Pointer,ct.c_int,ct.c_int]+[Pointer]*5+
                 [ct.c_double]*3+[Pointer]*2+[ct.c_int]*3+[Pointer,Pointer,Pointer,ct.POINTER(ct.c_ubyte)])
    dims=[len(body.buffers['field_'+q])-1 for q in 'xta']
    loc=np.zeros((*dims,6));v=np.zeros(10)
    # Q1 extrema lie at vertices. Strictly opposite signs select cells whose
    # checkpoint label interpolation crosses the bond inside the cell.
    corners=np.lib.stride_tricks.sliding_window_view(xi[...,0],(2,2,2))
    cutoff=body.x[(len(body.x)-1)*2//3]
    weld=np.ascontiguousarray((corners.min(axis=(-3,-2,-1))<cutoff)
                              &(corners.max(axis=(-3,-2,-1))>cutoff),dtype=np.uint8)
    rc=fn(ct.byref(body.data),pointer(x),len(x),pointer(a),len(a),len(t),
          *[pointer(ar) for ar in (*signed_t,*signed_y,chi)],*lengths,h,pointer(rule),pointer(weights),degree,
          0,len(body.a)-1,pointer(v),pointer(loc),pointer(t),weld.ctypes.data_as(ct.POINTER(ct.c_ubyte)))
    if rc or not np.isfinite(loc).all() or v[9]<=0: raise ArithmeticError(f'Matched kernel failure {rc}')
    if int(v[7])!=prior['shared_points'] or np.max(abs(v[4:6]-prior['forces_native_projected']))>1e-13:
        raise ArithmeticError(f'Original paired integration changed {v.tolist()}')
    sums=np.sum(loc,axis=(0,1,2),dtype=np.longdouble).astype(float)
    if max(abs(sums[:2]-v[4:6]))>1e-12: raise ArithmeticError(f'Cell sum mismatch {sums}, {v}')
    data=output/'cells.npz'
    np.savez_compressed(data,force_native=loc[...,0],force_spatial=loc[...,1],
                        cut_native=loc[...,2],cut_spatial=loc[...,3],straddling_native=loc[...,4],straddling_spatial=loc[...,5])
    record=dict(checkpoint_sha256=CHECKPOINT_SHA,archive_kernel_sha256=digest(old),
                label_sha256=digest(path),factor=level,degree=degree,h=h,original_pair=v[4:6].tolist(),
                original_minus_archived=(v[4:6]-prior['forces_native_projected']).tolist(),
                sum_native_spatial=sums[:2].tolist(),fixed_depth_split_sum_native_spatial=sums[2:4].tolist(),
                points=int(v[7]),min_sampled_Jinv=v[9],shape=dims,fields_sha256=digest(data))
    with (output/'readout.json').open('x') as io: json.dump(record,io,indent=2,allow_nan=False)
    print(json.dumps(record,allow_nan=False),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('area',type=Path);parser.add_argument('out',type=Path)
    parser.add_argument('level',choices=(1,2,4),type=int);parser.add_argument('degree',choices=(6,8),type=int)
    args=parser.parse_args();run(args.area.resolve(),args.out.resolve(),args.level,args.degree)
