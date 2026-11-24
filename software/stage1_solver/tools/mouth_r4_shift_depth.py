"""Held-depth work with independent Q1 depth grids for the labels and crystalline order.

The signed energy rule splits at both grids in both signed shapes. Regions refer
to the fixed CP10 χ front neighborhood on the native spatial chart.
"""
import ctypes as ct
from pathlib import Path
import json
import subprocess

import numpy as np

from mouth_cp10_unfolded_depth import H, LEVELS, DEGREES, KERNEL, control_context, read, write_new
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest
from mouth_free_shape_solid import Data, Pointer, pointer

ARCHIVED_KERNEL = '1bdeba71c1c584cb43046ce222434a503f756f3c38e963949eeeb3318d91d6a0'


def independent_order_kernel():
    if digest(KERNEL)!=ARCHIVED_KERNEL:raise ValueError('Signed-kernel provenance')
    source=KERNEL.read_text()
    def replace(old,new):
        nonlocal source
        if source.count(old)!=1:raise ValueError(old)
        source=source.replace(old,new)
    replace('double L;\n};','double L; Grid order_t;\n};')
    replace('auto t=basis(g.t,p[1]);','auto t=basis(g.t,p[1]); auto ot=basis(bank.order_t,p[1]);')
    replace('double lm=bank.chi[n],lp=bank.chi[n+g.a.n];',
            'auto on=((i+X)*bank.order_t.n+ot.i)*g.a.n+k+A; double lm=bank.chi[on],lp=bank.chi[on+g.a.n];')
    replace('order[X][A]=t.w[0]*lm+t.w[1]*lp;', 'order[X][A]=ot.w[0]*lm+ot.w[1]*lp;')
    replace('for(const Bank &b:banks)append_knots(zz,b.mesh.t,l,u);',
            'for(const Bank &b:banks){append_knots(zz,b.mesh.t,l,u);append_knots(zz,b.order_t,l,u);}')
    replace('double h,weld;Rows rows;', 'double h,weld;Rows rows; const unsigned char *seam; long double local[2][2]{};')
    replace('List parts{0.,1.};','const int at=(cell(d.field_x,PM[0])*(d.field_t.n-1)+cell(d.field_t,PM[1]))*(d.field_a.n-1)+cell(d.field_a,PM[2]); const bool front=seam[at]; List parts{0.,1.};')
    replace('rows.force[s]+=w*((energies[0][s]-energies[1][s])/(2*h));','const double df=w*((energies[0][s]-energies[1][s])/(2*h));rows.force[s]+=df;local[front][s]+=df;')
    replace('int first,int past,double *out){','int first,int past,double *out,const double *om,const double *op,const unsigned char *seam){')
    replace('ym},chi,lm},Bank{SpatialLabels{{nx,x},{nt,tp},{na,a},yp},chi,lp}},h,cut,{}};',
            'ym},chi,lm,{nt,om}},Bank{SpatialLabels{{nx,x},{nt,tp},{na,a},yp},chi,lp,{nt,op}}},h,cut,{},seam,{}};')
    replace('double high=std::max(tm[active_end],tp[active_end]);','double high=std::max(om[active_end],op[active_end]);')
    replace('out[c++]=r.min_j;','out[c++]=r.min_j;for(auto &reg:job.local)for(auto force:reg)out[c++]=double(force);')
    return source


def sector(root,area,k):
    ctx,solid=control_context(root,root/'r4_bucket_unfolded01')
    body=ctx[2]; body.data.offset=pointer(solid)
    if len(body.a)!=7 or k not in range(6):raise ValueError('Expected six angular intervals')
    work=area/f'sector_{k}'
    work.mkdir(exist_ok=False)
    cpp=work/'kernel.cpp';cpp.write_text(independent_order_kernel())
    lib=work/'kernel.so'
    subprocess.run(['g++','-std=c++17','-O3','-shared','-fPIC','-I',str(KERNEL.parent),str(cpp),'-o',str(lib)],check=True)
    fn=ct.CDLL(str(lib)).matched_depth;fn.restype=ct.c_int
    byte=ct.POINTER(ct.c_ubyte)
    fn.argtypes=([ct.POINTER(Data),Pointer,ct.c_int,Pointer,ct.c_int,ct.c_int]+[Pointer]*5+
                 [ct.c_double]*3+[Pointer]*2+[ct.c_int]*3+[Pointer,Pointer,Pointer,byte])
    masks=root/'r4_bucket_geometry01/masks.npz'
    record=dict(parent_sha256=CHECKPOINT_SHA,control_sha256=digest(root/'r4_bucket_unfolded01/control.npz'),
                masks_sha256=digest(masks),generated_kernel_sha256=digest(cpp),sector=k,rows={})
    if read(root/'r4_bucket_geometry01/geometry.json')['masks_sha256']!=digest(masks):raise ValueError('Unexpected front masks')
    with np.load(masks) as data:
        seam=np.ascontiguousarray(data['seam'],dtype=np.uint8)
    from mouth_cp10_spatial_crystal import signed_geometry
    lengths=[signed_geometry(ctx,0,s*H).data.length for s in (-1,1)]
    for level in LEVELS:
        for kind in ('unshifted','shifted'):
            path=area/f'labels_{level}'/(kind+'.npz')
            if digest(path)!=read(path.parent/'result.json')['cases'][kind]['sha256']:raise ValueError('Bad label source')
            with np.load(path) as bank:
                x,a,chi=[bank[c] for c in ('x','a','order')]
                ym,yp=[bank[c+'_labels'] for c in ('minus','plus')]
                om,op=[bank[c+'_order_t'] for c in ('minus','plus')]
                tm,tp=[bank[c+'_label_t'] for c in ('minus','plus')]
            for degree in DEGREES:
                q,w=np.polynomial.legendre.leggauss(degree)
                q,w=np.ascontiguousarray((1+q)*.5),np.ascontiguousarray(w*.5)
                result=np.zeros(14)
                rc=fn(ct.byref(body.data),pointer(x),len(x),pointer(a),len(a),len(tm),
                      *(pointer(ar) for ar in (tm,tp,ym,yp,chi)),*lengths,H,pointer(q),pointer(w),degree,k,k+1,
                      pointer(result),pointer(om),pointer(op),seam.ctypes.data_as(byte))
                if rc or not np.isfinite(result).all() or result[9]<=0 or result[7]<=0:raise ArithmeticError(str(result))
                f=result[4:6];local=result[10:].reshape(2,2)
                if max(abs(f-local.sum(axis=0)))>1e-12 or max(abs(f-(result[:2]-result[2:4])/(2*H)))>1e-8:
                    raise ArithmeticError('Signed force or local sum mismatch')
                record['rows'][f'{kind}_{level}_{degree}']=dict(energies=result[:4].reshape(2,2).tolist(),force=f.tolist(),
                    local=local.tolist(),min_Jinv=float(result[9]),points=int(result[7]))
    write_new(work/'result.json',record)
    print(json.dumps(record,allow_nan=False),flush=True)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path);p.add_argument('area',type=Path);p.add_argument('sector',type=int)
    a=p.parse_args();sector(a.root.resolve(),a.area.resolve(),a.sector)
