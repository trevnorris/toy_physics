"""Frozen control energy and order-motion force on a single identical material point set."""
import argparse
import ctypes as ct
from pathlib import Path
import subprocess

import numpy as np

from mouth_cp10_unfolded_depth import LEVELS, control_context, read, write_new
from mouth_exact_checkpoint_control import CHECKPOINT_SHA,digest
from mouth_free_shape_solid import Data,Grid,Pointer,pointer


def run(root,area,degree):
    prev=root/'r4_bucket_unfolded01'; ctx,solid=control_context(root,prev);body=ctx[2]
    with np.load(prev/'control.npz') as cp:full_order=cp['field'][2]
    q=body.pack(solid)
    case=area/f'static_{degree}';case.mkdir(exist_ok=False)
    source=Path(__file__).with_suffix('.cpp');lib=case/'frozen.so'
    subprocess.run(['g++','-O3','-shared','-fPIC','-std=c++17',str(source),'-o',str(lib)],check=True)
    fn=ct.CDLL(str(lib)).frozen_shift;byte=ct.POINTER(ct.c_ubyte)
    fn.argtypes=[ct.POINTER(Data),ct.c_void_p,ct.POINTER(Grid),ct.POINTER(Pointer),Pointer,byte,byte,Pointer,ct.POINTER(ct.c_longlong)]
    fn.restype=ct.c_int
    ids=[];vectors=[];axes=[];labels=[];source_hash=[]
    for kind in ('unshifted','shifted'):
        for level in LEVELS:
            path=area/f'labels_{level}'/(kind+'.npz')
            if digest(path)!=read(path.parent/'result.json')['cases'][kind]['sha256']:raise ValueError('Wrong source labels')
            source_hash.append(digest(path))
            with np.load(path) as f:
                for key in ('x','zero_label_t','a'):
                    v=f[key];vectors.append(v); axes.append(Grid(len(v),pointer(v)))
                labels.append(f['zero_labels'])
            ids.append(f'{kind}_{level}')
    avec=(Grid*len(axes))(*axes);lvec=(Pointer*len(labels))(*(pointer(v) for v in labels))
    mas=root/'r4_bucket_geometry01/masks.npz'
    if digest(mas)!=read(mas.parent/'geometry.json')['masks_sha256']:raise ValueError('Front mask source')
    with np.load(mas) as f:front=np.ascontiguousarray(f['seam'],dtype=np.uint8)
    t=body.buffers['field_t'];L=body.data.length
    td=np.where(t*L<L-.12,(t*L/(L-.12))**25,1.)-t
    last=np.where((full_order[body.ids]>0).any(axis=(0,1)))[0][-1]
    landed=body.t[:,None,None]+solid[...,1]
    lo=np.minimum.reduce([landed[j:j+len(body.t)-1,i:i+len(body.x)-1,k:k+len(body.a)-1]
                          for i in (0,1) for j in (0,1) for k in (0,1)])
    active=np.ascontiguousarray(lo<=t[min(last+1,len(t)-1)],dtype=np.uint8)
    gp,gw=np.polynomial.legendre.leggauss(degree);gp,gw=(gp+1)/2,gw/2
    body.data.gauss=Grid(degree,pointer(gp));body.data.weights=pointer(gw)
    out=np.zeros((2,7,2));count=ct.c_longlong()
    body.work=lambda d,r:fn(d,r,avec,lvec,pointer(td),active.ctypes.data_as(byte),front.ctypes.data_as(byte),pointer(out),ct.byref(count))
    body.evaluate(q,full_order,load=False,gradient=False)
    data=dict(parent_sha256=CHECKPOINT_SHA,control_sha256=digest(prev/'control.npz'),degree=degree,mask_sha256=digest(mas),
              axes_runs=['material',*ids],label_hashes=source_hash,
              integrals=['stored_energy','explicit_order_force'],regions=['outside_seam','seam'],raw_by_region=out.tolist(),
              ordered_points=count.value)
    summed=out.sum(axis=0,dtype=np.longdouble)
    old=read(prev/f'shape_{degree}/result.json')
    if count.value!=old['ordered_sample_count'] or max(abs(summed[:4,0]-np.asarray(old['crystal_values'])[:,0]))>3e-14:
        raise ArithmeticError('Stored energy or samples disagree with pre-splice reference')
    write_new(case/'result.json',data);print(data,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('area',type=Path);p.add_argument('degree',type=int,choices=(6,8))
    a=p.parse_args();run(a.root.resolve(),a.area.resolve(),a.degree)
