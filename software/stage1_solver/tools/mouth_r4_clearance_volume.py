"""Whole chart volumes independent of the orientation-preserving label diffeomorphism."""
import sys
sys.path.insert(0,'software/stage1_solver/tools')
import json
from pathlib import Path
import numpy as np
from mouth_cp10_label_grid import get_cp10
from mouth_supercore_shape import preparation
from mouth_r4_conforming_prisms import TRIANGLES
base=Path('_scratch/s11c')
body,s,axes=get_cp10(base)
ref=preparation(4)[0]
def integrate(order):
    g,w=np.polynomial.legendre.leggauss(order)
    cuts=np.unique(np.r_[0.,.06,.62,1.,np.linspace(0.,1.,body.data.spline_n+1),np.linspace(0.,1.,ref.reference.x.size)])
    t=(cuts[:-1,None]+(g[None,:]+1)/2*np.diff(cuts)[:,None]).ravel()
    weights=(np.diff(cuts)[:,None]*w[None,:]/2).ravel()
    smooth=np.maximum(0.,np.minimum(1.,(t-.06)/.56))
    bump=body.epsilon*16*smooth**2*(1-smooth)**2
    p=(3*((g+1)/2)**2-1)/2
    rw=body.sh.radius(t*body.sh.length)[0][:,None]+bump[:,None]*p
    current=body.sh.length/3*weights @ (rw**3 @ (w/2))
    ref_core=ref.length/3*weights@((body.cut*body.data.core_ratio*ref.radius(t*ref.length)[0])**3)
    verts=body.buffers['labels'].reshape(body.buffers['reg_t'].size,body.buffers['reg_e'].size,2).copy()
    verts[:,:,0]+=4
    ref_collar=0.
    for offsets in TRIANGLES:
        tri=np.array([verts[dj:dj+verts.shape[0]-1,di:di+verts.shape[1]-1] for dj,di in offsets])
        S=tri[...,0]
        c1,c2=tri[1]-tri[0],tri[2]-tri[0]
        weight=abs(c1[:,:,0]*c2[:,:,1]-c1[:,:,1]*c2[:,:,0])/2
        second=(S[0]**2+S[1]**2+S[2]**2+S[0]*S[1]+S[0]*S[2]+S[1]*S[2])/6
        ref_collar+=np.sum(weight*second)
    return dict(current_over_4pi=current,reference_over_4pi=ref_core+ref_collar,reference_collar_over_4pi=ref_collar,average_J=current/(ref_core+ref_collar))
record={q:integrate(q) for q in (8,12)}
assert abs(record[8]['average_J']-record[12]['average_J'])<1e-12
print(json.dumps(record,allow_nan=False))
