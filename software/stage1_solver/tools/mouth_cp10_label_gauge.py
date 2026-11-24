"""Independent physical-volume work and orientation checks on the CP10 spatial label grid."""
import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np

from mouth_sealed_shape_relax import Meridian
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
from mouth_last_rings_fluid import GRIDS


def profile(shape, t, dq):
    # Evaluate the stored reference Hermite pieces analytically so complex
    # tangents include displacement of the physical integration points.
    spline=shape.reference
    w=t*shape.unit_length
    i=np.clip(np.searchsorted(spline.x,w.real,side='right')-1,0,len(spline.x)-2)
    a,b,c,d=spline.c[:,i]
    s=w-spline.x[i]
    base=((a*s+b)*s+c)*s+d
    base_t=(3*a*s*s+2*b*s+c)*shape.unit_length
    e=1.5*np.expm1(shape.params[1]+.25*dq[1])
    m=shape.params[2]+.25*dq[2]
    asym=shape.params[3]+.25*dq[3]
    extra=shape.extra+.25*dq[4]
    v=t*(1-t)
    R=base+e*t+4*m*v+4*asym*v*(2*t-1)+16*extra*v*v
    Rt=base_t+e+4*m*(1-2*t)+4*asym*(-6*t*t+6*t-1)+32*extra*v*(1-2*t)
    return 4*R,4*Rt


def depth_grid(length,reference):
    nr,_,nz=GRIDS['double']
    end,origin=length-.12,reference-.12
    s=np.linspace(0.,1.,nr+1)
    lower=origin*s/(1+np.expm1(24*np.log(origin/end))*s**24)**(1/24)
    return np.r_[lower[:-1], end+.12*np.linspace(0,1,nz+1)]/length


VERTICES=np.array(list(product((0,1),repeat=3)))

def gradient_basis(p):
    factors=np.where(VERTICES,p,1-np.asarray(p))
    return np.column_stack([(2*VERTICES[:,d]-1)*np.prod(np.delete(factors,d,axis=1),axis=1)
                            for d in range(3)])


def work(shape,polar,reference,axes,labels,melt,dq):
    length=shape.length*np.exp(dq[0])
    x,t0,a=axes
    t=depth_grid(length,reference)
    if len(t)!=len(t0) or np.max(abs(depth_grid(shape.length,reference)-t0))>3e-15:
        raise ValueError('Incompatible CP10 depth grid')
    i,j,k=np.nonzero(melt)
    widths=np.array([np.diff(x)[i],np.diff(t)[j],np.diff(a)[k]]).T
    # The identity background is the analytic inverse of the normalized mouth
    # chart. Re-evaluate it on the moving knots with the same nodal xi.
    data=np.stack([labels[i+p,j+q,k+r]-
          np.stack([x[i+p],t[j+q],a[k+r]],axis=-1) for p,q,r in VERTICES],axis=1)
    total=0
    for p in product((.5-.5/np.sqrt(3),.5+.5/np.sqrt(3)),repeat=3):
        xx,tt,aa=np.array([x[i],t[j],a[k]])+np.asarray(p)[:,None]*widths.T
        b,bt=profile(shape,tt,dq)
        v=(tt-.06)/.56
        on=(v.real>0)&(v.real<1)
        bump=np.where(on,16*v*v*(1-v)**2,0)
        slope=np.where(on,32*v*(1-v)*(1-2*v)/.56,0)
        h=(polar+dq[5])*bump
        P=(3*np.cos(aa)**2-1)/2
        R=xx*b+xx**3*h*P
        Rx=b+3*xx**2*h*P
        Rt=xx*bt+xx**3*(polar+dq[5])*slope*P
        Ra=-3*xx**3*h*np.cos(aa)*np.sin(aa)
        deriv=np.einsum('svc,vd->scd',data,gradient_basis(p))/widths[:,None,:]
        gr=deriv[:,:,0]/Rx[:,None]
        gw=(deriv[:,:,1]-Rt[:,None]*gr)/length
        ga=(deriv[:,:,2]-Ra[:,None]*gr)/R[:,None]
        weight=4*np.pi/8*np.prod(widths,axis=1)*Rx*R**2*length*np.sin(aa)
        total+=.5*np.sum(weight*np.sum(gr*gr+gw*gw+ga*ga,axis=1))
    return total


def chart_volume_bounds(axes,labels,melt):
    # det(dxi/dp) is at most quadratic separately in each chart coordinate.
    # Its tensor Bernstein hull tests the entire cell, including faces. The
    # values here are floating-point diagnostics, not interval certificates.
    dims=tuple(len(x)-1 for x in axes)
    volumes=np.empty((*dims,3,3,3))
    corners=np.stack([labels[tuple(slice(v,v+n) for v,n in zip(vertex,dims))] for vertex in VERTICES],axis=-2)
    widths=np.broadcast_arrays(*np.ix_(*(np.diff(x) for x in axes)))
    widths=np.stack(widths,axis=-1)
    for ijk in product(range(3),repeat=3):
        deriv=np.einsum('...vc,vd->...cd',corners,gradient_basis(np.array(ijk)/2))/widths[...,None,:]
        volumes[(...,*ijk)]=np.linalg.det(deriv)
    bern=volumes.copy()
    for d in (3,4,5):
        v=np.moveaxis(bern,d,-1)
        v[...,1]=2*v[...,1]-.5*(v[...,0]+v[...,2])
    result={}
    for label,mask in [('active_cells',~melt),('melted_cells',melt)]:
        at=np.unravel_index(np.argmin(np.where(mask[...,None,None,None],volumes,np.inf)),volumes.shape)
        result[label]=dict(cell_count=int(mask.sum()),bernstein_lower_diagnostic=float(bern[mask].min()),
            nonpositive_hull_cells=int((np.min(bern[mask],axis=(1,2,3))<=0).sum()),
            probed_minimum=float(volumes[at]),minimum_cell=list(map(int,at[:3])),
            minimum_cell_fractions=[i/2 for i in at[3:]])
    return result


def measure(root,src,definition):
    from scipy.sparse import load_npz
    cp=root/'r4_melted_survival01/coarse/state_00010.npz'
    if digest(cp)!=CHECKPOINT_SHA:
        raise ValueError('Wrong R4 state')
    level=read(root/'r4_whole_mouth01/coarse/endpoint.json')['depth']
    with np.load(cp,allow_pickle=False) as c:
        shape=Meridian(4.,c['shape_params'],float(c['shape_extra']))
        polar=float(c['x'][5])
    with np.load(src,allow_pickle=False) as f:
        axes=tuple(f[k] for k in ('x','t','a'))
        xi=f['labels']
        melt=f['melted_cells']
        fixed=f['prescribed']
        groups=f['constraint_groups']
        if str(f['source_sha256'])!=CHECKPOINT_SHA:raise ValueError('Wrong spatial field')
        active=f['active_stencil']
    base=np.stack(np.meshgrid(*axes,indexing='ij'),axis=-1)
    K=load_npz(definition)
    v=xi-base
    reaction=(K@v.reshape(-1,3)).reshape(v.shape)
    free_norm=0.
    for c in range(3):
        gc=groups[...,c]
        summed=np.bincount(gc.ravel(),weights=reaction[...,c].ravel())
        free_norm=max(free_norm,float(np.max(abs(summed[gc[~fixed[...,c]]]))))
    unit_K=.5*float(np.sum(v*reaction))
    do=lambda d:work(shape,polar,level,axes,xi,melt,d)
    zero=np.zeros(6)
    unit=float(do(zero))
    if abs(unit-unit_K)>2e-10*max(1.,abs(unit)):
        raise ArithmeticError(f'Gauge quadratures disagree: {unit} / {unit_K}')
    eye=np.eye(6)
    analytic=np.array([-do(1j*1.e-25*e).imag/1.e-25 for e in eye])
    h=5.e-5
    signed=np.array([(do(-h*e)-do(h*e))/(2*h) for e in eye])
    return dict(source_sha256=CHECKPOINT_SHA,labels_sha256=digest(src),unit_work=unit,
        unit_independent_quadrature_delta=unit-unit_K,unit_max_free_residual=free_norm,
        min_signed_chart_determinant=chart_volume_bounds(axes,xi,melt),
        gauge=[dict(epsilon=e,energy=e*unit,max_abs_force_on_order_touching_nodes=e*float(np.max(abs(reaction[active]))),
                    six_gauge_forces=(e*analytic).tolist(),analytic_minus_fd=(e*(analytic-signed)).tolist())
               for e in (0.,1.e-10,1.e-8)])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('root','labels','operator','output'):p.add_argument(name,type=Path)
    args=p.parse_args()
    result=measure(args.root,args.labels,args.operator)
    with args.output.open('x') as f:json.dump(result,f,indent=2,allow_nan=False)
    print(json.dumps(result,allow_nan=False))
