"""Registered bicubic natural collar map on the held radius-four preparation.

Each coordinate interpolates every natural registry point. Clamped endpoint
derivatives are the adjacent segment slopes in the corresponding chart variable.
The weld is affine, so interpolation preserves its full fixed-label trace.
The other edges are straight with nonuniform registration; the geometric
traces and unit tangents persist while inter-knot material spacings change.
"""
import argparse
import json
from math import comb
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.signal import convolve2d

from mouth_free_shape_solid import pointer


def interpolate(axis, values):
    """Normalized ascending local power coefficients, indexed by cell then power."""
    width=np.diff(axis)
    ends=((1,(values[1]-values[0])/width[0]),(1,(values[-1]-values[-2])/width[-1]))
    c=CubicSpline(axis,values,axis=0,bc_type=ends).c[::-1].swapaxes(0,1)
    return c*width[:,None].reshape((len(width),1)+(1,)*(c.ndim-2))**np.arange(4).reshape(
        (1,4)+(1,)*(c.ndim-2))


def bernstein(c):
    """Power-to-Bernstein transformation for each of two unit-interval axes."""
    for axis in (0,1):
        n=c.shape[axis]-1
        b=np.array([[comb(i,j)/comb(n,j) if j<=i else 0 for j in range(n+1)] for i in range(n+1)])
        c=b@c if axis==0 else c@b.T
    return c


def coefficients(eta,t,registry):
    middle=(len(eta)-1)//2
    if eta[middle]!=0 or eta[-1]!=.5 or t[0]!=0 or t[-1]!=1:
        raise ValueError('Invalid radius-four collar chart extent')
    e=eta[middle:]
    nodes=np.array(registry[:,middle:],copy=True)
    # The first stored component is S-4. Keep that offset out of the fit.
    along_e=interpolate(e,nodes.transpose(1,0,2)).transpose(2,0,1,3)
    poly=interpolate(t,along_e).transpose(0,2,4,1,3).copy()
    return e,poly


def certify(e,t,poly):
    d_bounds=[np.inf,-np.inf]
    first=[np.inf,np.inf]
    for j,row in enumerate(poly):
        for i,(s,z) in enumerate(row):
            def de(v):
                return v[:,1:]*np.arange(1,4)[None,:]/(e[i+1]-e[i])
            def dt(v):
                return v[1:,:]*np.arange(1,4)[:,None]/(t[j+1]-t[j])
            det=bernstein(convolve2d(de(s),dt(z))-convolve2d(dt(s),de(z)))
            d_bounds[0]=min(d_bounds[0],det.min())
            d_bounds[1]=max(d_bounds[1],det.max())
            first[0]=min(first[0],bernstein(de(s)).min())
            first[1]=min(first[1],bernstein(dt(z)).min())
    if min(d_bounds[0],*first)<1.e-12:
        raise ValueError(f'Could not certify ordered orientation of the reference map: {d_bounds}, {first}')
    return dict(registered_map_det_bernstein_bounds=d_bounds,
                dS_deta_dT_dt_bernstein_lower=first)


def enable_smooth_collar(body):
    """Opt in a single geometry to the registered natural map; preserve all held fields."""
    e,poly=coefficients(body.buffers['reg_e'],body.buffers['reg_t'],body.buffers['labels'])
    proof=certify(e,body.buffers['reg_t'],poly)
    body.buffers['collar_coeff']=np.ascontiguousarray(poly,dtype=np.float64)
    body.data.collar_coeff=pointer(body.buffers['collar_coeff'])
    return proof


def measure(root):
    from mouth_r4_collar_boundary_geometry import audit
    cert=audit(root)
    with np.load(root/'supercore_joined_front01/adopt309/state.npz',allow_pickle=False) as state:
        z=state['labels']
    t=np.linspace(0,1,z.shape[0]); eta=np.linspace(-.5,.5,z.shape[1])
    e,poly=coefficients(eta,t,z)
    cert.update(certify(e,t,poly))
    nm=e.size-1
    nodes=z[:,-e.size:]
    err=0.; discontinuities=np.zeros(3)
    for j,row in enumerate(poly):
        for i,p in enumerate(row):
            for a in (0,1):
                for b in (0,1):
                    got=p[:,0 if a==0 else slice(None),0 if b==0 else slice(None)].reshape(2,-1).sum(axis=1)
                    err=max(err,np.max(abs(got-nodes[j+a,i+b])))
            for axis,neighbor,h in ((0,poly[j+1,i] if j+1<len(poly) else None,np.diff(t)),
                                     (1,poly[j,i+1] if i+1<nm else None,np.diff(e))):
                if neighbor is None:
                    continue
                ix=(j,i)[axis]
                a=p.swapaxes(1,axis+1); b=neighbor.swapaxes(1,axis+1)
                for k in range(3):
                    scale=np.array([[1,1,1,1],[0,1,2,3],[0,0,2,6]][k])
                    left=np.einsum('r,vrq->vq',scale,a)/h[ix]**k
                    right=b[:,k,:]*scale[k]/h[ix+1]**k
                    discontinuities[k]=max(discontinuities[k],np.max(np.sum(abs(left-right),axis=-1)))
    weld=poly[:,-1].sum(axis=-1)
    target=np.zeros_like(weld);target[:,0,0]=.5
    target[:,1,0]=t[:-1];target[:,1,1]=np.diff(t)
    cert.update(method='tensor clamped cubic through all registered collar nodes, exact endpoint secants',
                knot_max_abs=float(err),weld_polynomial_max_sum_abs=float(np.max(np.sum(abs(weld-target),axis=-1))),
                full_face_jump_coefficient_l1_bounds_for_normal_order_0_1_2=discontinuities.tolist(),
                cells=int(len(t[:-1])*nm))
    return cert


def boundary_extrema(root):
    """Stationary points of cubic position errors and quadratic tangent errors on every edge segment."""
    from numpy.polynomial import Polynomial
    with np.load(root/'supercore_joined_front01/adopt309/state.npz',allow_pickle=False) as state:
        z=state['labels']
    t=np.linspace(0,1,z.shape[0]); eta=np.linspace(-.5,.5,z.shape[1])
    e,c=coefficients(eta,t,z); knots=z[:,-len(e):]
    traces={
        'weld':(t,knots[:,-1],c[:,-1].sum(axis=-1),1),
        'lateral_wall':(t,knots[:,0],c[:,0,:,:,0],1),
        'entrance':(e,knots[0],c[0,:,:,0,:],0),
        'exit':(e,knots[-1],c[-1].sum(axis=-2),0)}
    out={}
    for name,(q,nodes,trace,axis) in traces.items():
        pos=(0.,0.); derivative=(0.,0.); perpendicular=0.
        for i,row in enumerate(trace):
            dq=q[i+1]-q[i]
            other=row[1-axis].copy();other[0]-=nodes[i,1-axis]
            perpendicular=max(perpendicular,float(abs(other).sum()))
            err=Polynomial(row[axis])-Polynomial([nodes[i,axis],nodes[i+1,axis]-nodes[i,axis]])
            for case,poly in enumerate([err,err.deriv()/dq]):
                candidates=[0.,1.]+[float(a.real) for a in poly.deriv().roots() if a.imag==0 and 0<a.real<1]
                best=max((abs(float(poly(a))),q[i]+dq*a) for a in candidates)
                if case==0:pos=max(pos,best)
                else:derivative=max(derivative,best)
        out[name]={'max_along_boundary_position_shift':pos[0],'position_witness_label':pos[1],
                   'max_parameterized_tangent_shift':derivative[0],'tangent_witness_label':derivative[1],
                   'max_transverse_shift_coefficient_bound':perpendicular}
    return dict(bounds='along-boundary extrema from derivative roots on every bicubic boundary segment',
                boundaries=out)


if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('root',type=Path);a.add_argument('destination',type=Path)
    a.add_argument('--boundary-extrema',action='store_true')
    opts=a.parse_args()
    record=boundary_extrema(opts.root) if opts.boundary_extrema else measure(opts.root)
    opts.destination.parent.mkdir(parents=True,exist_ok=True)
    opts.destination.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    print(json.dumps(record,allow_nan=False))
