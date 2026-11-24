"""Bracket the physical volume-ratio infimum on the open set chi > 0 at R4.

The upper endpoint is a positive-hat witness, with the polar, meridional
and registered-reference Jacobians evaluated in interval arithmetic. The
lower endpoint uses a uniform physical multiplier bound and a Bernstein
subdivision certificate. It need not be attained at a point of positive chi.
"""
import ctypes as ct
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from pathlib import Path
import subprocess

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import numpy as np

from mouth_active_prism_orientation import (
    _Support, _coefficients, _dot_vertices, _cross, _sub, _add, _down, _up,
    _possible_fold, _quadratic_min, _radial_bound, _split, _exact_array, UNIT)
from mouth_free_shape_solid import Data, Pointer, pointer
from mouth_r4_fold_constitutive import Prism, native_callback, point

IV=MPIntervalContext()
IV.dps=50


def iv(x):
    q=Q(x)
    return IV.mpf(q.numerator)/q.denominator


@lru_cache(maxsize=4)
def kernels(folder):
    import mouth_free_shape_solid
    include=Path(mouth_free_shape_solid.__file__).resolve().parent
    folder=Path(folder).resolve()/'ordered_J_reader'
    folder.mkdir(exist_ok=True,parents=True)
    library=folder/'ordered_J.so'
    source=Path(__file__).with_suffix('.cpp')
    subprocess.run(['g++','-O3','-std=c++17','-shared','-fPIC','-I',str(include),
                    str(source),'-o',str(library)],check=True)
    fn=ct.CDLL(str(library)).ordered_j_candidates
    fn.argtypes,fn.restype=[ct.POINTER(Data),Pointer],None
    return fn,native_callback(include,folder)


def bind(trial,q,u):
    b=trial.b
    offsets=b.unpack(q)
    order=np.ascontiguousarray(u[2,b.ids])
    b.data.offset,b.data.correction=pointer(offsets),pointer(order)
    return offsets,order


def interpolants(pts,bary,y):
    w=np.array((1-Q(bary[1])-Q(bary[2]),Q(bary[1]),Q(bary[2])),object)
    y=Q(y)
    return w@((1-y)*pts[0]+y*pts[1]),w,y


def nodes(prism):
    labels,offsets=np.moveaxis(np.array(prism.nodes,float),2,0)
    ref=_exact_array(labels).transpose(1,0,2)
    return ref,ref+_exact_array(offsets).transpose(1,0,2)


def profile(body,t,ref):
    s=body.buffers['spline']; n=len(s)
    k=min(n-1,max(0,int(n*t))); h=t-Q(k,n)
    co=[Q(v) for v in s[k]]
    shape=body.buffers['shape'][4:] if ref else body.buffers['shape'][:4]
    e,m,a,quart=map(Q,shape); prod=t*(1-t)
    return 4*((((co[0]*h+co[1])*h+co[2])*h+co[3])+e*t+4*prod*(m+a*(2*t-1))+16*quart*prod**2)


def reference_factors(body,x,t,core):
    """Intervals for all consistent reference-triangle choices including edge ties."""
    d=body.data
    if core:
        sx=profile(body,t,True)*Q(d.core_ratio)
        return [(iv(x*sx),iv(sx*Q(d.ref_length)))]
    er,tr=body.buffers['reg_e'],body.buffers['reg_t']
    labels=body.buffers['labels'].reshape(len(tr),len(er),2)
    e=-IV.ln(iv(x))/iv(d.pressure)
    def cell(grid,z):
        return min(len(grid)-2,max(0,int(np.searchsorted(grid,z,side='right'))-1))
    j=cell(tr,float(t))
    if j and t<Q(tr[j]): j-=1
    if j<len(tr)-2 and t>=Q(tr[j+1]): j+=1
    dt=Q(tr[j+1])-Q(tr[j]); v=(t-Q(tr[j]))/dt
    possible=[]
    for i in range(cell(er,_down(float(e.a))),cell(er,_up(float(e.b)))+1):
        de=Q(er[i+1])-Q(er[i]); u=(e-iv(er[i]))/iv(de)
        for up in (False,True):
            if (up and u+iv(v)<1) or (not up and u+iv(v)>1):continue
            z0, ze, zt=(_exact_array(labels[tj,ei]) for tj,ei in
                        ((j+up,i+up),(j+up,i+1-up),(j+1-up,i+up)))
            s=1-2*up
            u1,v1=(1-u,1-v) if up else (u,v)
            S=iv(4+z0[0])+u1*iv(ze[0]-z0[0])+iv(v1*(zt[0]-z0[0]))
            se,te=s*(ze-z0)/de; st,tt=s*(zt-z0)/dt
            det=abs((se*tt-st*te)/(Q(d.pressure)*x))
            if S<=0 or det<=0:raise ValueError('Invalid reference factor')
            possible.append((S,iv(det)))
    return possible


def witness(prism,bary,angle,body,support):
    ref,pts=nodes(prism)
    X,w,y=interpolants(pts,bary,angle)
    x,t,a=interpolants(ref,bary,angle)[0]
    if x<=0 or not (w>=0).all() or not 0<=y<=1: return None
    chi=support.at(X)
    if chi<=0: return None
    scale=np.prod([Q(hi)-Q(lo) for lo,hi in prism.grid]); det=w@_coefficients(pts)@np.array([(1-y)**2,2*y*(1-y),y*y])/scale
    if det<=0:raise ArithmeticError('Readout found positive order at a verified nonpositive chart determinant')
    xx,tt,aa=X
    # Exact binary constants match the PDE's native normalized geometry.
    h=(tt-Q(.06))/Q(.56)
    f=16*h*h*(1-h)**2*Q(body.epsilon) if 0<h<1 else Q(0)
    p=(3*IV.cos(iv(aa))**2-1)/2
    b=profile(body,tt,False)
    R=iv(xx*b)+iv(xx**3*f)*p
    rx=iv(b)+iv(3*xx**2*f)*p
    if a==0:
        angular=iv(w@(pts[1,:,2]-pts[0,:,2])/(Q(prism.grid[2][1])-Q(prism.grid[2][0])))
    else:
        angular=IV.sin(iv(aa))/IV.sin(iv(a))
    multiplier=rx*iv(body.data.length)*R**2*angular
    vals=[multiplier*iv(det)/(S**2*det2) for S,det2 in
          reference_factors(body,x,t,prism.i<body.cut_node)]
    if not vals:raise ValueError('Missing reference-triangle interval')
    return dict(cell_t_x_angle=[prism.j,prism.i,prism.k],triangle=prism.triangle,barycentric=list(map(float,w)),
                angular_fraction=float(y),chi=float(chi),chart_J=float(det),placement_x_t_angle=list(map(float,X)),
                J_interval=[float(_down(min(float(v.a) for v in vals))),float(_up(max(float(v.b) for v in vals)))])


def positive_multiplier(body,off):
    """Uniform lower bound for rx L (R/S)^2 sin(a')/(sin(a) det_ref).

    Ratios x'/x and a'/a are bounded by nodal ratios (convex weights).
    The estimate remains valid in the axis limits. A deliberately loose
    sine inequality and upper bounds on the reference chart only widen it.
    """
    d=body.data
    lo=Q(_radial_bound(body)); ref=body.buffers['shape'][4:]; c=body.buffers['spline']; n=len(c)
    if lo<=0:raise ArithmeticError('No positive physical chart bound')
    m=max(sum(abs(Q(z))/n**k for k,z in enumerate(row[::-1])) for row in c)
    # 4t(1-t) and 16t²(1-t)² never exceed one on [0,1].
    core=4*(m+sum(abs(Q(z)) for z in ref))*Q(d.core_ratio)
    er,tr=body.buffers['reg_e'],body.buffers['reg_t']
    if er[0]>0 or tr[0]!=0 or tr[-1]!=1:raise ValueError('Unexpected registration domain')
    cut=Q(body.x[body.cut_node]); p=Q(d.pressure); emin=min(Q(b)-Q(a) for a,b in zip(er,er[1:]))
    if not -IV.ln(iv(cut))/iv(p)<iv(Q(er[-1])+emin):raise ValueError('Unsupported registration extrapolation')
    labels=body.buffers['labels'].reshape(len(tr),len(er),2)
    # With less than one e-cell of extrapolation, the sum of absolute
    # barycentric weights is at most five even on the outermost triangle.
    S_over_x=max(core,(4+5*max(abs(Q(v)) for v in labels[:,:,0].flat))/cut)
    le=_exact_array(labels); de=np.diff(_exact_array(er)); dt=np.diff(_exact_array(tr))
    se=np.diff(le,axis=1)/de[None,:,None]; st=np.diff(le,axis=0)/dt[:,None,None]
    det_max=max(core*Q(d.ref_length),max(np.max(abs(se[:-1,:,0]*st[:,:-1,1]-se[:-1,:,1]*st[:,:-1,0])),
            np.max(abs(se[1:,:,0]*st[:,1:,1]-se[1:,:,1]*st[:,1:,0])))/(p*cut))
    ratios=[]
    for values,delta in [(body.x[None,1:,None],off[:,1:,:,0]),(body.a[None,None,1:],off[:,:,1:,2])]:
        # Conservative outward-rounded nodal division vs exact known positive denominator.
        ratios.append(Q(float(_down(1+_down(delta/values)).min())))
    if np.any(off[:,0,:,0]) or np.any(off[:,:,0,2]) or min(ratios)<=0:
        return Q(0)
    # sin(a') >= a'/2 and sin(a)<=a throughout the physical angle range.
    return lo*Q(d.length)*(lo*ratios[0]/S_over_x)**2*ratios[1]/(2*det_max)


def remaining_prisms(body,off,threshold):
    """Vectorized outward-rounded Bernstein screening for det(chart)>threshold."""
    axes=(body.x,body.t,body.a)
    grid=np.stack(np.broadcast_arrays(body.x[None,:,None],body.t[:,None,None],body.a[None,None,:]),axis=-1)
    coords=grid+off
    denom=np.ones((len(body.t)-1,len(body.x)-1,len(body.a)-1))
    for delta in (_up(np.diff(body.x))[None,:,None],_up(np.diff(body.t))[:,None,None],_up(np.diff(body.a))[None,None,:]):
        denom=_up(denom*delta)
    cutoff=_up(float(threshold)*denom)[...,None]
    for tri,offsets in enumerate((((0,0),(0,1),(1,0)),((1,1),(1,0),(0,1)))):
        v=np.stack([coords[j:j+len(body.t)-1,i:i+len(body.x)-1] for j,i in offsets],axis=-2)
        bottom=(_down(v[:,:,:-1]),_up(v[:,:,:-1])); top=(_down(v[:,:,1:]),_up(v[:,:,1:]))
        edges=[[_sub(tuple(v[...,i,:] for v in plane),tuple(v[...,0,:] for v in plane)) for i in (1,2)]
               for plane in (bottom,top)]
        (x0,t0),(x1,t1)=edges
        da=_sub(top,bottom)
        mixed=[_add(a,b) for a,b in zip(_cross(x0,t1),_cross(x1,t0))]
        ok=np.all((_dot_vertices(_cross(x0,t0),da)[0]>cutoff)&(_dot_vertices(_cross(x1,t1),da)[0]>cutoff)&
                  (_dot_vertices(mixed,da)[0]>2*cutoff),axis=-1)
        for j,i,k in np.argwhere(~ok):
            yield Prism(body,off,list(map(int,(j,i,k))),tri)


def chart_margin(body,off,support):
    # An undecided positive margin merely widens the readout. Zero remains
    # a lower bound because caller certifies strict J>0 on actual chi>0.
    for exponent in (6,12,24,40,54,68):
        delta=Q(1,2**exponent); good=True
        for prism in remaining_prisms(body,off,delta):
            ref,pts=nodes(prism)
            if not support.positive[support.box(pts)].any():continue
            scale=np.prod([Q(hi)-Q(lo) for lo,hi in prism.grid])
            coeff=_coefficients(pts)/scale-delta
            ref=np.array([[list(b)+[Q(y)] for b in UNIT] for y in (0,1)],object)
            stack=[(pts,coeff,ref,0)]
            while stack:
                p,c,r,d=stack.pop()
                minima=[_quadratic_min(row)[0] for row in c]
                if min(minima)>0:continue
                possible=_possible_fold(minima)
                box=support.box(possible@p)
                if not support.positive[box].any():continue
                if d==10:
                    good=False;break
                for child in _split(p,c,r,support,box):
                    stack.append((*child,d+1))
            if not good:break
        if good:return delta
    return Q(0)


def ordered_prism_point(trial,q,u,folder,witness=None):
    """Evaluate the retained physical law at the old or a supplied material witness."""
    b=trial.b
    off,order=bind(trial,q,u)
    site=witness or dict(cell_t_x_angle=[5,18,0],triangle=0,
                         barycentric=[11/24,1/3,5/24],angular_fraction=43/64)
    prism=Prism(b,off,site['cell_t_x_angle'],site['triangle']); bary=site['barycentric'];y=site['angular_fraction']
    _,fn=kernels(Path(folder).resolve())
    # Keep arrays backing b.data alive until the native read finishes.
    if order.shape!=b.ids.shape:raise ValueError('Malformed crystal order')
    ref,pts=nodes(prism);loc,_,_=interpolants(pts,bary,y)
    chi=_Support(b,u[2]).at(loc)
    with mp.workdps(70):
        if prism.ref_point(bary,y)[0]==0.:
            return dict(**site,chi=float(chi),J=None,W=None,chi_W=None,scope='exact central-axis limit unresolved by point probe')
        row=point(prism,bary,y,trial,u,fn)
    return dict(cell_t_x_angle=list(site['cell_t_x_angle']),triangle=site['triangle'],
                barycentric=bary,angular_fraction=y,chi=float(chi),J=row['J'],W=row['W'],
                chi_W=float(chi)*row['W'] if row['W'] is not None else None,
                physical_r_w_theta=row['physical_r_z_theta'],principal_stretches=row['principal_stretches'])


def ordered_j_readout(trial,q,u,folder):
    """Certify strict orientation; return a bracket for inf J over positive order."""
    b=trial.b
    if trial.certificate(q,u[2],False)['status']!='certified':
        raise ArithmeticError('Cannot report active-J infimum on infeasible state')
    off,order=bind(trial,q,u)
    support=_Support(b,u[2])
    fn,_=kernels(Path(folder).resolve())
    samples=np.empty((len(b.t)-1,len(b.x)-1,len(b.a)-1,2,4))
    fn(ct.byref(b.data),pointer(samples))
    # Keep backing storage live through the native call.
    if order.shape!=b.ids.shape:raise ValueError('Invalid order storage')
    values=samples[...,0]
    ranked=np.argsort(values,axis=None)
    best=None
    for index in ranked[:16]:
        j,i,k,tri=map(int,np.unravel_index(index,values.shape))
        if not np.isfinite(values[j,i,k,tri]):break
        _,b1,b2,y=samples[j,i,k,tri]
        prism=Prism(b,off,[j,i,k],tri)
        # Native candidate search is only a locator. Very small interpolation
        # weights can be affected by rounding at a knot. Verify this concrete
        # witness in exact rational Q1 and enclose its physical multiplier.
        candidates=[(b1,b2,y)]
        for h in (2.**-36,2.**-18):
            candidates.extend((b1+e*h,b2+f*h,y+g*h) for e,f,g in product((-1,0,1),repeat=3))
        for b1p,b2p,yp in candidates:
            if min(b1p,b2p,yp)<0 or b1p+b2p>1 or yp>1:continue
            row=witness(prism,[1-b1p-b2p,b1p,b2p],yp,b,support)
            if row:
                if best is None or row['J_interval'][1]<best['J_interval'][1]:best=row
                break
    if best is None:raise ArithmeticError('No verified ordered physical sample')
    margin=chart_margin(b,off,support)
    multiplier=positive_multiplier(b,off)
    lower=max(0.,float(_down(float(margin*multiplier))))
    upper=best['J_interval'][1]
    if not 0<=lower<=upper:raise ArithmeticError('Inconsistent active J bracket')
    local=ordered_prism_point(trial,q,u,folder)
    return dict(infimum_bounds=[lower,upper],certified_chart_margin=float(margin),
                physical_multiplier_lower_bound=max(0.,float(_down(float(multiplier)))),
                upper_witness=best,prism_5_18_0=local,machine_zero_witness=bool(upper<=64*np.finfo(float).eps),
                machine_zero_threshold=64*np.finfo(float).eps,
                scope='Continuum of each discrete prism where interpolated chi > 0, including one-sided boundary limits. Lower bound is conservative. Upper bound is a verified physical sample. An attained minimum is not assumed.')
