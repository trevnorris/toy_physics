"""Enclose the physical inverse Jacobian of the held R4 cubic map on every Q1 cell.

The kernel certifies det(dxi/dchart). In four physical dimensions its metric
multiplier is |det d(S,T)/d(xi_r,xi_t)| S^2 sin(xi_a)/(L R_x R^2 sin(a)).
The following bounds use the same binary coefficients as the solver.
"""
import argparse
import ctypes as ct
from fractions import Fraction as F
from math import comb
from pathlib import Path
from decimal import Decimal, localcontext
import json

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

from mouth_free_shape_solid import Grid, pointer
from mouth_r4_eulerian_solid import _library

# Extend each normalized polynomial patch slightly past both end faces. This
# includes rounded knot spacings, the rounded weld radius and the chart
# acceptance margin. Coefficient and Bernstein arithmetic below is rational.
PAD = F(1, 1024)


def basis_to_bernstein(n):
    """Power coefficients on [-PAD, 1+PAD] to Bernstein on the unit interval."""
    affine = np.array([[comb(k,j)*(-PAD)**(k-j)*(1+2*PAD)**j if j<=k else F(0)
                        for k in range(n+1)] for j in range(n+1)], dtype=object)
    bern = np.array([[F(comb(i,j),comb(n,j)) if j<=i else F(0)
                      for j in range(n+1)] for i in range(n+1)],dtype=object)
    return bern @ affine


def positive_profile(body, ref):
    """Enclose exactly the core-radius polynomial on each production depth piece."""
    spline = body.buffers['spline']
    count = len(spline)
    exit_shift, width, asymmetry, waist = map(F, body.buffers['shape'][4:] if ref else body.buffers['shape'][:4])
    width_t = F(1,count)
    bound = basis_to_bernstein(4)
    low, high = [], []
    for k,row in enumerate(spline):
        t0 = F(float(k)/count)
        # Work directly in the normalized local variable. Polynomial object
        # evaluation would introduce its floating domain/window transform.
        t = np.array([t0,width_t],dtype=object)
        v = np.array([t0*(1-t0),width_t*(1-2*t0),-width_t**2],dtype=object)
        linear = np.array([width+asymmetry*(2*t0-1),2*width_t*asymmetry],dtype=object)
        co = 16*waist*np.convolve(v,v)
        co[:4] += [F(s)*width_t**p for p,s in enumerate(row[::-1])]
        co[:4] += 4*np.convolve(v,linear)
        co[:2] += exit_shift*t
        co *= 4
        if not all(isinstance(v,F) for v in co):
            raise TypeError('Inexact meridional polynomial coefficients')
        hull = bound @ co
        low.append(min(hull)); high.append(max(hull))
    return min(low),max(high)


def multiply(a,b):
    c = np.full(tuple(np.array(a.shape)+b.shape-1),F(0),dtype=object)
    for (i,j),v in np.ndenumerate(a):
        c[i:i+b.shape[0],j:j+b.shape[1]] += v*b
    return c


def collar_bounds(body):
    eta = body.buffers['reg_e']
    eta = eta[(len(eta)-1)//2:]
    depth = body.buffers['reg_t']
    coeff = body.buffers['collar_coeff']
    det_bound, position_bound = basis_to_bernstein(5), basis_to_bernstein(3)
    ds, positions = [], []
    for j,i in np.ndindex(coeff.shape[:2]):
        s,z = np.vectorize(F,otypes=[object])(coeff[j,i])
        de,dt = F(eta[i+1]-eta[i]),F(depth[j+1]-depth[j])
        power = np.array([F(1),F(2),F(3)],dtype=object)
        det = (multiply(s[:,1:]*power,z[1:]*power[:,None])
               -multiply(z[:,1:]*power,s[1:]*power[:,None]))/(de*dt)
        ds.append(np.min(det_bound @ det @ det_bound.T))
        positions.append(4+np.min(position_bound @ s @ position_bound.T))
    return min(ds),min(positions)


def round_down(value):
    return np.nextafter(float(value),-np.inf)


def round_up(value):
    return np.nextafter(float(value),np.inf)


def cell_values(values):
    return sliding_window_view(values,(2,2,2),axis=(0,1,2))


def chart_by_cell(axes,labels):
    blocks=np.ascontiguousarray(cell_values(labels).transpose(0,1,2,4,5,6,3))
    grids=[[Grid(2,pointer(a[i:i+2])) for i in range(len(a)-1)] for a in axes]
    active=np.ones(1,np.uint8)
    out=np.empty(11)
    fn=_library().r4_chart_orientation
    bounds=np.empty(blocks.shape[:3])
    for i,j,k in np.ndindex(bounds.shape):
        rc=fn(grids[0][i],grids[1][j],grids[2][k],pointer(blocks[i,j,k]),active.ctypes.data,pointer(out))
        if rc or not out[0]>0:
            raise ArithmeticError(('Uncertified chart cell',i,j,k,rc,out.tolist()))
        bounds[i,j,k]=out[0]
    return bounds


def physical_certificate(body,labels,chart_lower):
    axes=tuple(body.buffers['field_'+c] for c in 'xta')
    if ('collar_coeff' not in body.buffers or labels.shape!=tuple(len(a) for a in axes)+(3,)
            or axes[0][0]!=0 or axes[0][-1]!=1 or axes[2][0]!=0
            or np.max(axes[2])>1.6 or np.any(~np.isfinite(labels))):
        raise ValueError('Unrecognized held physical chart')
    x,t,a = np.moveaxis(labels,-1,0)
    if (np.min(x)<0 or np.min(a)<0 or np.max(a)>1.6 or np.min(t)<-2e-13 or np.max(t)>1+2e-13
            or np.max(x)>1+2e-13 or np.any(x[0]!=0) or np.any(a[:,:,0]!=0)):
        raise ValueError('Endpoint falls outside the axis-limit or polynomial enclosure')
    data=body.data
    with localcontext() as c:
        c.prec=50
        emax = -Decimal.from_float(body.cut).ln()/Decimal.from_float(data.pressure)
        emin = -Decimal.from_float(float(np.max(x))).ln()/Decimal.from_float(data.pressure)
    axes_rounding_pad=1e-12
    if (float(emax)>0.5+axes_rounding_pad or float(emin)<-axes_rounding_pad or
            min(np.min(np.diff(body.buffers['reg_e'])),np.min(np.diff(body.buffers['reg_t'])),
                1/data.spline_n)*float(PAD) <= 4*axes_rounding_pad):
        raise ValueError('Polynomial enclosure does not contain the used knot intervals')
    bmin,bmax=positive_profile(body,False)
    refmin,_=positive_profile(body,True)
    amp=abs(F(data.epsilon))
    if bmin<=3*amp or refmin<=0:
        raise ValueError('Held physical radial metric is unresolved')
    det_min,s_min=collar_bounds(body)
    if min(det_min,s_min)<=0:
        raise ValueError('Cubic natural metric is unresolved')
    # R=x(b+x^2 B P); 0<=16y^2(1-y)^2<=1, |P|<=1.
    vol_hi=F(data.length)*(bmax+3*amp)*(bmax+amp)**2
    core=round_down((F(data.core_ratio)*refmin)**3*F(data.ref_length)/vol_hi)
    collar=round_down(s_min*s_min*det_min/(F(data.pressure)*F(float(x.max()))*vol_hi))
    if min(core,collar)<=0: raise ValueError('Reference volume lower bound is not positive')
    # A coordinate minus its vertexwise minimum proportional value is Q1
    # with nonnegative vertices, including at the zero axis. These nodal
    # ratio bounds apply by limits at r=0 and theta=0.
    rx=np.full_like(x,np.inf); rx[1:]=np.nextafter(x[1:]/axes[0][1:,None,None],-np.inf)
    ra=np.full_like(a,np.inf); ra[:,:,1:]=np.nextafter(a[:,:,1:]/axes[2][None,None,1:],-np.inf)
    r_ratio=cell_values(rx).min(axis=(-3,-2,-1))
    a_ratio=cell_values(ra).min(axis=(-3,-2,-1))
    # Concavity gives sin(a')/sin(a)>=min(1,a'/a) on [0,pi/2].
    # The decimal cutoff is strictly below pi/2. If a label exceeds it,
    # sin(a')>=a'(1-a'^2/6) on [0,1.6] supplies a safe alternative.
    cutoff=F(1570796326794896619,10**18)
    if F(float(max(a.max(),axes[2][-1])))<=cutoff:
        sine=np.minimum(1.,a_ratio)
        angular_rule='concavity'
    else:
        scale=round_down(1-F(float(a.max()))**2/6)
        sine=np.nextafter(scale*a_ratio,-np.inf)
        angular_rule='sine Taylor lower bound'
    if min(r_ratio.min(),sine.min())<=0: raise ValueError('Radial or angular ratio has no positive enclosure')
    # Collar volume already bounds S directly, so only x_phys <= 1 is
    # required there. The normalized core volume uses (xi_r/x_phys)^2.
    quad=np.nextafter(core*r_ratio,-np.inf)
    core_cell=np.nextafter(quad*r_ratio,-np.inf)
    rings=cell_values(x)
    interface=body.cut
    factor=np.where(rings.max(axis=(-3,-2,-1))<=interface,core_cell,
                    np.where(rings.min(axis=(-3,-2,-1))>interface,collar,np.minimum(core_cell,collar)))
    chart=chart_by_cell(axes,labels)
    if float(chart.min())!=chart_lower:raise ArithmeticError('Independent per-cell chart pass differs from whole mesh')
    limits=np.nextafter(np.nextafter(chart*factor,-np.inf)*sine,-np.inf)
    idx=np.unravel_index(np.argmin(limits),limits.shape)
    return dict(physical_global_lower=float(limits[idx]),cell=list(map(int,idx)),cell_count=int(limits.size),
                chart_global_lower=float(chart.min()),limiting_chart_lower=float(chart[idx]),
                limiting_radial_metric_lower=float(factor[idx]),limiting_sine_ratio_lower=float(sine[idx]),
                global_xi_r_over_x_lower=float(r_ratio.min()),global_sine_ratio_lower=float(sine.min()),
                angular_rule=angular_rule,collar_patch_determinant_lower=round_down(det_min),
                collar_S_lower=round_down(s_min),physical_profile_bounds=[round_down(bmin),round_up(bmax)],
                reference_core_profile_lower=round_down(refmin),core_metric_lower=core,collar_metric_lower=collar,
                coverage='Both one-sided weld limits; regular coordinate-axis limits by Q1 nodal inequalities.',
                interpretation='Conservative physical inverse-determinant lower bound over the entire held body.')


if __name__=='__main__':
    from mouth_cp10_label_grid import get_cp10
    from mouth_r4_eulerian_solid import chart_bounds
    from mouth_r4_smooth_collar_reference import enable_smooth_collar
    from mouth_exact_checkpoint_control import digest
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path)
    parser.add_argument('state',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    body,cp,axes=get_cp10(args.root)
    with cp, np.load(args.state,allow_pickle=False) as s:
        if int(s['radius'])!=4 or int(s['winding'])!=1 or not np.array_equal(cp['x'],s['x']):
            raise ValueError('Physical bound is valid only on the held R4 shape')
        labels=s['labels']
    enable_smooth_collar(body)
    active=np.ones(tuple(len(a)-1 for a in axes),np.uint8)
    chart=chart_bounds(axes,labels,active)['chart_det_global_lower']
    result=physical_certificate(body,labels,chart)
    result.update(state_sha256=digest(args.state),script_sha256=digest(Path(__file__)))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False)
    print(json.dumps(result,allow_nan=False),flush=True)
