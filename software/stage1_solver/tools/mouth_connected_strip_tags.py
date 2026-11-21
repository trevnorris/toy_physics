"""Direct weak tag offsets on a flowing, locally homogeneous, connected Morse sheet."""
import argparse

import numpy as np
from scipy.integrate import solve_ivp, quad

from mouth_connected_passage import WAIST, arc_speed, natural_pitch
from mouth_middle_tag_grip import cubic_trace
from mouth_order_wave_pulse import emit
from neutral_patch_sites import WITHIN, morse


def experiment(mesh, start, end, wave, cycles, shape):
    trace, _ = cubic_trace(mesh)
    p = natural_pitch()
    u0 = arc_speed(trace, WAIST)
    a0 = trace.geometry.radius(WAIST)[0]
    # Columns at half-pitch intervals alternate between adjacent hoop rows.
    size = 2*wave*cycles
    index = np.arange(size)
    ref = np.cos(np.pi*index/wave)
    mirrored_distance = np.minimum(index,size-index)/2
    if shape=='localized':
        ref = np.exp(-.5*(mirrored_distance/1.5)**2)
        ref -= np.mean(ref)
    ref /= np.max(ref)-np.min(ref)
    far = wave if shape=='wave' else size//2
    # Include opposing neighbor pairs explicitly; all sites in this strip
    # are mobile in both directions. Periodic hoop translations identify sites
    # at equal material x, not the meridional perturbation or its harmonics.
    cols, rows = np.mgrid[-18:19,-8:9]
    keep=((cols-rows)%2==0)&((cols!=0)|(rows!=0))
    shifts=cols[keep]
    ex=p*shifts/2
    ey=p*rows[keep]*np.sqrt(3)/2
    take=(index[:,None]+shifts[None,:])%size
    depth,alpha,rest=WITHIN
    def deformation(w):
        return float(arc_speed(trace,w)/u0),float(trace.geometry.radius(w)[0]/a0)
    def drift(w, flat):
        x,y = flat.reshape(2,size)
        sx,sy=deformation(w)
        bx,by=ex*sx,ey*sy
        dx,dy=sx*(x[take]-x[:,None]),sy*(y[take]-y[:,None])
        r0=np.hypot(bx,by)
        r=np.hypot(bx+dx,by+dy)
        dr=(dx*(2*bx+dx)+dy*(2*by+dy))/(r+r0)
        weight=np.exp(-alpha*(r0-rest))
        dw=weight*np.expm1(-alpha*dr)
        initial=morse(r0,*WITHIN)[1]/r0
        change=(2*depth*alpha*dw*(1-2*weight-dw)-initial*dr)/r
        side=initial+change
        vx=(bx*change+dx*side).sum(axis=1)/sx
        vy=(by*change+dy*side).sum(axis=1)/sy
        return np.r_[vx,vy]/trace.speed(w)
    cuts=np.r_[start,trace.edges[(trace.edges>start)&(trace.edges<end)],end]
    def propagate(amplitude):
        state=np.r_[np.zeros(size),amplitude*ref]
        for left,right in zip(cuts[:-1],cuts[1:]):
            # The carrier derivative is piecewise smooth. Select the within-cell
            # trace explicitly at each endpoint of an integration interval.
            def rhs(w,z):
                return drift(np.clip(w,np.nextafter(left,right),np.nextafter(right,left)),z)
            result=solve_ivp(rhs,(left,right),state,method='DOP853',rtol=3e-9,atol=2e-14)
            if not result.success: raise ArithmeticError(result.message)
            state=result.y[:,-1]
        xy=state.reshape(2,size)
        if np.max(np.abs(xy.mean(axis=1)))>1e-13:
            raise ArithmeticError('Internal forces moved the mean strip offset')
        return xy
    for h in (1e-6,5e-7):
        plus,minus=propagate(h),propagate(-h)
        odd=(plus-minus)/(2*h)
        relative=float(odd[1,0]-odd[1,far])
        scale=deformation(end)[1]/deformation(start)[1]
        row=dict(test='actual_tags',mesh=mesh,depth=[start,end],transit_time=trace.transit(start,end,8),
                 shape=shape,period_in_pitches=wave*cycles,reference_tag_distance=far/2,
                 amplitude=h,normalized_relative_offset=relative,
                 physical_offset_over_initial=relative*scale,max_absolute_shift=float(np.max(abs(plus))),
                 maximum_even_axial_shift=float(np.max(abs((plus+minus)[0]/2))))
        if shape=='wave':
            def log_rate(w):
                sx,sy=deformation(w)
                bx,by=ex*sx,ey*sy
                r=np.hypot(bx,by)
                t=morse(r,*WITHIN)[1]
                b=np.exp(-alpha*(r-rest))
                radial=2*depth*alpha**2*b*(2*b-1)
                yy=t/r+(radial-t/r)*(by/r)**2
                return -(yy@(1-np.cos(np.pi*shifts/wave)))/trace.speed(w)
            log_gain=sum(quad(log_rate,left,right,epsabs=1e-10)[0] for left,right in zip(cuts[:-1],cuts[1:]))
            expected=np.exp(log_gain)
            row.update(linear_gain=expected,relative_error=abs(relative/expected-1))
            if abs(relative/expected-1)>1e-4:
                raise ArithmeticError('Traveling tags disagree with linear mode')
        emit(**row)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,default=34)
    p.add_argument('--shape',choices=['wave','localized'],required=True)
    p.add_argument('--start',type=float,required=True)
    p.add_argument('--end',type=float,required=True)
    p.add_argument('--cycles',type=int,default=4)
    args=p.parse_args()
    experiment(args.mesh,args.start,args.end,8,args.cycles,args.shape)
