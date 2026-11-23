"""Direct material work across the registered core-to-collar join at fixed outer shape."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss,legval,legder
from scipy.interpolate import RegularGridInterpolator

from mouth_internal_collar import CollarPair,CUT
from mouth_internal_collar_check import read
from mouth_poly_bonded_core import _cut,core_stretches
from mouth_supercore_shape import preparation


def window(t,span):
    start,stop=span
    inside=(t>start)&(t<stop)
    phase=np.pi*(t-start)/(stop-start)
    val=np.where(inside,np.sin(phase)**2,0.)
    derivative=np.where(inside,np.pi*np.sin(2*phase)/(stop-start),0.)
    return val,derivative


def check(state,quad,width,span):
    saved,labels,inner,outer=read(state)
    shape=saved.shape
    ref,_,ref_co=preparation(18)
    chart=CollarPair(shape,saved.co,saved.na,saved.nt,quad)
    # In the collar use reference measure, B=F F^T, and both hoop stretches.
    nodes=labels.reshape(-1,2)[chart.id]; A=(chart.B*nodes[:,:,0]).sum(axis=1)
    Ae,Te=np.einsum('qi,qik->kq',chart.E,nodes)
    At,Tt=np.einsum('qi,qik->kq',chart.T,nodes)
    d=Ae*Tt-At*Te
    f11=(chart.re*Tt-chart.rt*Te)/d; f12=(-chart.re*At+chart.rt*Ae)/d
    f21=-chart.L*Te/d; f22=chart.L*Ae/d
    hoop=chart.r/(shape.neck+A);J=-chart.re*chart.L*hoop**2/d
    mask=chart.inside
    material=[tuple(v[mask] for v in (chart.r,chart.z,chart.qw*(shape.neck+A)**2*d,
                                      f11,f12,f21,f22,hoop,J))]
    # In the core integrate the same physical virtual over the radial band
    # defining its support; no material outside that band can do work.
    q,qw=leggauss(quad)
    t=((np.arange(chart.nt)[:,None]+(q+1)/2)/chart.nt).ravel()
    tw=np.tile(qw/(2*chart.nt),chart.nt)
    tc=np.repeat(t,len(q));b,cut,dt=_cut(shape,chart.co,tc)
    rad=cut-width*np.tile((q+1)/2,len(t))
    qcore,shear,k,Jc=core_stretches(shape,chart.co,ref,ref_co,tc,rad)
    weight=4*np.pi*rad**2*chart.L*np.repeat(tw,len(q))*np.tile(qw*width/2,len(t))/Jc
    material.append((rad,tc,weight,qcore,shear,np.zeros_like(rad),np.full_like(rad,k),qcore,Jc))
    interp=RegularGridInterpolator(inner[:2],inner[2])
    def energy(domain,sign,step,axis):
        r,t,measure,f11,f12,f21,f22,h,J=domain
        b,c,dt=_cut(shape,chart.co,t)
        z=(r-c)/width
        lip=np.maximum(0,1-z*z)
        cut_z=-6*z*lip**2/width
        envelope,rate=window(t,span)
        # Compact spatial displacement, tangent/position continuous through the join.
        val=lip**3*envelope
        dv_r=cut_z*envelope
        dv_w=-cut_z*c*dt*envelope/chart.L+lip**3*rate/chart.L
        step*=sign
        radial=axis=='r'; rr=r+step*radial*val;tt=t+step*(1-radial)*val/chart.L
        m11=1+step*radial*dv_r;m12=step*radial*dv_w
        m21=step*(1-radial)*dv_r;m22=1+step*(1-radial)*dv_w
        hoop=h*rr/r
        Jn=J*(m11*m22-m12*m21)*(rr/r)**2
        n11=m11*f11+m12*f21;n12=m11*f12+m12*f22
        n21=m21*f11+m22*f21;n22=m21*f12+m22*f22
        W=.5*(n11*n11+n12*n12+n21*n21+n22*n22+2*hoop**2-4)-np.log(Jn)+.25*(Jn-1)**2
        bb=shape.profile(tt)[0]
        chi=interp(np.stack((tt,rr/bb),axis=1))
        return (measure*chi)@W
    result={}
    for delta in [1e-4,5e-5]:
        result[str(delta)]={}
        for mode in ('r','w'):
            parts=[];positive=[];negative=[]
            for domain in material:
                pos=energy(domain,1,delta,mode); neg=energy(domain,-1,delta,mode)
                zero=energy(domain,0,delta,mode)
                parts.append(float((pos-neg)/(2*delta)))
                positive.append(float(pos-zero));negative.append(float(neg-zero))
            result[str(delta)][mode]=dict(collar=parts[0],core=parts[1],sum=sum(parts),
                                                 plus_energy=sum(positive),minus_energy=sum(negative),
                                                 plus_parts=positive,minus_parts=negative)
    return result


def trace(state,quad,span):
    """Two one-sided physical Cauchy tractions on the same material join."""
    chart,labels,inner,outer=read(state)
    shape=chart.shape
    ref,_,ref_co=preparation(18)
    q,qw=leggauss(quad)
    t=((np.arange(chart.nt)[:,None]+(q+1)/2)/chart.nt).ravel()
    weights=np.tile(qw/(2*chart.nt),chart.nt)
    j=np.minimum((t*chart.nt).astype(int),chart.nt-1);v=t*chart.nt-j
    b,db,_=shape.profile(t)
    p=np.exp(legval(2*t-1,chart.co))
    pt=2*legval(2*t-1,legder(chart.co))
    for tag,left in [('bonded_join',False),('receiver_edge',True)]:
        e=-.5 if left else .5
        r=b*np.exp(-e*p);slope=r*(shape.length*db/b-e*p*pt)/shape.length
        i=0 if left else 2*chart.na-1;sgn=1. if left else -1.
        tri=chart.triangles[2*(j*2*chart.na+i)+int(not left)]
        nodes=labels.reshape(-1,2)[tri]
        Ae,Te=((nodes[:,1]-nodes[:,0])*(sgn*2*chart.na)).T
        At,Tt=((nodes[:,2]-nodes[:,0])*(sgn*chart.nt)).T
        d=Ae*Tt-At*Te;re=-r*p;rt=slope*shape.length
        f11=(re*Tt-rt*Te)/d;f12=(rt*Ae-re*At)/d
        f21=-shape.length*Te/d;f22=shape.length*Ae/d
        hoop=r/(shape.neck+e);J=-re*shape.length*hoop**2/d
        points=np.stack((t,(r-b)/(CUT-b) if left else r/b),axis=1)
        grid=outer if left else inner
        chi=RegularGridInterpolator(grid[:2],grid[2])(points)
        # n*dS = 4 pi r^2 (1,-slope) L dt points from smaller radius to larger.
        # The sin^2 trace is the same one used in the two-sided material virtual.
        factor=weights*4*np.pi*r**2*shape.length*chi*window(t,span)[0]
        def stress(f11,f12,f21,f22,J):
            Q=-1+.5*J*(J-1);cross=f11*f21+f12*f22
            normal_r=(f11**2+f12**2+Q-cross*slope)/J
            normal_w=(cross-(f21**2+f22**2+Q)*slope)/J
            return np.array([factor@normal_r,factor@normal_w])
        column=stress(f11,f12,f21,f22,J)
        if left:
            yield tag,dict(collar=column.tolist(),order=[float(chi.min()),float(chi.max())])
        else:
            qc,s,k,jc=core_stretches(shape,chart.co,ref,ref_co,t,r)
            core=stress(qc,s,0,k,jc)
            yield tag,dict(collar=column.tolist(),core=core.tolist(),core_minus_collar=(core-column).tolist(),
                           order=[float(chi.min()),float(chi.max())])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('states',nargs='+',type=Path)
    p.add_argument('--quad',type=int,nargs='+',default=[5,8])
    p.add_argument('--width',type=float,nargs='+',default=[.15,.22])
    p.add_argument("--window",nargs="+",choices=["full","below_front"],default=["full","below_front"])
    args=p.parse_args()
    for path in args.states:
        for width in args.width:
            for quad in args.quad:
                for name in args.window:
                    span=(0.,1.) if name=="full" else (.1,.65)
                    print(json.dumps(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                                          width=width,quad=quad,depth_span=span,
                                          signed=check(path,quad,width,span),traction=dict(trace(path,quad,span))),
                                     allow_nan=False),flush=True)
