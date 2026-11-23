"""Angular virtual work on a two-sided wound mouth with zero-mean side-wall displacement."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_exchange_carrier import mapped_nodes
from mouth_live_waist_fields import ORDER_INERTIA
from mouth_sealed_shape_relax import Meridian
from mouth_winding_raw_work import basis, topology
from mouth_wound_carrier import order_slopes, straight_vortex


def angular(a, harmonic):
    c = np.cos(a)
    if harmonic == 0:
        return np.ones_like(a), np.zeros_like(a)
    return (3*c*c-1)/2, -3*c*np.sin(a)


def geometry(shape, outer, side, x, w, a, harm, eps):
    b, bw = shape.radius(w)
    t = w/shape.length
    g = 16*(t*(1-t))**2
    gw = 32*t*(1-t)*(1-2*t)/shape.length
    h, dh = angular(a, harm)
    k, kx = (1-x, -1.) if side else (x**3, 3*x*x)
    raw = b+(outer-b)*x if side else b*x
    width = outer-b if side else b
    return (raw+eps*k*g*h, width+eps*kx*g*h,
            bw*(1-x if side else x)+eps*k*gw*h, eps*k*g*dh)


def components(raw, r, rx, rw, ra):
    gradx, gradw, grada = raw
    radial = gradx/rx
    return radial, gradw-rw*radial, (grada-ra*radial)/r


def quad(n):
    x,w = leggauss(n)
    return list(zip((x+1)/2,w/2))


def shape_work(source, degrees, steps):
    path = Path(source)
    g=json.loads((path.parent/'outcome.json').read_bytes())
    limits=json.loads((path.parent/'limit-validation.json').read_bytes())
    if g['exitCode'] or g['childOutcome']['exitCode'] or not limits['verified'] or int(limits['actual']['memory.swap.max']):
        raise ValueError('Unvalidated carrier fields')
    with np.load(path) as z:
        if z['radius']!=4:
            raise ValueError('Expected entrance radius four')
        shape = Meridian(4.,z['shape_params'],float(z.get('shape_extra', 0.)))
        n,k,outer,roof=(z[key].item() for key in ('mesh','polar','outer','roof'))
        N=int(z['winding'])
        order=z['order']; amp=z['amplitude_ratio']; phase=z['phase']
    print(json.dumps(dict(kind='input',sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        shape_params=shape.params.tolist(),extra=shape.extra,radius=4,winding=N,mesh=n,polar=k)),flush=True)
    grids,did=topology(shape,n,k,outer,roof)
    xs,zs=mapped_nodes(n,shape.length,outer,roof)
    angles=np.pi/2*np.linspace(0.,1.,k+1)**1.35
    neutral=straight_vortex(N,outer)
    for nq in degrees:
        points=quad(nq)
        for harmonic in (0,2):
            for h in steps:
                forces=[]
                for eps in (-h,h):
                    energies=[]
                    for side in (0,1):
                        grid=grids[1,side]
                        t,j,l=np.indices((n,n,k))
                        cells=np.stack([grid[t+iz,j+ix,l+ip] for ix,iz,ip in itertools.product((0,1),repeat=3)],-1).reshape(-1,8)
                        v=amp[did[cells]]; p=phase[cells]
                        av=v.mean(axis=1)
                        x0,dx=xs[side][j].ravel(),np.diff(xs[side])[j].ravel()
                        w0,dw=zs[1][t].ravel(),np.diff(zs[1])[t].ravel()
                        a0,da=angles[l].ravel(),np.diff(angles)[l].ravel()
                        x,w,a=x0+dx/2,w0+dw/2,a0+da/2
                        if not side:
                            b,bw=shape.radius(w)
                            # Pull the order field with the varied wall; the uniform
                            # reference order has zero initial polar gradient.
                            gr,gz=order_slopes(order,x,w/shape.length,b,bw,shape.length)
                            order_chart=(gr*b,gz+x*bw*gr,np.zeros_like(gr))
                            grads=components(order_chart,*geometry(shape,outer,side,x,w,a,harmonic,eps))
                        else:
                            grads=np.zeros((3,len(x)))
                        gsq=sum(u*u for u in grads)
                        vol=np.zeros_like(x); rest=np.zeros_like(x)
                        grad_energy=np.zeros_like(x); angleproj=np.zeros_like(x)
                        vgrad=np.zeros_like(x); pot=np.zeros_like(x)
                        for (xa,wx),(ta,wt),(aa,wa) in itertools.product(points,repeat=3):
                            x,w,a=x0+xa*dx,w0+ta*dw,a0+aa*da
                            r,rx,rw,ra=geometry(shape,outer,side,x,w,a,harmonic,eps)
                            if np.min(rx)<=0:
                                raise ValueError('Inverted grid')
                            f,_=neutral(r*np.sin(a))
                            measure=4*np.pi*r*r*np.sin(a)*rx*dx*da*dw*wx*wa*wt
                            vol+=measure
                            rest+=measure*f*f
                            _,df=basis([xa,ta,aa])
                            def physical(field):
                                return components([field@df[0]/dx,field@df[1]/dw,field@df[2]/da],r,rx,rw,ra)
                            vv,pp=physical(v),physical(p)
                            mass=measure*f*f
                            vgrad+=mass*sum(vv1**2 for vv1 in vv)
                            grad_energy+=mass*sum(pp1**2 for pp1 in pp)
                            angleproj+=mass*sum(q*u for q,u in zip(grads,pp))**2
                            pot+=mass*f**8
                        amp2=av*av
                        flow=.5*amp2*(grad_energy-ORDER_INERTIA*angleproj/(rest/vol*amp2+ORDER_INERTIA*gsq))
                        cellpot=pot*(amp2-1)**2*(amp2**3+2*amp2**2+3*amp2+4)/20
                        # f0 grad(f0) on a nonspherical side has a polar normal
                        # component. The sealed inside and outside have independent v.
                        face=grid[:,0 if side else -1,:]
                        radial_index=0. if side else 1.
                        ti,li=np.indices((n,k))
                        vertices=amp[did[np.stack([face[ti+iz,li+ia] for iz,ia in itertools.product((0,1),repeat=2)],-1)]]
                        surf=0.
                        for (ta,wt),(aa,wa) in itertools.product(points,repeat=2):
                            w=zs[1][ti]+ta*np.diff(zs[1])[ti]
                            a=angles[li]+aa*np.diff(angles)[li]
                            r,_,_,ra=geometry(shape,outer,side,radial_index,w,a,harmonic,eps)
                            sn=np.sin(a)
                            f,df=neutral(r*sn)
                            vv=vertices@basis((ta,aa))[0]
                            area=4*np.pi*r*r*sn*np.diff(zs[1])[ti]*np.diff(angles)[li]*wt*wa
                            surf+=np.sum((1-2*side)*.5*area*(vv*vv-1)*f*df*(sn-(ra/r)*np.cos(a)))
                        energies.append([vgrad.sum()/2,cellpot.sum(),flow.sum(),surf])
                    forces.append(energies)
                pieces=(np.array(forces[0])-forces[1])/(2*h)
                print(json.dumps(dict(kind='side_work',mesh=n,polar=k,winding=N,harmonic=harmonic,step=h,
                            quad_order=nq,inside_pieces=pieces[0].tolist(),outside_pieces=pieces[1].tolist(),
                            net_force=pieces.sum(),energy_minus=forces[0],energy_plus=forces[1]),allow_nan=False),flush=True)
                if harmonic==2 and N==0 and nq>=3 and abs(pieces.sum())>1.e-4:
                    raise ArithmeticError('Angularly uniform control has a residual P2 force')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path)
    parser.add_argument('--quad',type=int,nargs='+',default=[3,5])
    parser.add_argument('--step',type=float,nargs='+',default=[.001,.0005])
    args=parser.parse_args()
    shape_work(args.source,args.quad,args.step)
