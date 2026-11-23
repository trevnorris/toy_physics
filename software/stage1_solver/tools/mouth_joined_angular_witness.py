"""Angular virtual work of a resolved two-bath state, supported upstream of the front."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_joined_front_audit import Front
from mouth_supercore_front import check_guard
from mouth_breathing_mobile_patch import ORDER_INERTIA,KAPPA,WELL
from mouth_wound_carrier import straight_vortex


def bump(w,length):
    # Zero before t=.10 and past t=.70. Both displacement and slope vanish.
    u=np.clip((w/length-.10)/.60,0.,1.)
    return 16*(u*(1-u))**2,32*u*(1-u)*(1-2*u)/(.60*length)


def deform(model,r,z,angle,side,step):
    b,bw=model.shape.radius(z)
    width=28.-b if side else b
    x=(r-b)/width if side else r/b
    old_slope=bw*(1-x if side else x)
    k,kx=(1-x,-1.) if side else (x**3,3*x*x)
    g,gz=bump(z,model.shape.length)
    c=np.cos(angle)
    P=(3*c*c-1)/2
    d=P*step*g
    rnew=r+k*d
    stretch=1+kx*d/width
    new_slope=old_slope+step*k*gz*P
    theta_slope=-3*step*k*g*c*np.sin(angle)
    return rnew,stretch,old_slope,new_slope,theta_slope


def check(directory, inputs=None, grid=None):
    live=directory/'live.npz'
    rows=(json.loads(line) for line in (directory/'stdout').read_text().splitlines())
    if inputs is None:
        meta=next(row for row in rows if row.get('kind')=='preparation')
        source=Path(meta['state'])
    else:
        meta=next(row for row in rows if row.get('kind')=='same_lattice')
        source=inputs
        meta['grid']=grid
        emit=next(row for row in rows if row.get('kind')=='fields')
        if hashlib.sha256(live.read_bytes()).hexdigest()!=emit['sha256']:
            raise ValueError('Mismatched archived reservoir fields')
    if hashlib.sha256(source.read_bytes()).hexdigest()!=meta['source_sha256']:
        raise ValueError('Material and resolved state differ')
    check_guard(live)
    front=Front(source,tuple(meta['grid']),meta['winding'])
    m=front.joined(front.sh,False)
    with np.load(live) as f:
        u=f['fields']
    wt,r,ang,z,F,Fp,s,cut,liveband,area,elastic=m.D.T
    nt=len(m.z[1])
    t=m.V@u[0]
    grads=np.array([np.stack([G@field for G in m.G],-1) for field in u])
    chi=(m.V@u[2])+m.ref[:,0]
    grads[2]+=m.ref[:,1:]
    vortex=straight_vortex(m.winding,28.)
    if inputs is None:
        base,gr=m.evaluate(u)
        # Test a nonstationary trial direction, so a near-zero gradient at equilibrium
        # cannot disguise wrong signs in derivatives of amplitude, phase or order.
        rng=np.random.default_rng(18371)
        trial=u+.0002*rng.normal(size=u.shape)*m.masks
        direction=rng.normal(size=u.shape)*m.masks
        value,g=m.evaluate(trial)
        inc=.000001
        gradient_diff=(m.evaluate(trial+inc*direction)[0]-m.evaluate(trial-inc*direction)[0])/(2*inc)
        analytic=(direction*g).sum()
        print(json.dumps(dict(kind='adjoint',grid=meta['grid'],finite=gradient_diff,analytic=analytic,
                              input_hashes=[meta['source_sha256'],hashlib.sha256(live.read_bytes()).hexdigest()])),flush=True)
        if abs(gradient_diff-analytic)>3.e-5*(1+abs(analytic)):
            raise ArithmeticError('Trial energy/gradient mismatch')

    def moved(step):
        energy=0.
        for side in (0,1):
            box=m.ids[1,side]
            # The Joined matrix places all eight integration-point blocks
            # contiguously for each region.
            ncell=np.prod(np.subtract(box.shape,1))
            _,qstart=m.blocks[1,side]
            start=sum(8*np.prod(np.subtract(m.ids[k].shape,1)) for k in list(m.ids)[:qstart//8])
            ids=slice(int(start),int(start+8*ncell))
            a,rr,w=ang[ids],r[ids],z[ids]
            rn,stretch,oldw,neww,at=deform(m,rr,w,a,side,step)
            def transport(g):
                radial=g[:,0]/stretch
                return np.column_stack([radial,(rr*g[:,1]-at*radial)/rn,
                                        g[:,2]+oldw*g[:,0]-neww*radial])
            gv,gp,gc=[transport(g[ids]) for g in grads]
            reference=transport(m.ref[ids,1:])
            fn,_=vortex(rn*np.sin(a))
            v=t[ids];c=chi[ids];c0=m.ref[ids,0]
            n=(fn*v)**2
            square=lambda x:np.einsum('ij,ij->i',x,x)
            dot=np.einsum('ij,ij->i',gp,gc)
            E=.5*fn**2*square(gv)+fn**10*(v*v-1)**2*(v**6+2*v**4+3*v*v+4)/20
            E+=.5*n*(square(gp)-ORDER_INERTIA*dot**2/(n+ORDER_INERTIA*square(gc)))
            E+=.5*KAPPA*(square(gc)-square(reference))+WELL*((c*(1-c))**2-(c0*(1-c0))**2)+elastic[ids]*(c-c0)
            if np.max(abs(chi-m.ref[:,0])[(z/m.shape.length<.71)&(z>0)])>1e-10:
                raise ArithmeticError('Order correction enters the varied part of the collar')
            energy+=(wt[ids]*(rn/rr)**2*stretch)@E
            face=box[0 if side else -1]
            zz=m.z[1]; az=m.a
            nodes=np.stack([face[i:i+len(az)-1,j:j+nt-1] for i in (0,1) for j in (0,1)],-1)
            for h in (1/2-1/np.sqrt(12),1/2+1/np.sqrt(12)):
                for k in (1/2-1/np.sqrt(12),1/2+1/np.sqrt(12)):
                    a=az[:-1,None]+h*np.diff(az)[:,None];w=zz[None,:-1]+k*np.diff(zz)[None,:]
                    rb=m.shape.radius(w)[0]
                    rn,_,_,_,at=deform(m,rb,w,a,side,step)
                    ff,dh=vortex(rn*np.sin(a))
                    amp=u[0][nodes]@np.kron([1-h,h],[1-k,k])
                    surf=(1-2*side)/2*4*np.pi*rn**2*np.sin(a)*(np.sin(a)-at/rn*np.cos(a))
                    energy+=np.sum(surf*ff*dh*(amp*amp-1)*np.diff(zz)[None,:]*np.diff(az)[:,None]/4)
        return float(energy)
    for h in (.001,.0005):
        plus,minus=moved(h),moved(-h)
        print(json.dumps(dict(kind='joined_shape_work',winding=m.winding,grid=meta['grid'],
                              step=h,force=(minus-plus)/(2*h))),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('state',type=Path)
    parser.add_argument('--material',type=Path)
    parser.add_argument('--grid',type=int,nargs=3)
    arg=parser.parse_args()
    check(arg.state,arg.material,arg.grid)
