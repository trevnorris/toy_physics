"""Reintegrate fixed-wall angular energies and stretches in four dimensions."""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import RegularGridInterpolator

from mouth_bonded_front_adjoint import get_material_load
from mouth_bonded_material_motion import Bonded
from mouth_joined_front_audit import Front
from mouth_joint_polar_shape import Trial
from mouth_polar_exit_readout import accepted
from mouth_supercore_front import check_guard
from mouth_supercore_joined_front import emit,file_sha


def angular_material(body,d,e,s,report=False):
    sh=body.chart.shape
    local=d.reshape(-1,2)[body.id]
    pos=np.einsum('ij,ijk->ik',body.basis,local)
    r,z=body.r+pos[:,0], body.w+pos[:,1]
    b,bz=sh.radius(z)
    x=r/b
    v=(z/sh.length-.06)/.56
    win=(0<v)&(v<1)
    h=np.where(win,16*v**2*(1-v)**2,0.)
    hw=np.where(win,32*v*(1-v)*(1-2*v)/(.56*sh.length),0.)
    B=np.eye(2)+np.stack([np.einsum('ij,ijk->ik',grad,local) for grad in (body.dr,body.dw)],axis=-1)
    K=B@body.K@B.transpose(0,2,1)
    mer=np.linalg.cholesky(K)
    root=np.zeros((len(r),4,4))
    root[:,:2,:2]=mer
    root[:,2,2]=root[:,3,3]=body.h*r/body.r
    if np.min(body.j*np.linalg.det(B))<=0:
        raise ValueError('Inverted meridional material chart')
    interp=RegularGridInterpolator(tuple(np.linspace(0,1,n) for n in body.grid.shape),body.grid)
    g=x**3*(e+s*(1-x))
    gx=x*x*(3*e+3*s-4*s*x)
    mu,weight=leggauss(10)
    mu,weight=(mu+1)/2,weight/2
    if report:
        mu=np.r_[mu,0.,1.];weight=np.r_[weight,0.,0.]
    energy=0.; minimum=minimum_tip=np.inf
    maximum=0.; jmin=np.inf
    for a,w in zip(mu,weight):
        p=(3*a*a-1)/2
        physical=r+g*h*p
        xs=x.copy()
        for _ in range(6):
            residual=b*xs+e*xs**3*h*p-physical
            xs-=residual/(b+3*e*xs**2*h*p)
        if max(abs(b*xs+e*xs**3*h*p-physical))>3e-13:
            raise ValueError('Order chart inversion failed')
        chi=interp((z/sh.length,xs))
        pull=np.zeros_like(root)
        pull[:,0,0]=1+h*gx*p/b
        pull[:,0,1]=p*(g*hw-gx*h*x*bz/b)
        pull[:,0,2]=g*h*(-3*a*np.sqrt(1-a*a))/r
        pull[:,1,1]=1
        pull[:,2,2]=pull[:,3,3]=physical/r
        F=pull@root
        J=np.linalg.det(F)
        # Native solid energy is per unloaded crystalline reference volume.
        W=(np.sum(F**2,axis=(1,2))-4)/2 - np.log(J) + (J-1)**2/4
        energy += w*(body.measure@(chi*W))
        if report:
            sv=np.linalg.svd(F,compute_uv=False)
            minimum=min(minimum,sv[:,-1].min())
            maximum=max(maximum,sv[:,0].max())
            minimum_tip=min(minimum_tip,sv[(z>0.95*sh.length)&(x>.9),-1].min())
            jmin=min(jmin,J.min())
    return energy,dict(minimum=minimum,maximum=maximum,minimum_tip=minimum_tip,Jmin=jmin)


def run(args):
    check_guard(args.endpoint)
    prep,row,x,d0,u0=accepted(args.state,args.snapshot)
    base=Front(args.state,tuple(prep['field_grid']),1)
    n,k=(d0.shape[1]-1)//3,d0.shape[0]-1
    trial=Trial(base,args.state,x,n,k,False)
    with np.load(args.endpoint) as f:
        if not np.array_equal(x,f['x']) or str(f['input_sha'])!=row['sha256'] or str(f['state_sha'])!=file_sha(args.state):
            raise ValueError('Mismatch between angular shape and native material preparation')
        d1=f['displacement'];u1=f['field'];s=float(f['internal_amplitude'])
        recorded=json.loads(str(f['readout']))
    body=Bonded(args.state,n,k,shape=trial.shape,full_depth=True)
    m=trial.model
    result=[]
    for label,d,u,amp in [('initial',d0,u0,0.),('relaxed',d1,u1,s)]:
        solid,limits=angular_material(body,d,x[5],amp,True)
        m.material_load=get_material_load(body,d[body.mask]*body.scale,m)
        result.append(solid+m.evaluate(u)[0])
        emit(kind='material_review',shape=row['sha256'],source=file_sha(args.endpoint),
             state=label,s=s if label=='relaxed' else 0.,solid=solid,total=result[-1]+trial.static,**limits)
        for step in (4e-5,2e-5):
            lo=angular_material(body,d,x[5],amp-step)[0]
            hi=angular_material(body,d,x[5],amp+step)[0]
            emit(kind='signed_angular',state=label,step=step,gradient=(hi-lo)/(2*step),
                 curvature=(hi+lo-2*solid)/step**2)
    drop=result[0]-result[1]
    emit(kind='energy_drop',value=drop,native=recorded)
    if abs((result[0]+trial.static)-row['energy'])>1e-7:
        raise ArithmeticError('Independent initial energy disagrees with accepted value')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('state','snapshot','endpoint'):
        p.add_argument(key,type=Path)
    run(p.parse_args())
