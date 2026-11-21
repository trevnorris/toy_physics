"""Finite-offset force test of the degree-one candidate in a spatially localized burst.

Retains degrees 1, 2 and 3 in the deformed guide. Degrees above three first
affect the degree-one frequency beyond second order in the displacement.
"""
import argparse

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

from mouth_native_wave_bend import (Bend, angular_stiffness,assemble,bend_forms,emit,find_pole,
                                   impedance,momentum_and_mass,prepared,truncate)


def run(a):
    dom,chi,intake,_=prepared(a)
    dom,chi,cut=truncate(dom,chi,a.cut)
    points,blocks,_,parts,feed,_=assemble(dom,chi,a.sub,1,a.feed_spacing)
    masses,_=momentum_and_mass(dom,blocks,chi,len(points),0.)
    angular=sum(angular_stiffness(dom,blocks,chi,len(points)))
    outlet,_,_=impedance(dom,blocks,points,chi,cut,0.)
    mass=masses.sum(axis=0)
    labels=(1,2,3)
    ix=[(mass>0)&((points[:,0]>0)|(l==1)) for l in labels]
    scale=[sparse.diags(1/np.sqrt(mass[m])) for m in ix]
    def block(m,i,j):
        return scale[i]@m[ix[i]][:,ix[j]]@scale[j]
    straight=[block(sum(parts)+(l*(l+1)-2)*angular,i,i) for i,l in enumerate(labels)]
    bends=[Bend(bend_forms(dom,blocks[2],chi,len(points),s)) for s in ((.55,1.15),(1.2,2.))]
    f=a.collar_fraction
    def bend_operator(quadratic):
        rows=[]
        for i,l in enumerate(labels):
            row=[]
            for j,k in enumerate(labels):
                forms=[(b.square(l,k) if quadratic else b.linear(l,k)) for b in bends]
                power=2 if quadratic else 1
                mat=None if forms[0] is None else f**power*forms[0]+(1-f)**power*forms[1]
                row.append(block(mat,i,j) if mat is not None else None)
            rows.append(row)
        return sparse.bmat(rows,format='csr')
    first=bend_operator(False)
    second=bend_operator(True)
    kk=sparse.block_diag(straight,format='csr')+a.offset*first+a.offset**2*second
    derivative=first+2*a.offset*second
    damp=np.concatenate([(feed+outlet)[m]/mass[m] for m in ix])
    r,w=points[ix[0]].T
    behind=w>.5000000001
    eig,mode=eigsh(straight[0][behind][:,behind],k=1,sigma=0,v0=np.ones(behind.sum()),tol=1e-12)
    seed=np.zeros(len(r)); seed[behind]=mode[:,0]
    base,_,_=find_pole(straight[0],damp[:len(r)],np.sqrt(eig[0])+1e-4j,seed)
    source=np.zeros(kk.shape[0]); source[:len(r)]=seed
    pole,shape,err=find_pole(kk,damp,base+1e-5j,source)
    ed=(abs(pole)**2+np.vdot(shape,kk@shape).real)/4
    emit(test='offset',intake=intake,offset=a.offset,collar_fraction=f,sub=a.sub,cut=cut,
         omega0=base.real,omega=pole.real,loss=-pole.imag,residual=err,
         frequency_curvature=2*(pole.real-base.real)/(base.real*a.offset**2),
         configurational_flare_force_per_energy=-np.vdot(shape,derivative@shape).real/(4*ed))
    if a.no_burst: return
    drive=3.184
    duration=12*np.pi/drive
    source.fill(0)
    source[:len(r)]=np.sqrt(mass[ix[0]])*np.exp(-(r/.7)**2)*np.maximum(0.,1-((w-.95)/.18)**2)**3
    source/=np.linalg.norm(source)
    forced=np.flatnonzero(source)
    sample=.12
    stride=int(np.ceil(sample*np.sqrt(np.max(np.asarray(abs(kk).sum(axis=1))))/1.45))
    dt=sample/stride
    update=dt/(1+.5*dt*damp)
    decay=(1-.5*dt*damp)/(1+.5*dt*damp)
    u=np.zeros_like(damp); vel=np.zeros_like(damp)
    work=lost=0.
    trace=[]
    for i in range(round(a.end/sample)*stride):
        t=i*dt
        rhs=kk@u
        amp=1e-4*np.sin(np.pi*t/duration)**2*np.cos(drive*(t-duration/2)) if t<duration else 0.
        if amp: rhs[forced]-=amp*source[forced]
        vn=decay*vel-update*rhs
        un=u+dt*vn
        vm=.5*(vel+vn)
        if amp: work+=dt*amp*(source[forced]@vm[forced])
        lost+=dt*(damp@vm**2)
        if (i+1)%stride==0:
            mid=.5*(u+un)
            energy=.5*(vn@vn+u@(kk@un))
            trace.append([(i+.5)*dt-duration,energy,-.5*mid@(derivative@mid),energy+lost-work])
        u,vel=un,vn
    t,energy,force,balance=np.array(trace).T
    error=np.max(abs(balance))/work
    if error>1e-8 or not np.isfinite(trace).all(): raise ArithmeticError('Release ledger failed')
    emit(test='burst',input_work=work,dt=dt,energy_balance=error)
    for lo,hi in ((10,30),(30,54),(54,78)):
        if t[-1]<hi: continue
        mask=(t>=lo)&(t<=hi)
        weight=np.sin(np.pi*(t[mask]-lo)/(hi-lo))**2
        ratio=(weight@force[mask])/(weight@energy[mask])
        emit(test='force',window=[lo,hi],flare_config_force_per_energy=ratio,
             restoring_slope_per_energy=-ratio/a.offset)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.set_defaults(radius=1.,state=None,mesh=34)
    p.add_argument('--offset',type=float,required=True)
    p.add_argument('--collar-fraction',type=float,default=0.)
    p.add_argument('--outer',type=float,default=15.)
    p.add_argument('--end',type=float,default=91.)
    p.add_argument('--sub',type=int,default=1)
    p.add_argument('--cut',type=float,default=2.15)
    p.add_argument('--feed-spacing',type=float,default=.1)
    p.add_argument('--no-burst',action='store_true')
    a=p.parse_args()
    if not a.offset or not 0<=a.collar_fraction<=1: p.error('Specify a finite displacement and a fraction in [0,1].')
    run(a)
