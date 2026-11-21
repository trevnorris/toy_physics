"""Cross-check wide finite-bias cases in physical coordinates by stationary action.

Convergence uses the size of the Newton correction, nodal-volume normalized
equation residual, agreement of the two independently integrated bath currents
and a centered on-shell energy derivative with respect to the imposed phase.
"""
import argparse
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from mouth_carrier_lateral_response import Dipole, emit, weighted_gram
from mouth_carrier_lateral_inertia import WINDOWS, volume
from mouth_carrier_weak_radius import PhysicalDomain


def at_bias(dom,phase,bias,guess):
    K=dom.stiff*.5
    free=dom.free
    basis=dom.N[:,free]
    k=K[free][:,free]
    # The diagonal of the consistent mass matrix gives a positive local scale.
    weight=basis.power(2).T@dom.weight
    if weight.min()<=0:
        raise ValueError('Diagonal nodal volume has an invalid sign')
    z=guess.copy() if guess is not None else np.expm1(1j*bias*phase)
    z[dom.bulk]=np.expm1(1j*bias)
    z[dom.feed]=0.
    for step in range(14):
        u=dom.N@z
        dn=2*u.real+abs(u)**2
        n=1+dn
        f=1+u
        pot=dom.weight*dn*(n+1)*(n*n+1)/4
        residual=K@z+dom.N.T@(pot*f)
        derivative=2*dom.weight*n**3
        xy=weighted_gram(basis,derivative*f.real*f.imag)
        mat=sparse.bmat([[k+weighted_gram(basis,pot+derivative*f.real**2),xy],
                         [xy,k+weighted_gram(basis,pot+derivative*f.imag**2)]],format='csc')
        delta=spsolve(mat,-np.r_[residual[free].real,residual[free].imag])
        err=max(abs(residual[free])/weight)
        correction=max(abs(delta))
        if err<2e-8 and correction<2e-12:
            break
        re,im=np.split(delta,2)
        z[free]+=re+1j*im
    else:
        raise ArithmeticError(f'Unsettled stationary action: {err=}, {correction=}')
    field=1+z
    current=2*(field.conj()*residual).imag
    left,right=-sum(current[dom.feed]),sum(current[dom.bulk])
    kinetic=np.vdot(z,K@z).real
    energy=kinetic+dom.weight@(dn**2*(n**3+2*n**2+3*n+4)/20)
    emit(test='stationary_action',bias=bias,feed=left,bulk=right,energy=energy,
         max_correction=correction,max_volume_normalized_residual=err,raw_residual=max(abs(residual[free])),
         iterations=step,min_density=min(n))
    if abs(left-right)>2e-8*max(abs(left),.01) or min(n)<.5:
        raise ArithmeticError('Signed current closure or density test failed')
    return z,(left+right)/2,energy


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('radius','depth','outer','roof','bias'):
        p.add_argument('--'+key,type=float,required=True)
    p.add_argument('--mesh',type=int,required=True)
    a=p.parse_args()
    dom=PhysicalDomain(a.radius,a.depth,a.mesh,a.outer,a.roof,1.)
    emit(test='geometry',**vars(a))
    phase=dom.bulk.astype(float)
    mask=dom.free
    phase[mask]=spsolve(dom.stiff[mask][:,mask],-(dom.stiff@phase)[mask])
    halfstep=.0005
    z,rate,energy=at_bias(dom,phase,a.bias,None)
    _,_,minus=at_bias(dom,phase,a.bias-halfstep,z)
    _,_,plus=at_bias(dom,phase,a.bias+halfstep,z)
    derivative=(plus-minus)/(2*halfstep)
    emit(test='energy_derivative',finite_difference=derivative,flux=rate,halfstep=halfstep)
    if abs(derivative-rate)>2e-6*rate:
        raise ArithmeticError('Energy derivative failed independent current comparison')
    field=1+z
    model=Dipole(dom,field)
    feed,_=model.steady(1.,0.)
    bulk,_=model.steady(0.,1.)
    tow,_=model.steady(-1.,-1.)
    src=dom.N.T@(dom.weight*(1j*(dom.N@tow)-model.rq*(dom.N@field)))
    h,_=model.solve(np.zeros_like(field),src)
    for w in np.array(WINDOWS)*[a.radius,a.radius,a.depth,a.depth]:
        if w[1]>=a.outer or w[3]>=a.roof:
            raise ValueError('Interior window intersects a driven bath')
        C,B,T=[volume(dom,np.c_[field,u],w,False)['virtual_wall_force'] for u in (feed,bulk,tow)]
        M=-volume(dom,np.c_[field,h,tow],w,True)['virtual_wall_force']
        emit(test='matched_radius',radius=a.radius,depth=a.depth,mesh=a.mesh,window=w.tolist(),
             rate=rate,C=C,B=B,T=T,M=M,ratio=C/M)
        if abs(C+B+T)>5e-8 or M<=0:
            raise ArithmeticError('Independent tow or mass check failed')
