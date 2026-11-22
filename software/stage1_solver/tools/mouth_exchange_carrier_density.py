"""Density relaxation of the coherent sharp-front carrier at finite drive.

The amplitude is continuous across the reversible thin layer; the phase jump
energy scales as the square of the local density. Cell-average quadrature is
used for the density multiplying the Q1 phase energy and for the potential.
"""
import argparse
import json
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from mouth_exchange_carrier import assemble, BIAS
from mouth_exchange_geometry import ExchangeShape
from mouth_exchange_interior import InteriorShape


def solve(step, n, bias, geometry):
    pieces,volumes,_,stiff,feed,bulk,free=assemble(step,n,.04/6,28.,22.,geometry,return_system=True)
    index=np.arange(stiff.shape[0])
    for name,ids,_,_ in pieces:
        if name.endswith('layer'):
            index[ids[:,2:].ravel()]=ids[:,:2].ravel()
    _, density_ids = np.unique(index,return_inverse=True)
    density_count=max(density_ids)+1
    density_free=np.ones(density_count,bool)
    density_free[density_ids[np.r_[feed,bulk]]]=False
    phase=np.zeros(len(index));phase[bulk]=bias
    f=np.ones(density_count)
    cache=[]
    for name,ids,K,dK in pieces:
        # Across an element only phase differences enter the local quadratic form.
        row=np.broadcast_to(ids[:,:,None],K.shape).ravel()
        col=np.broadcast_to(ids[:,None,:],K.shape).ravel()
        did=density_ids[ids]
        drow,dcol=density_ids[row],density_ids[col]
        cache.append((name,ids,K,dK,row,col,did,drow,dcol))
    rows=np.concatenate([a[4] for a in cache])
    cols=np.concatenate([a[5] for a in cache])
    dr,dc=density_ids[rows],density_ids[cols]
    ones=np.ones((4,4))/16
    previous=np.inf
    for it in range(50):
        pvals=[]
        for name,ids,K,_,_,_,did,_,_ in cache:
            exponent=4 if name.endswith('layer') else 2
            avg=f[did].mean(axis=1)
            pvals.append((avg[:,None,None]**exponent*K).ravel())
        A=sparse.coo_matrix((np.concatenate(pvals),(rows,cols)),shape=stiff.shape).tocsr()
        phase[free]=0.
        phase[free]=spsolve(A[free][:,free],-(A@phase)[free])
        grad=np.zeros_like(f)
        fvals=[]
        energy={}
        work={}
        flux={}
        for name,ids,K,dK,_,_,did,_,_ in cache:
            local=phase[ids]-phase[ids][:,:1]
            avg=f[did].mean(axis=1)
            stored=.5*np.einsum('ni,nij,nj->n',local,K,local)
            W=.5*np.einsum('ni,nij,nj->n',local,dK,local)
            layer=name.endswith('layer')
            k=4 if layer else 2
            e=avg**k*stored
            v=avg**k*W
            g=k*avg**(k-1)*stored
            H=k*(k-1)*avg**(k-2)*stored
            if not layer:
                v0,dv=volumes[name]
                localf=f[did]-f[did][:,:1]
                rho=avg**2
                potential=(rho-1)**2*(rho**3+2*rho**2+3*rho+4)/20
                e+=v0*potential+.5*np.einsum('ni,nij,nj->n',localf,K,localf)
                v+=dv*potential+.5*np.einsum('ni,nij,nj->n',localf,dK,localf)
                g+=v0*avg*(rho**4-1)/2
                H+=v0*(9*rho**4-1)/2
                np.add.at(grad,did,np.einsum('nij,nj->ni',K,localf))
                hess=K+H[:,None,None]*ones
            else:
                flux[name]=-float((avg[:,None]**4*np.einsum('nij,nj->ni',K[:,:2,:],local)).sum())
                hess=H[:,None,None]*ones
            np.add.at(grad,did,np.broadcast_to(g[:,None]/4,did.shape))
            energy[name]=sum(e);work[name]=sum(v)
            fvals.append(hess.ravel())
        H=sparse.coo_matrix((np.concatenate(fvals),(dr,dc)),shape=(density_count,density_count)).tocsr()
        change=spsolve(H[density_free][:,density_free],-grad[density_free])
        max_change=np.max(abs(change)); value=sum(energy.values())
        if max_change<1e-10:
            break
        if value>previous+3e-9 or min(f)<.6:
            raise ArithmeticError('Carrier density left near-unit lowering branch')
        f[density_free]+=change
        previous=value
    else:
        raise ArithmeticError('Density and phase have not reached a common stationary field')
    currents=A@phase
    D=sum(currents[bulk])
    if max(abs(sum(currents[feed])+D),abs(sum(flux.values())-D))>2e-7*abs(D):
        raise ArithmeticError('Finite-drive current fails two-sided closure')
    return dict(step=step,mesh=n,bias=bias,energy=value,work=sum(work.values()),
                sector_work=work,sector_energy=energy,current=D,current_by_cut=flux,
                minimum_amplitude=min(f),iterations=it+1,density_correction=max_change)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--mesh',type=int,required=True)
    p.add_argument('--step',type=float,required=True)
    p.add_argument('--interior',action='store_true')
    p.add_argument('--bias-fraction',type=float,default=1.)
    a=p.parse_args()
    print(json.dumps(solve(a.step,a.mesh,BIAS*a.bias_fraction,
                           InteriorShape if a.interior else ExchangeShape),allow_nan=False))
