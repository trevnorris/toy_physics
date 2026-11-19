"""Test a position-dependent conversion energy with no direct twist bias.

The prior preparation gives mean deposited quanta e_i=(1-r^6)(1-n_i,z)/2
after six encounters with a local axis n_i. A proposed static energy shift
beta*e_i on the transiting parcel is spin independent and switches off at
both ports. A closed-cycle test is independent of common energy offsets.
"""
import numpy as np
from scipy.linalg import expm_frechet
from scipy.spatial.transform import Rotation

from mouth_collar_order_preparation import calibration
from mouth_conversion_shared_links import prepared,RATE
from mouth_conversion_twist import contact,adj,RHO_UP,H0,tr_product
from mouth_material_links import section


def force_and_tag(p,shift,f):
    h,dh=contact(p,f)
    h[...,np.arange(6),np.arange(6)] += np.repeat(shift,2,axis=-1)
    vals,u=np.linalg.eigh(h)
    r=adj(u)@RHO_UP@u
    phase=np.exp(-1j*vals)
    out=u@(r*phase[..., :,None]*phase.conj()[...,None,:])@adj(u)
    gap=vals[..., :,None]-vals[...,None,:]
    rmean=r*np.exp(-.5j*gap)*np.sinc(gap/(2*np.pi))
    return -tr_product(rmean,adj(u)@dh@u),tr_product(H0,out-RHO_UP)


def phases(p,shift,angles):
    impulse=0.
    tag=0.
    for k in range(3):
        j,e=force_and_tag(np.roll(p,k,axis=1),np.roll(shift,k,axis=1),angles)
        impulse+=j/3
        tag+=e/3
    return impulse,tag


def period(p,b,steps):
    phase=(np.arange(steps)+.5)*2*np.pi/steps
    # Batch along phase to keep peak allocation small on the finer sphere.
    total=np.zeros(len(b))
    for segment in np.array_split(phase,steps//8):
        j,_=phases(p,b,segment[:,None])
        total+=j.sum(axis=0)
    return total*(2*np.pi/steps)


def direct_cycle(p,b,n=48):
    # Independent Frechet evaluation on the two faces of a selected edge.
    value=np.zeros(len(b))
    for k in range(3):
        pk,bk=np.roll(p,k,axis=1),np.roll(b,k,axis=1)
        for f in (np.arange(n)+.5)*2*np.pi/n:
            h,dh=contact(pk,f)
            for i in range(len(b)):
                s,ds=expm_frechet(-1j*(h[i]+np.diag(np.repeat(bk[i],2))),-1j*dh[i])
                value[i]+=(-1j*np.trace(RHO_UP@adj(s)@ds)).real
    return value*2*np.pi/(3*n)


def run(level,r):
    pair,w,near,sign,e=prepared(level,r)
    normals,triangles=section(level)
    normals=Rotation.from_rotvec([.26,-.19,.17]).apply(normals)
    emitted=(1-r**6)*(1-normals[triangles,2])/2
    # Independent material observable: trend in tagged spin-independent
    # site energy with that site's prior released excitation density.
    for beta in (0.,.5,1.,-1.):
        cycles=[]
        for eta,p,b in zip((1,-1),pair,(beta*emitted,beta*emitted[:,::-1])):
            mean=period(p,b,32)
            rate=RATE*w
            edge_work=np.sum((rate*mean)[near]*sign,axis=1)
            i=np.argmax(abs(edge_work))
            direct=direct_cycle(p[near[i]],b[near[i]]) if beta==1 and eta==1 else None
            ferr=0. if direct is None else np.max(abs(direct-mean[near[i]]))
            idx=near[e]
            _,high=phases(p[idx],b[idx],sign[e]*.08)
            _,low=phases(p[idx],b[idx],-sign[e]*.08)
            odd=.5*rate[idx]@(high-low)
            print(f"shift level={level} eta={eta:+} beta={beta:+g} min_face_period={min(mean):.12g} "
                  f"max_face_period={max(mean):.12g} max_shared_cycle_abs={abs(edge_work[i]):.12g} "
                  f"argmax_edge={i} test_edge_work={edge_work[e]:.12g} test_odd_tag={odd:.12g} "
                  f"max_weighted_single_return={max(abs(rate*mean)):.12g} "
                  f"independent_face_error={ferr if direct is not None else 'not_evaluated'}",flush=True)
            assert ferr<1e-10
            cycles.append(edge_work)
        err=np.max(abs(cycles[0]+cycles[1]))
        assert err<1e-12
        if beta==0: assert max(abs(cycles[0]))<1e-13


if __name__=="__main__":
    r=calibration()
    for level in (2,3):
        run(level,r)
