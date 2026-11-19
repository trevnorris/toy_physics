"""Energy and angular impulse for directed, number-conserving conversion.

One parcel visits a single three-site circuit. A forward dispersion v*p
supplies autonomous drift with the same flux velocity for all internal
channels. The calculation uses hbar=1, v=1, incoming kinetic energy 10,
gap 1, and one time unit in the circuit unless specified. The measured
bulk tag accepts the exported excitation; transverse propagation in the
bulk is unnecessary. Link controls have marked orientations.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm_frechet

from mouth_collar_exchange_triangle import ring


H0=np.diag([0.,1.]*3)
RHO_UP=np.diag([1.,0.]*3)/3
FULL=np.eye(6)/6
GAP=np.array([0.,1.]*3)


def adj(x):
    return x.conj().swapaxes(-1,-2)


def tr_product(a,b):
    return np.einsum('...ij,...ji->...',a,b).real


def contact(p,f,j=2.,t=1.):
    """h and dh/dF, with the whole marked twist on directed matrix element 20."""
    p=np.asarray(p)
    base=np.stack([ring(v,j) for v in p.reshape(-1,3,3)]).reshape(p.shape[:-2]+(6,6))
    for site in (0,1):
        base[...,2*site:2*site+2,2*site+2:2*site+4]=-t*np.eye(2)
        base[...,2*site+2:2*site+4,2*site:2*site+2]=-t*np.eye(2)
    shape=np.broadcast_shapes(p.shape[:-2],np.shape(f))
    h=np.broadcast_to(H0+base,shape+(6,6)).copy()
    dh=np.zeros_like(h)
    closing=-t*np.exp(-1j*np.asarray(f))
    for s in (0,1):
        h[...,4+s,s]=closing
        h[...,s,4+s]=closing.conj()
        dh[...,4+s,s]=-1j*closing
        dh[...,s,4+s]=1j*closing.conj()
    return h,dh


def read(p,f,duration=1.,rho=RHO_UP,j=2.,t=1.):
    h,dh=contact(p,f,j,t)
    lam,u=np.linalg.eigh(h)
    phase=np.exp(-1j*duration*lam)
    eig_rho=adj(u)@rho@u
    out=u@(phase[..., :,None]*eig_rho*phase.conj()[...,None,:])@adj(u)
    de=tr_product(H0,out-rho)
    entrance=-tr_product(h-H0,rho)
    exit_work=tr_product(h-H0,out)
    gap=duration*(lam[..., :,None]-lam[...,None,:])
    time_average=duration*np.exp(-.5j*gap)*np.sinc(gap/(2*np.pi))
    impulse=-tr_product(eig_rho*time_average,adj(u)@dh@u)
    return de,impulse,entrance,exit_work


def held_trials():
    triad=np.eye(3)
    for name,p,rho,t in [
        ("ordered",triad,RHO_UP,1.),("mirror",triad[::-1],RHO_UP,1.),
        ("unpolarized",triad,FULL,1.),("no_hopping",triad,RHO_UP,0.),
        ("identical_pins",np.tile([1.,0.,0.],(3,1)),RHO_UP,1.)
    ]:
        for f in (0.,.08,.04,.02):
            de,imp,win,wout=read(p,np.array([f,-f]),rho=rho,t=t)
            H,dH=contact(p,f,t=t)
            S,dS=expm_frechet(-1j*H,-1j*dH)
            check=(-1j*np.trace(rho@adj(S)@dS)).real
            residual=np.max(abs(de+win+wout))
            print(f"held model={name} F={f:+.4f} tag_plus={de[0]:.12g} tag_minus={de[1]:.12g} "
                  f"tag_odd={(de[0]-de[1])/2:.12g} impulse_plus={imp[0]:.12g} "
                  f"impulse_minus={imp[1]:.12g} drift_entry_plus={win[0]:.12g} "
                  f"drift_exit_plus={wout[0]:.12g} energy_residual={residual:.3g} "
                  f"impulse_derivative_error={abs(check-imp[0]):.3g}",flush=True)
            assert residual<1e-12 and abs(check-imp[0])<1e-12
    # A physically added twist-sensitive scalar conversion energy b(F)=b0*sin F
    # through the occupied zone. Its phase is measured with fixed channel frames.
    f=.08
    h,dh=contact(triad,f)
    for b0 in (0.,.15):
        S,dS=expm_frechet(-1j*(h+b0*np.sin(f)*np.eye(6)),
                          -1j*(dh+b0*np.cos(f)*np.eye(6)))
        de=tr_product(H0,S@RHO_UP@adj(S)-RHO_UP)
        imp=(-1j*np.trace(RHO_UP@adj(S)@dS)).real
        print(f"scalar_shift b0={b0:g} tag={de:.12g} impulse={imp:.12g} "
              f"shift_prediction={-b0*np.cos(f):.12g}")


def moving_ledger():
    p=np.eye(3)
    for speed in (.04,-.04,.02,-.02):
        fi=.08
        def ode(t,y):
            rho=y[:-1].reshape(6,6)
            h,dh=contact(p,fi+speed*t)
            return np.r_[(-1j*(h@rho-rho@h)).ravel(), speed*tr_product(rho,dh)]
        y0=np.r_[RHO_UP.astype(complex).ravel(),0.]
        sol=solve_ivp(ode,(0.,1.),y0,method="DOP853",rtol=2e-11,atol=2e-13)
        assert sol.success
        end=sol.y[:,-1]
        out=end[:-1].reshape(6,6)
        hi,_=contact(p,fi)
        hf,_=contact(p,fi+speed)
        drift=-tr_product(hi-H0,RHO_UP)+tr_product(hf-H0,out)
        tag=tr_product(H0,out-RHO_UP)
        actuator=end[-1].real
        # Force on the twist is the negative of the work supplied by its driver.
        residual=drift+tag-actuator
        print(f"moving rate={speed:+.3f} tag={tag:.12g} drift={drift:.12g} "
              f"work_on_parcel={actuator:.12g} work_to_twist={-actuator:.12g} "
              f"avg_impulse={-actuator/speed:.12g} residual={residual:.3g}",flush=True)
        assert abs(residual)<3e-10


if __name__=="__main__":
    held_trials()
    moving_ledger()
