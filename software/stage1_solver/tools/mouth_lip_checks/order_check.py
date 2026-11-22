"""Independent P1 minimized curved patch and analytic response at the symmetric flat front."""
import json
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from scipy.integrate import quad
from scipy.special import expit

def energy(momentum, count, direction):
    # Normal distance in units of delta_B, energy per middle area in a_B delta_B.
    x = np.linspace(-16,16,count+1)
    h = np.diff(x)
    q,w = np.polynomial.legendre.leggauss(5); q=(q+1)/2; w=w/2
    s=x[:-1,None]+q*h[:,None]
    delta=np.asarray(direction)[None,None,:]*s[:,:,None]
    area=np.prod(1+momentum*delta, axis=-1)
    darea=area*np.sum(delta/(1+momentum*delta),axis=-1)
    v=expit(-x); v[0]=1; v[-1]=0
    live=np.setdiff1d(np.arange(len(v)),[0,count//2,count])
    weight=h[:,None]*w[None,:]*area
    for _ in range(12):
        slope=np.diff(v)/h
        c=v[:-1,None]*(1-q)+v[1:,None]*q
        wp=2*c*(1-c)*(1-2*c); wpp=2*(1-6*c+6*c*c)
        dg=(weight*2*slope[:,None]/h[:,None]).sum(axis=1)
        g=np.r_[(weight*wp*(1-q)).sum(axis=1)-dg,0.] + np.r_[0.,(weight*wp*q).sum(axis=1)+dg]
        L=(weight*2/(h*h)[:,None]).sum(axis=1)
        diag=np.r_[(weight*wpp*(1-q)**2).sum(axis=1)+L,0.]+np.r_[0.,(weight*wpp*q*q).sum(axis=1)+L]
        off=(weight*wpp*q*(1-q)).sum(axis=1)-L
        A=diags([off,diag,off],[-1,0,1],format='csr')[live,:][:,live]
        if max(abs(g[live]))<1e-12: break
        v[live]-=spsolve(A,g[live])
    assert max(abs(g[live]))<1e-12,'P1 patch stationarity'
    e=slope[:,None]**2+(c*(1-c))**2
    return (weight*e).sum(), (h[:,None]*w[None,:]*darea*e).sum()

def main():
    # With y=chi_0', write the linear order correction as y h. h(0)=0,
    # h' = integral_s^infty y^2 / y^2 on the positive half. The second
    # variation at fixed pin supplies the negative relaxation coefficient.
    alpha=-2*quad(lambda t:t*(.5-t/3)**2/(1-t)**3,0,.5,epsabs=1e-13)[0]
    # The mixed-curvature metric term comes directly from the undeformed kink.
    beta=4*quad(lambda s:s*s*(expit(-s)*expit(s))**2,0,60,epsabs=1e-13)[0]
    print(json.dumps(dict(analytic_H2=alpha,analytic_pair=beta,flat_tension=1/3)),flush=True)
    for n in (512,1024):
        for direction in ((1,0,0),(1,-1,0),(1,1,1),(-.3,.7,.7)):
            d=np.array(direction); H=d.sum(); P=(H**2-d@d)/2
            expected=2*(alpha*H**2 + beta*P)
            for t in (.01,.005):
                vals=[energy(z,n,direction) for z in (-t,0,t)]
                E,M=np.array(vals).T
                slope=(E[-1]-E[0])/(2*t)
                stiffness=(M[-1]-M[0])/(2*t)
                print(json.dumps(dict(n=n,d=list(direction),delta_curvature=t,energies=E.tolist(),moments=M.tolist(),slope=slope,stiffness_from_energy=(E[0]+E[-1]-2*E[1])/t**2,stiffness_from_moment=stiffness,predicted=expected)),flush=True)
                assert abs(slope)<1e-10 and abs(M[1])<1e-10
                assert abs(stiffness-expected)<max(.0004*abs(expected),.0001) if (t==.005 and n==1024) else True

if __name__=='__main__': main()
