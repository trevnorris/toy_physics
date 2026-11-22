"""Independent P1 crystal amplitude and finite-volume two-coordinate phase conduction."""
import numpy as np
from scipy.sparse import coo_matrix,diags
from scipy.sparse.linalg import spsolve
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import interp1d

def c(r):
    return 100*np.sin(r/100)
def jac(r):
    return c(r)*np.cos(r/100)

def amplitudes(n):
    r=np.linspace(0,8,n+1); h=8/n
    x,w=leggauss(3)
    t=(x+1)/2; w=w*h/2
    pts=r[:-1,None]+h*t; N=np.array([1-t,t]); G=np.array([-1/h,1/h])
    wg=w*jac(pts); f=np.tanh(r*.6); f[-1]=1
    for it in range(14):
        fv=np.column_stack([f[:-1],f[1:]])@N
        dp,d2p=fv/c(pts)**2+fv*(fv**2-1), 1/c(pts)**2+3*fv**2-1
        deriv=np.diff(f)/h
        z=np.zeros_like(f); diag=np.zeros_like(f); off=np.zeros(n)
        for j in [0,1]:
            sl=slice(j,n+j)
            z[sl]+=np.sum(wg*(G[j]*deriv[:,None]+N[j]*dp),axis=1)
            diag[sl]+=np.sum(wg*(G[j]**2+N[j]**2*d2p),axis=1)
        off[:]=np.sum(wg*(G[0]*G[1]+N[0]*N[1]*d2p),axis=1)
        err=max(abs(z[1:-1]))
        if err<2e-12: break
        mat=diags([off[1:-1],diag[1:-1],off[1:-1]],[-1,0,1],format='csc')
        f[1:-1]-=spsolve(mat,z[1:-1])
    assert err<2e-12
    E=2*np.pi*np.sum(wg*(.5*deriv[:,None]**2+.5*(fv/c(pts))**2+.25*(1-fv**2)**2))
    return interp1d(r,f),E,np.interp(.5,f,r)

def phase(nr,ns):
    f,E,half=amplitudes(nr*2)
    dr,ds=8/nr,32/ns
    r=(np.arange(nr)+.5)*dr
    s=(np.arange(ns)+.5)*ds-16
    rface=(np.arange(1,nr))*dr
    area=jac(r)*dr
    ids=np.arange(nr*ns).reshape(nr,ns)
    sc=np.r_[-16,s,16]; lo,hi=sc[:-1],sc[1:]
    x,w=leggauss(8)
    probe=(lo+hi)[:,None]/2+(hi-lo)[:,None]*x/2
    ring=f(r)**2
    normal=np.empty((nr,ns+1))
    for i,a in enumerate(ring):
        resist=((1+probe/100)**-3/(1-a/np.cosh(probe/2)**2))@w*(hi-lo)/2
        normal[i]=area[i]/resist
    lateral=jac(rface[:,None])*(1+s/100)*ds/dr
    a=np.r_[ids[:,:-1].ravel(),ids[:-1,:].ravel()]; b=np.r_[ids[:,1:].ravel(),ids[1:,:].ravel()]
    k=np.r_[normal[:,1:-1].ravel(),lateral.ravel()]
    bc=np.r_[ids[:,0],ids[:,-1]]
    M=coo_matrix((np.r_[k,k,-k,-k,normal[:,0],normal[:,-1]],(np.r_[a,b,a,b,bc],np.r_[a,b,b,a,bc])),shape=(ids.size,ids.size)).tocsr()
    rhs=np.zeros(ids.size);rhs[ids[:,-1]]=normal[:,-1]
    theta=spsolve(M,rhs).reshape(nr,ns)
    q0=np.sum(normal[:,0]*theta[:,0]);q1=np.sum(normal[:,-1]*(1-theta[:,-1]))
    total=np.sum(area); r0=50*((1-16/100)**-2-(1+16/100)**-2)
    lower=sum(1./np.sum(1/normal,axis=1))/total*r0
    print(dict(n=nr,depth_mesh=ns,crystal_energy=E,half_amplitude_radius=half,reduced_conductance=q0/total*r0,relative_flux_disagreement=(q0-q1)/q0,parallel_path_bound=lower))
    assert abs(q0-q1)<2e-10 and q0/total*r0>=lower-1e-7
for nr,ns in [(48,128),(96,256)]:
    phase(nr,ns)
