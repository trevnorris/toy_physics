import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.linalg import eigh
def mode(R):
    p=30/(19*R)
    def rhs(s,z):
        b,w,a=z
        return np.array([np.sin(a),np.cos(a),2*np.cos(a)/b-p])
    def end(s,z): return z[2]
    end.direction=-1;end.terminal=True
    y=solve_ivp(rhs,[0,12*R],[R,0,0],events=end,rtol=1.5e-12,atol=1e-13,dense_output=True)
    S=y.t_events[0][0];s=np.linspace(0,S,5001); b,w,alpha=y.sol(s)
    weight=4*np.pi*b*b; curvature=2*(np.cos(alpha)/b)**2+rhs(s,y.sol(s))[-1]**2
    for n in (1,2,4,8,12):
        k=np.arange(1,n+1)[:,None]*np.pi/S
        f=np.cos(k*s)
        d=-k*np.sin(k*s)
        f-=(simpson(weight*f,x=s)/simpson(weight,x=s))[:,None]
        mass=simpson(weight*f[:,None]*f[None,:],x=s)
        op=simpson(weight*(d[:,None]*d[None,:]-curvature*f[:,None]*f[None,:]),x=s)
        lam,z=eigh(op,mass)
        best=z[:,0]@f
        dV=simpson(best*weight,x=s)
        assert abs(dV)<1e-11
        print('arc',R,n,lam[:2],dV)
    k=np.arange(1,5)[:,None]*np.pi/w[-1]
    f=np.cos(k*w);d=-k*np.sin(k*w)*np.cos(alpha)
    f-=(simpson(weight*f,x=s)/simpson(weight,x=s))[:,None]
    mass=simpson(weight*f[:,None]*f[None,:],x=s)
    op=simpson(weight*(d[:,None]*d[None,:]-curvature*f[:,None]*f[None,:]),x=s)
    print('axial',R,eigh(op,mass,eigvals_only=True)[:2])
for R in (1,2):mode(R)
