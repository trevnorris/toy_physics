import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit

k=2.
def test(N,shift=0):
    x=np.linspace(-16,16,N+1); h=x[1]-x[0]
    phase=.4510084282857874
    def data(y):
        f=np.r_[1,y[:N-1],1]
        q=np.r_[1,y[N-1:],0]
        F=(f[:-1]+f[1:])/2; Q=(q[:-1]+q[1:])/2
        df=np.diff(f)/h; dq=np.diff(q)/h
        return f,q,F,Q,df,dq
    def energy(y):
        f,q,F,Q,df,dq=data(y); n=F**2
        V=(n-1)**2*(n**3+2*n*n+3*n+4)/20
        work=np.sum(.5*df*df+V+k*.5*dq*dq+Q**2*(1-Q)**2)*h
        cap=h*np.sum(1/n+dq**2/n**2); i=phase/cap
        HF=h*(2*F*((n-1)*(n+1)*(n**2+1)/4) + i*i*(1/F**3+2*dq**2/F**5))/2
        Hq=h*Q*(1-Q)*(1-2*Q)
        Dw=(k-i**2/F**4)*dq
        gF=(HF[:-1]+df[:-1])+(HF[1:]-df[1:])
        gq=(Hq[:-1]+Dw[:-1])+(Hq[1:]-Dw[1:])
        return work+.5*phase*i, np.r_[gF,gq]
    seed=np.r_[np.ones(N-1),expit(-x[1:-1]+shift)]
    opt=minimize(energy,seed,jac=True,method='L-BFGS-B',options={'maxiter':20000,'gtol':6e-10,'ftol':2e-16,'maxls':40})
    f,q,F,Q,df,dq=data(opt.x)
    E,g=energy(opt.x);n=F*F
    i=phase/(h*np.sum(1/n+dq**2/n**2))
    front=np.interp(.5,q[::-1],x[::-1])
    at_mid=np.interp(front,.5*(x[:-1]+x[1:]),dq)
    print(dict(cells=N,seed_offset=shift,minimum_location=front,relative_current=i,lowest_density=min(f)**2,order_slope=at_mid,stationarity=max(abs(g)),energy=E))
    assert max(abs(g))<6e-7 and abs(i-.01401962575665)<2e-7
    return E
for cells in [256,512]:
    E=test(cells); F=test(cells,1.75)
    print('free-translation energy change:',abs(E-F))
    assert abs(E-F)<2e-8
