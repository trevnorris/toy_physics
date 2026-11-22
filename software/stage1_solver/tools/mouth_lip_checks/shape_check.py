"""Independent P1 carrier/backing, closed cavity, and physical frusta (local, radius 1 and 2 only)."""
import argparse, json, sys
from functools import lru_cache
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp, simpson
from scipy.interpolate import CubicSpline
sys.path.insert(0,'software/stage1_solver/tools')
sys.path.insert(0,'software/stage1_solver/tools/mouth_exchange_checks')
import duct, wave
from mouth_exchange_geometry import ExchangeShape

@lru_cache(1)
def shape0():
 def ode(s,y):
    b,w,t=y
    return np.sin(t),np.cos(t),2*np.cos(t)/b-30/19
 def end(s,y): return y[2]
 end.direction=-1; end.terminal=True
 sol=solve_ivp(ode,(0.,10.),(1.,0.,0.),events=end,rtol=1e-12,atol=1e-13,dense_output=True)
 t=np.linspace(0.,sol.t_events[0][0],1001)
 b,w,a=sol.sol(t)
 return w[-1],CubicSpline(w/w[-1],b,bc_type=((1,0.),(1,0.)))

class Shape:
 def __init__(self,R,x):
    self.neck=R;self.x=x
    L,self.b=shape0()
    self.length=R*L*np.exp(x[0])
    self.exit=1.5*R*np.exp(x[1])
 def radius(self,w):
    t=np.clip(w/self.length,0.,1.); e=self.exit/self.neck-1.5; m,a=self.x[2:]
    p=self.b(t)+e*t+4*m*t*(1-t)+4*a*t*(1-t)*(2*t-1)
    dp=self.b(t,1)+e+4*m*(1-2*t)+4*a*(-6*t*t+6*t-1)
    return p*self.neck,np.where((w<0)|(w>self.length),0.,dp*self.neck/self.length)

def fields(sh,n):
 z=np.linspace(0.,sh.length,4001);r,s=sh.radius(z)
 mer=-np.gradient(s,z,edge_order=2)/(1+s*s)**1.5
 hoop=1/(r*np.hypot(1.,s))
 metric=4*np.pi*r**2*np.hypot(1.,s)
 # Area with the two rotational principal curvatures appropriate to four-space.
 wall=simpson(metric,x=z)+4*np.pi*sh.exit**3/3
 corr=simpson(metric*.04**2*(-5/12*(mer+2*hoop)**2+(np.pi**2/3-2)*(2*hoop*mer+hoop**2)),x=z)
 return np.array([wall,corr,duct.solve(sh,n),wave.clamped(sh,n)[0]])

def run():
 p=argparse.ArgumentParser(); p.add_argument('kind',choices=['side','axial']);p.add_argument('R',type=int,choices=(1,2));p.add_argument('n',type=int);p.add_argument('--x',nargs=4,type=float,default=[0,0,0,0]);p.add_argument('--direction',nargs=4,type=float,default=[1,0,0,0]);v=p.parse_args()
 if v.kind=='side':
    assert v.R==1
    h=.005
    E=[duct.solve(ExchangeShape(q),v.n) for q in (-h,0,h)]
    print(json.dumps(dict(n=v.n,h=h,E=E,terminal_linear_current=2*E[1]/.4510084282857874,slope=(E[-1]-E[0])/(2*h),stiffness=(E[-1]+E[0]-2*E[1])/h**2)))
    return
 h=.0025; x=np.array(v.x); d=np.array(v.direction)
 F=np.array([fields(Shape(v.R,x+q*d),v.n) for q in (-h,0,h)])
 g=(F[-1]-F[0])/(2*h)
 print(json.dumps(dict(R=v.R,n=v.n,x=x.tolist(),direction=d.tolist(),labels=['area_including_cap','order_correction','sealed_weak_carrier','closed_wave'],reference=F[1].tolist(),gradient=g.tolist(),total_gradient=g.sum())))
if __name__=='__main__':run()
