"""Independent meridional circumference/frustum audit against the archived 10-cosine direction."""
from pathlib import Path
import hashlib, json, sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.interpolate import CubicSpline
sys.path.insert(0,'software/stage1_solver/tools')
from mouth_exchange_geometry import ExchangeShape
archive=json.loads(Path('research/pde_ledger_v3/_measurements/mouth_conversion_shape.json').read_text())['primary']['exchange']
assert hashlib.sha256(archive['literal_stdout'].encode()).hexdigest()==archive['stdout_sha256']
ref,=[row for row in json.loads(archive['literal_stdout']) if row['waist']==1 and row['exit_to_waist']==1.5 and row['cosines']==10]

def ode(s, y):
    b, w, theta=y
    return np.sin(theta),np.cos(theta),2*np.cos(theta)/b-30/19

def end(s,y):return y[2]
end.direction=-1; end.terminal=True
sol=solve_ivp(ode,[0,12.],[1.,0.,0.],events=end,rtol=2e-13,atol=2e-14,max_step=.015,dense_output=True)
s=np.linspace(0.,sol.t_events[0][0],12001)
b,w,angle=sol.sol(s)
co=np.asarray(ref['coefficients'])
g=co@np.cos(np.pi*np.arange(1,1+len(co))[:,None]*w[None,:]/w[-1])
g+=1.-g[0]

def frustum(eps,c):
    h=eps*g+c
    B=b+h*np.cos(angle);W=w-h*np.sin(angle)
    r0,r1=B[:-1],B[1:]
    A=(4*np.pi/3)*np.hypot(np.diff(B),np.diff(W))@(r0**2+r0*r1+r1**2)
    V=(np.pi/3)*np.diff(W)@(r0**3+r0**2*r1+r0*r1**2+r1**3)
    return A,V,B,W
A,V,_,_=frustum(0.,0.)
rows=[]
for eps in [0.,-.01,-.005,.005,.01]:
    c=0 if not eps else brentq(lambda c:frustum(eps,c)[1]-V,-abs(eps),abs(eps),xtol=2e-14)
    a,v,B,W=frustum(eps,c)
    test=ExchangeShape(eps)
    x=np.linspace(0., min(test.length,W[-1]),1001)
    predicted=CubicSpline(W,B)
    btest,dtest=test.radius(x)
    err=max(abs(predicted(x)-btest)); derr=max(abs(predicted(x,1)-dtest))
    rows.append(dict(step=eps,waist=B[0],flare=B[-1],length=W[-1],offset=c,area_change=a-A,
                     volume_roundoff=(v-V)/V,profile_error=err,derivative_error=derr))
    assert err<3e-7 and derr<1e-5 and abs(test.length-W[-1])<1e-9
print(json.dumps(rows,indent=2))
