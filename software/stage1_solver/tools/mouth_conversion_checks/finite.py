import numpy as np, json, hashlib
from pathlib import Path
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

root=Path('_scratch/s11c/conversion_routes_shape01')
def raw(d):
    d=root/d
    o=json.loads((d/'outcome.json').read_text())
    assert o['exitCode']==o['childOutcome']['exitCode']==0 and o['limitsVerified']
    return json.loads((d/'stdout').read_text())

for case in raw('exchange'):
    R=case['waist'];g=case['exit_to_waist']
    if case['cosines']!=6 and R==1:continue
    p=3*(g+1)/(R*(g*g+g+1))
    def rhs(_,y):
        b,w,a=y
        return [np.sin(a),np.cos(a),2*np.cos(a)/b-p]
    def event(s,y):return y[2]
    event.terminal=True;event.direction=-1
    sol=solve_ivp(rhs,[0,16*R],[R,0.,0.],events=event,rtol=1e-12,atol=1e-13,dense_output=True)
    s=np.linspace(0,sol.t_events[0][0],20001); b,w,a=sol.sol(s); L=w[-1]
    vec=np.array(case['coefficients']);k=(1+np.arange(len(vec)))[:,None]*np.pi/L
    z=vec@np.cos(k*w)
    # Normalize the profile at the waist and match the named cosine offsets.
    z+=1-z[0]
    def surface(eps,offset):
        dr=eps*z+offset;B=b+dr*np.cos(a); W=w-dr*np.sin(a);b0,b1=B[:-1],B[1:]
        seg=np.hypot(np.diff(B),np.diff(W))
        A=(4*np.pi/3)*sum(seg*(b0*b0+b0*b1+b1*b1))
        V=(np.pi/3)*sum(np.diff(W)*(b0**3+b0*b0*b1+b0*b1*b1+b1**3))
        return A,V
    A,V=surface(0,0)
    assert abs(A/case['area']-1)<1e-8
    for ref in case['controls']:
        e=ref['fraction']*R;change=[]
        for sign in [-1,1]:
            c=brentq(lambda t:surface(sign*e,t)[1]-V,-e,e,xtol=1e-14)
            at,vt=surface(sign*e,c)
            assert abs(vt/V-1)<2e-12
            change.append(at-A)
        print('finite frusta',R,g,e,change)
        np.testing.assert_allclose(change,ref['area_changes'],rtol=2e-6,atol=2e-11)
