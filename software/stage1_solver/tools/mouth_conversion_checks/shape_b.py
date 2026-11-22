from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp,simpson
pi4=4*np.pi
def shape(R):
    p=30/(19*R)
    def fun(s,y):
        b,w,a=y
        return [np.sin(a),np.cos(a),2*np.cos(a)/b-p]
    def event(s,y): return y[2]
    event.terminal=True;event.direction=-1
    res=solve_ivp(fun,(0,12*R),(R,0.,0.),events=event,rtol=1.5e-12,atol=1e-13,dense_output=True)
    s=np.linspace(0,res.t_events[0][0],4501)
    b,w,a=res.sol(s); c=np.cos(a); da=2*c/b-p
    L=w[-1]; alpha=L**2*pi4*(R**2-p*R**3/3.)
    def work(d):
        bt=b+d*c
        A=pi4*simpson(bt**2*(1-d*da),x=s)
        V=pi4/3*simpson(bt**3*c*(1-d*da),x=s)
        return A-p*V+alpha/L
    h=.0007*R; second=(work(h)+work(-h)-2*work(0))/h**2
    theory=-pi4*simpson(b*b*(da**2+2*(c/b)**2),x=s)
    print('arc evolution',R,L,b[-1],alpha,second,theory)
    assert abs(b[-1]/R-1.5)<1e-9 and abs(second/theory-1)<2e-6
for r in (1,2):shape(r)
base=Path('_scratch/s11c/downhill_front01')
def rows(name,code):
    d=base/name
    for f in ['limit-validation','duration-validation']:
        assert json.loads((d/(f+'.json')).read_text())['verified']
    assert json.loads((d/'limit-validation.json').read_text())['actual']['memory.swap.max']=='0'
    assert json.loads((d/'duration-validation.json').read_text())['wallDeadlineSeconds'] is None
    g=json.loads((d/'outcome.json').read_text())
    assert g['exitCode']==g['childOutcome']['exitCode']==code and g['limitsVerified']
    return [json.loads(x) for x in (d/'stdout').read_text().splitlines()]
for row in rows('planar',0):
    h=row['epsilon']
    # Infinite planar limit with n=1, delta=1, mobility 1/2. Integrated
    # squared logistic gradient 1/6, so energy loss per area is 2j²/6.
    assert abs(row['ordered_current']-3*h)<2e-10 and abs(row['conversion_loss']-3*h*h)<1e-10
for fine,raw,inner in [('small_refined','small_simpson',100),('large_refined','large',200)]:
    ref=rows(fine,0)[0]; drive=rows(raw,0)
    assert len(drive)==5 and abs(ref['half_order']-drive[0]['half_order'])<5e-6
    for row in drive:
        q=row['front_current']*(row['half_order']/inner)**3
        b=row['epsilon']+.5*q*q*(1-(inner/400)**6)
        assert abs(q*b-row['conversion_loss'])<2e-9
        if row['external_holder']:
            assert row['velocity_on_release']*(row['half_order']-ref['half_order'])>0
    for full,half in zip(drive[1:3],drive[3:]):
        assert abs(full['velocity_on_release']/half['velocity_on_release']-2)<1e-5
    print('driven B',inner,'front',ref['half_order'],'head',b,'velocities', [v['velocity_on_release'] for v in drive[1:]])
assert rows('small',1)[0]['max_total_normal_balance']>2e-7
