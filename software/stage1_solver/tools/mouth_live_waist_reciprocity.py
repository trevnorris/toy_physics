"""Independently re-solve energies after a signed order perturbation on the saved waist."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_live_waist_fields import CoupledLoads, BIAS
from mouth_sealed_shape_relax import Meridian

ROOT=Path(__file__).resolve().parents[3]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--radius',type=int,choices=(1,2),required=True)
    p.add_argument('--stage',choices=('late','tight'),required=True)
    p.add_argument('--mesh',type=int,required=True)
    a=p.parse_args()
    record=ROOT/'research/pde_ledger_v3/_measurements/mouth_sealed_shape_relax.json'
    domain=json.loads(record.read_bytes())['training'][str(a.radius)]
    state=domain['states'][-2] if a.stage=='late' else domain['extension_states'][-1]
    shape=Meridian(a.radius,state['coordinates'])
    nr=4*a.mesh; nz=2*nr
    x,t=np.meshgrid(np.linspace(0,1,nr+1),np.linspace(0,1,nz+1))
    order=.5+.49*np.tanh((1-x)*shape.profile(t)[0]/.08)*np.tanh((1-t)*shape.length/.08)
    eta=np.sin(2*np.pi*t)*np.cos(np.pi*x/2)**2
    loads=CoupledLoads(shape,nr,nz,a.mesh,4.2,target_frequency=state['wave_frequency'])
    total,g,reading=loads.evaluate(order)
    base_seed=loads.wave.seed
    cc,cg,extra=loads.evaluate(order,action=0)
    ce=float((cg*eta).sum()); we=float(((g-cg)*eta).sum())
    shape_info=dict(radius=a.radius,stage=a.stage,field_mesh=a.mesh,order_grid=[nz+1,nr+1])
    print(json.dumps(dict(**shape_info,full=reading,carrier_work=ce,wave_work=we)),flush=True)
    tests=[]
    for h in (.002,.001):
        sides=[]
        for sign in (-1,1):
            loads.wave.seed=base_seed
            _,_,info=loads.evaluate(order+sign*h*eta)
            sides.append(info)
        results=dict(half_step=h,carrier_from_energy=(sides[1]['carrier_energy']-sides[0]['carrier_energy'])/(2*h),
                     wave_from_energy=(sides[1]['wave_frequency']-sides[0]['wave_frequency'])/(2*h))
        print(json.dumps(results),flush=True)
        tests.append(results)
    reverse=loads.evaluate(order,action=0,bias=-BIAS)
    zero=loads.evaluate(order,action=0,bias=0)
    print(json.dumps(dict(reverse_current=reverse[2]['current'],reverse_energy=reverse[0],
                          zero_energy=zero[0],zero_gradient=max(abs(zero[1].ravel())))),flush=True)
    assert abs(reverse[0]-cc)<1.e-8 and abs(zero[0])<1.e-11 and abs(reverse[2]['current']+extra['current'])<1.e-7
    for name,work in [('carrier',ce),('wave',we)]:
        numerical=(4*tests[1][name+'_from_energy']-tests[0][name+'_from_energy'])/3
        if abs(numerical-work)>2e-6*max(.01,abs(work)):
            raise ArithmeticError(f'{name} order-work test failed: {numerical} vs {work}')


if __name__=='__main__':
    main()
