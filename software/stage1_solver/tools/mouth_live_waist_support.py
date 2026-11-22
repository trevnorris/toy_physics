"""Bounded-order continuation with the transverse pencil restricted to its supported nodes."""
from collections import Counter
import argparse
from pathlib import Path

import numpy as np
from scipy.ndimage import label

from mouth_live_waist_active_wave import ActiveWave
from mouth_live_waist_order import BIAS, Order, Meridian, emit, load_for, state
from mouth_live_waist_bounded import bounded


def active_load(*args):
    loads=ActiveWave.install(load_for(*args))
    solve=loads.evaluate
    errors=Counter()
    def trial(chi,action,bias):
        old=loads.wave.seed,loads.wave.target,loads.carrier.f.copy()
        try:
            return solve(chi,action,bias)
        except (ValueError,ArithmeticError,RuntimeError) as exc:
            loads.wave.seed,loads.wave.target,loads.carrier.f=old
            errors[type(exc).__name__+': '+str(exc).split(':')[0]]+=1
            # Reject this trial. The physical energy is undefined without its occupied mode.
            return 1.e8,np.zeros_like(chi),{'undefined_energy':True}
    loads.evaluate=trial
    loads.rejected=errors
    return loads


def release(order,loads,chi,energy):
    groups,count=label((chi<=1e-12)&order.free)
    changes=[]
    for i in range(1,count+1):
        region=groups==i
        for h in (.0001,.0002):
            c=chi+h*region
            e0,_=order.energy_gradient(c)
            ef,_,info=loads.evaluate(c,1.,BIAS)
            if 'undefined_energy' in info:
                changes.append(dict(group=i,step=h,undefined=True))
            else:
                changes.append(dict(group=i,step=h,energy_slope=(ef+e0-energy)/h))
    return changes


def run(a):
    s=state(1,a.point)
    params=np.array(s['coordinates'])
    shape=Meridian(1,params); order=Order(shape,a.mesh)
    loads=active_load(shape,a.mesh,a.field_mesh,s['wave_frequency'])
    chi,flow=order.condition(order.seed,loads,0.,BIAS)
    emit(kind='flow',**flow)
    for action in (.125,.25,.375,.5,.625,.75,.875,1.):
        chi,metrics=bounded(order,chi,loads,action)
        if 'undefined_energy' in metrics or metrics['total_inside_and_drives']>1e7:
            raise ArithmeticError('No accepted wave/flow state at the candidate')
        emit(kind='coupled',action=action,rejections=dict(loads.rejected),**metrics)
    np.savez_compressed(a.save,order=chi,frequency=metrics['wave_frequency'],
                        vector=loads.wave.seed[1],pole=loads.wave.seed[0])
    if metrics.get('fraction_at_zero',0):
        checks=release(order,loads,chi,metrics['total_inside_and_drives'])
        emit(kind='opening_check',changes=checks)
        if any('undefined' in q or q['energy_slope']<-.001 for q in checks):
            raise ArithmeticError('One-sided reordering test did not resolve a balanced state')
    for coord in ('waist','depth'):
        for h in (-.002,-.001,.001,.002):
            p=params.copy(); p[2 if coord=='waist' else 0]+=h
            q=Meridian(1,p); order2=Order(q,a.mesh)
            field=active_load(q,a.mesh,a.field_mesh,s['wave_frequency'],loads.wave.seed)
            c,shifted=bounded(order2,chi,field,1.)
            if 'undefined_energy' in shifted or shifted['total_inside_and_drives']>1e7:
                raise ArithmeticError('Signed state is undefined')
            emit(kind='signed',coordinate=coord,offset=h,**shifted)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mesh',type=int,required=True)
    ap.add_argument('--field-mesh',type=int,required=True)
    ap.add_argument('--point',choices=('mid','tight'),required=True)
    ap.add_argument('--save',type=Path,required=True)
    run(ap.parse_args())
