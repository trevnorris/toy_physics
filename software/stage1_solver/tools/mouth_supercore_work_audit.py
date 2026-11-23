"""Independent full-field signed variations of the radius-four wound material and fixed collar."""
import argparse
import contextlib
import hashlib
import io
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import map_coordinates

from mouth_exchange_carrier import BIAS
from mouth_poly_bonded_core import core_source
from mouth_poly_shallow_release import CUT, Order, emit, exterior, source
from mouth_sealed_shape_relax import Meridian
from mouth_wound_carrier import Wound


def state(path):
    guard=json.loads((path.parent/'outcome.json').read_bytes())
    limits=json.loads((path.parent/'limit-validation.json').read_bytes())
    if guard['exitCode'] or guard['childOutcome']['exitCode'] or not limits['verified'] or int(limits['actual']['memory.swap.max']):
        raise ValueError('Unverified material input')
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    if digest not in (path.parent/'stdout').read_text():
        raise ValueError('Unidentified material input')
    with np.load(path,allow_pickle=False) as z:
        if int(z['radius'])!=4:
            raise ValueError('Wrong entrance')
        return z['shape_params'],z['order'],z['registration'],digest,float(z['shape_extra']) if 'shape_extra' in z else 0.


def energy(shape,chi,registration,ref,ref_registration,mesh,polar):
    # Each of the two perturbed sides gets its own converged carrier and
    # receiving material. Keep ordered labels fixed on the mapped shape.
    n=chi.shape[1]-1
    t,x=np.linspace(0.,1.,2*n+1),np.linspace(0.,1.,n+1)
    elastic=source(shape,registration,t,x,False)
    bonded=core_source(shape,registration,ref,ref_registration,t,x)
    receiver=exterior(shape,n,CUT,side_order_work=lambda w,s:
                      source(shape,registration,w/shape.length,s,True))
    bulk=Wound(shape,chi,1,mesh,polar,28.,22.)
    with contextlib.redirect_stdout(io.StringIO()):
        _,_,carrier=bulk.solve(BIAS)
    terms=dict(order=Order(shape,n).energy_gradient(chi)[0],inner=np.sum(chi*elastic),
               bonded=np.sum(chi*bonded),outer=receiver['energy']+receiver['static_work'],
               carrier=carrier['excess_energy'])
    if receiver['residual']>1e-5 or carrier['residual']>1e-7:
        raise ArithmeticError('Unresolved carrier or receiving material')
    return terms,carrier['total_current']


def contract_field(chi):
    # Move outward from the center of the prepared pocket to sample order.
    # Raising order in the part that used to be its rim shrinks the pocket.
    nt,nr=np.subtract(chi.shape,1)
    iz,ix=np.indices(chi.shape,dtype=float)
    t,x=iz/nt,ix/nr
    c=map_coordinates(chi,[nt*(.5+1.03*(t-.5)),nr*1.03*x],order=1,mode='nearest')
    shape=np.maximum(0.,1-((t-.5)/.39)**2)**2*np.maximum(0.,1-(x/.65)**2)**2
    v=np.maximum(c-chi,0)*shape
    v[[0,-1]]=0.
    v[:,-1]=0.
    return v


def run(current,initial,mesh,polar,tests,steps):
    params,chi,co,sha,reference_extra=state(current)
    q0,_,refco,refsha,initial_extra=state(initial)
    base,ref=Meridian(4.,params,reference_extra),Meridian(4.,q0,initial_extra)
    emit(kind='input',source=str(current),reference=str(initial),material_sha256=sha,reference_sha256=refsha,
         order_mesh=chi.shape[1]-1,carrier_mesh=mesh,polar=polar,relative_steps=steps)
    if 'heal' in tests:
        v=contract_field(chi)
        emit(kind='contraction_support',positive_nodes=int((v>0).sum()),max_shift=float(np.max(v)))
        for h in [0.,.25,.5]:
            e,D=energy(base,chi+h*v,co,ref,refco,mesh,polar)
            if not h:zero=e
            delta={k:v-zero[k] for k,v in e.items()}
            emit(kind='contracting_old_pocket',fraction=h,terms=delta,energy_change=sum(delta.values()),current=D)
    for test in (test for test in tests if test!='heal'):
        for h in steps:
            states=[]
            for sign in (-1,1):
                moved=params.copy()
                c=co.copy()
                extra=reference_extra
                if test in ('axial','exit'):
                    i=int(test=='exit')
                    moved[i]+=np.log1p(sign*h)
                    size=h*(base.exit if i else base.length)
                elif test=='collar':
                    c[0]+=sign*h
                    size=h
                elif test=='quartic':
                    extra+=sign*h/4
                    size=h
                shape=Meridian(4.,moved,extra)
                e,D=energy(shape,chi,c,ref,refco,mesh,polar)
                states.append(e)
                emit(kind='perturbation',test=test,h=h,sign=sign,terms=e,current=D)
            outward={k:(v-states[1][k])/(2*size) for k,v in states[0].items()}
            emit(kind='signed_force',test=test,h=h,force=sum(outward.values()),components=outward)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('current',type=Path)
    p.add_argument('reference',type=Path)
    p.add_argument('--mesh',type=int,default=9)
    p.add_argument('--polar',type=int,default=7)
    p.add_argument('--tests',nargs='+',choices=['heal','axial','exit','quartic','collar'],required=True)
    p.add_argument('--steps',nargs='+',type=float,default=[.001])
    a=p.parse_args()
    run(a.current,a.reference,a.mesh,a.polar,a.tests,a.steps)
