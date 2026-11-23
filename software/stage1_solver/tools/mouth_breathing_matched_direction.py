"""Complete missing wave-empty force readings from Rayleigh slopes on a recorded central basis."""
import argparse
import hashlib
import json

import numpy as np

from mouth_action_breathing import volume_derivatives
from mouth_breathing_matched import displacement
from mouth_exit_witness import inputs
from mouth_live_waist_fields import Carrier
from mouth_poly_bonded_core import core_source
from mouth_poly_shallow_release import BIAS, CUT, ROOT, Order, Trial, emit, exterior, source
from mouth_sealed_shape_relax import Meridian


def run(radius):
    mesh=72
    folder=ROOT/f'_scratch/s11c/mouth_breathing_matched01/r{radius}n{mesh}'
    raw=(folder/'stdout').read_bytes()
    guard=json.loads((folder/'outcome.json').read_bytes())
    if guard['exitCode']!=1 or not guard['limitsVerified'] or b'Receiving order search failed' not in (folder/'stderr').read_bytes():
        raise ValueError('Expected an interrupted off-center receiver minimization')
    archive=list(map(json.loads,raw.splitlines()))
    record=archive[0]
    available={(r['relative_exit'],r['shape_step']):r for r in archive if r['kind']=='spectrum'}
    partial={(r['index'],r['exit_step']) for r in archive if r['kind']=='signed'}
    best=min(k[1] for k in available if k[0]==0.)
    center=available[0.,best]
    path=ROOT/f'_scratch/s11c/mouth_exit_release01/r{radius}n{mesh}/endpoint.npz'
    trial=Trial(radius,mesh,bonded_core=True)
    shifts,chi,_,digest,status,*_=inputs(path,radius,mesh,trial)
    if status!='balanced' or digest!=record['material_sha256']:
        raise ValueError('Different mouth preparation')
    polar_parent=ROOT/'_scratch/s11c/mouth_action_breathing01/run/stdout'
    parent_raw=polar_parent.read_bytes()
    if hashlib.sha256(parent_raw).hexdigest()!=record['mass_sha256']:
        raise ValueError('Unexpected registered inertia')
    inertia={r['relative_exit']:np.array(r['mass']) for r in map(json.loads,parent_raw.splitlines())
             if r['kind']=='basis_state' and r['n']==mesh and r['radius']==radius}
    shape=Meridian(radius,trial.x+shifts[:4])
    b,init=shape.exit,record['original_exit_scale']
    t=np.linspace(0,1,2*mesh+1)
    x=np.linspace(0,1,mesh+1)
    co=trial.co.copy(); co[0]+=shifts[4]
    cr=trial.cr.copy(); cr[0]+=shifts[5]
    emit(kind='input',radius=radius,mesh=mesh,input_sha256=digest,partial_stdout_sha256=hashlib.sha256(raw).hexdigest(),
         central_mesh_step=best,retained=sorted(partial),original_exit_scale=init,material_sha256=digest)
    emit(**center)

    def energy(offset,q):
        config=shifts[:4].copy()
        config[0]+=np.log1p(q[0])
        config[1]+=np.log1p(offset+q[1]*init/b)
        config[2]+=q[2]
        moved=Meridian(radius,trial.x+config)
        inside=source(moved,co,t,x,False)
        inside+=core_source(moved,co,*trial.core_reference,t,x)
        outside=exterior(moved,mesh,CUT,side_order_work=lambda w,y:source(moved,cr,w/moved.length,y,True))
        ef,_,info=Carrier(moved,mesh,2*mesh,trial.mesh).response(chi,BIAS)
        if abs(info['terminal_current']-info['current'])>1e-6*abs(info['current']):
            raise ArithmeticError('Native current mismatch')
        return (Order(moved,mesh).energy_gradient(chi)[0]+np.sum(chi*inside)
                +outside['energy']+outside['static_work']+ef)

    for step,h in ((.012,.002),(.006,.001)):
        offsets=(-step,step)
        ev=[energy(u,np.zeros(3)) for u in offsets]
        load=(ev[1]-ev[0])/(2*step*b)
        for idx,vec in enumerate(center['mass_normalized_vectors']):
            if (idx,step) in partial:
                continue
            u=np.array(vec)
            v=u/max(abs(u))
            freq=[]
            for j,offset in enumerate(offsets):
                E0=ev[j]
                em,ep=energy(offset,-h*v),energy(offset,h*v)
                mass=v@inertia[offset]@v
                curv=((ep-E0)+(em-E0))/h**2
                if curv<=0:raise ValueError('Projected curvature is nonpositive')
                freq.append(np.sqrt(curv/mass))
                emit(kind='work',index=idx,relative_exit=offset,step=h,
                     minus=em-E0,plus=ep-E0,curvature=curv,mass=mass)
            omega=center['frequency'][idx]
            force=-(freq[1]-freq[0])/(2*step*b)
            report=dict(kind='signed',index=idx,exit_step=step,omega=omega,outward_force_per_action=force,
                        current_closing_load=load,volume_change_at_unit_modal_coordinate=volume_derivatives(shape,init)@u,
                        frozen_central_vector=True,rayleigh_frequencies=freq)
            if force>0:
                J=load/force
                report.update(conditional_action=J,displacement=displacement(shape,cr,init,u,omega,J))
            emit(**report)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('radius',type=int,choices=(1,2))
    run(parser.parse_args().radius)
