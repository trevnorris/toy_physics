"""Distinguish angular energy from re-settling and verify mobile-order bounds on R4."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from mouth_bonded_front_adjoint import set_front
from mouth_fixed_wall_angular import AngularBonded, equilibrate_scalar
from mouth_joined_front_audit import Front
from mouth_joint_polar_shape import Trial
from mouth_polar_exit_readout import accepted
from mouth_supercore_joined_front import check_guard, emit, file_sha


def run(state,seed,end,short_steps=False):
    check_guard(end)
    prep,_,x,d,f = accepted(state,seed)
    with np.load(end) as saved:
        if saved['input_sha'].item()!=file_sha(seed) or not np.array_equal(saved['x'],x):
            raise ValueError('Mismatched fixed-wall end state and seed')
        dn,u = saved['displacement'],saved['field']
        chosen = float(saved['internal_amplitude'])
    base=Front(state,tuple(prep['field_grid']),1)
    n,k=(d.shape[1]-1)//3,d.shape[0]-1
    t=Trial(base,state,x,n,k,False)
    m=t.model
    b=AngularBonded(state,n,k,shape=t.shape,epsilon=x[5])
    q=d[b.mask]*b.scale
    qn=dn[b.mask]*b.scale
    if not short_steps:
        audit_readouts(t,m,b,base,f,u,d,dn,q,qn,chosen,n,k,seed,end)
    verify_direction(b,qn,(4e-6,1e-6,2.5e-7) if short_steps else (2e-5,1e-5),
                     short_steps)


def audit_readouts(t,m,b,base,f,u,d,dn,q,qn,chosen,n,k,seed,end):
    set_front(b,m,f)
    b.energy(q)
    offset = b.angular_energy_delta
    s=equilibrate_scalar(b,q)
    b.energy(q)
    scalar_only = b.angular_energy_delta-offset
    shape_energy_before=t.energy(d,f)
    baseline_after=t.energy(dn,u)
    b.front=None
    b.internal_amplitude=chosen
    b.energy(qn)
    new_s=b.angular_energy_delta
    b.internal_amplitude=0.
    b.energy(qn)
    angular_final=new_s-b.angular_energy_delta
    complete = baseline_after-shape_energy_before+angular_final
    emit(kind='energy_partition',input_sha=file_sha(seed),end_sha=file_sha(end),grid=[3*n,k],
         scalar_minimum_at_original_d_and_f=s,fixed_d_f_angular_drop=scalar_only,
         final_angular_drop_at_final_d=angular_final,
         baseline_material_field_energy_shift=baseline_after-shape_energy_before,
         all_terms_change=complete,correction_to_frozen_meridian=complete-scalar_only,
         peak_internal_radial_motion=float(abs(chosen)*27/256),
         peak_pole_minus_equator=float(abs(chosen)*81/512),
         peak_field_change=float(abs(f-u).max()),peak_material_nodal_change=float(abs(dn-d).max()))
    interp=RegularGridInterpolator((np.linspace(0,1,len(base.c)),np.linspace(0,1,base.c.shape[1])),base.c)
    mobile = []
    for block in (1,2):
        xn,_,zn=np.meshgrid(m.rr[0],m.a,m.z[block],indexing='ij')
        original=(interp((zn/t.shape.length,xn)) if block==1
                  else m.terminal_profile(xn,zn-t.shape.length))
        ids=m.ids[block,0]
        mask=m.masks[2,ids]
        mobile.extend((original+u[2,ids])[mask].tolist())
    emit(kind='physical_order',mobile_terminal_node_range=[min(mobile),max(mobile)],
         mobile_terminal_incident_and_receiving_trace_samples=len(mobile),
         shaft_reference_node_range=[float(base.c.min()),float(base.c.max())],
         deep_unperturbed_axis_order=float(base.c[len(base.c)//2,0]),
         receiver_native_range=[t.receiver['minimum_order'],t.receiver['maximum_order']])
    if min(mobile)<0 or max(mobile)>1:raise ArithmeticError('Unphysical mobile order node')


def verify_direction(b,qn,steps,short_steps):
    # Isolate the angular correction from the spherical material to avoid
    # differencing large extensive energies at this weak displacement.
    tt,xx=np.meshgrid(b.t,b.x,indexing='ij')
    smooth=xx*(1-xx)*np.sin(np.pi*tt)
    v=np.stack((smooth,smooth*np.cos(np.pi*xx)),axis=-1)[b.mask]*b.scale
    b.internal_amplitude=.002
    _,grad=b.energy(qn)
    native=b.native.copy()
    from mouth_bonded_material_motion import Bonded
    Bonded.energy(b,qn)
    projected=(native-b.native)[b.mask] @ (v/b.scale)
    for h in steps:
        vals=[]
        for sign in (-1,1):
            b.energy(qn+sign*h*v)
            vals.append(b.angular_energy_delta)
        signed=(vals[1]-vals[0])/(2*h)
        emit(kind='correction_meridional_work',s=b.internal_amplitude,h=h,
             signed=signed,analytic=projected,difference=signed-projected)
        if abs(signed-projected)>(2e-8 if short_steps and h==steps[-1] else 1.e-6):
            raise ArithmeticError('Independent angular meridional work does not close')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path);p.add_argument('pairs',type=Path,nargs='+')
    p.add_argument('--short-steps',action='store_true')
    a=p.parse_args()
    for start,end in zip(a.pairs[::2],a.pairs[1::2],strict=True):
        run(a.state,start,end,a.short_steps)
