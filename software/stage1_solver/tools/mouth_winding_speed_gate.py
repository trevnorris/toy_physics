"""Kinematic limits for a conserved phase winding through the two saved terminal three-balls."""
import hashlib
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad

from mouth_curved_steady import SteadyMouth
from mouth_exit_witness import inputs
from mouth_poly_shallow_release import ROOT, Trial, emit
from mouth_sealed_shape_relax import Meridian


def preparation(radius):
    ledger=ROOT/'research/pde_ledger_v3/_measurements'
    record=json.loads((ledger/'mouth_breathing_exit_support.json').read_bytes())
    trial=Trial(radius,72,bonded_core=True)
    path=ROOT/f'_scratch/s11c/mouth_exit_release01/r{radius}n72/endpoint.npz'
    state,order,_,digest,status,*_=inputs(path,radius,72,trial)
    if status!='balanced' or digest!=record['physical_preparations'][f'r{radius}n72']:
        raise ArithmeticError('Different held material state')
    load_case='r1n72_direction' if radius==1 else 'r2n72'
    run=record['runs']['mouth_breathing_matched01/'+load_case]
    raw=run['literal_stdout'].encode()
    if run['guard_exit'] or run['child_exit'] or hashlib.sha256(raw).hexdigest()!=run['stdout_sha256']:
        raise ArithmeticError('Unverified matched closing load')
    load,=[r['current_closing_load'] for r in map(json.loads,raw.splitlines())
           if r['kind']=='signed' and r['index']==0 and r['exit_step']==.006]
    return Meridian(radius,trial.x+state[:4]),order,digest,load


def capacity(b,m):
    """Max subsonic normal mass current with n<=1, including both polar hemispheres."""
    a=abs(m)
    volume=4*np.pi*b**3/3
    if a>=b:
        return 0.,0.,0.
    annular_fraction=(1-(a/b)**2)**1.5
    # Cylindrical radius s. The remaining z interval is 2 sqrt(b²−s²).
    value,err=quad(lambda s:4*np.pi*np.sqrt((b*b-s*s)*(s*s-a*a)),
                   a,b,epsabs=2e-13,epsrel=2e-12,limit=200)
    # Independent integral at constant z: primitive of s sqrt(1−a²/s²).
    def section(z):
        s=np.sqrt(b*b-z*z)
        if a==0:
            return np.pi*s*s
        return np.pi*(s*np.sqrt(max(0.,s*s-a*a))-a*a*np.arccosh(max(1.,s/a)))
    check,_=quad(section,-np.sqrt(b*b-a*a),np.sqrt(b*b-a*a),
                 epsabs=2e-13,epsrel=2e-12)
    if abs(value-check)>2e-10*volume or err>2e-10*volume:
        raise ArithmeticError('Independent three-disk capacity integrals disagree')
    return annular_fraction,value,check


def unwound_one(shape,order,n):
    state=SteadyMouth(shape,order,n)
    edges=state.radial[0]*shape.exit
    roots,weights=leggauss(4)
    r=edges[:-1,None]+np.diff(edges)[:,None]*(1+roots)/2
    area=4*np.pi*r*r*np.diff(edges)[:,None]*weights/2
    j=state.disk(r)
    density=state.details('core',r,np.full_like(r,shape.length))[0]
    normal=j/density**3
    flux=(area*j).sum()
    emit(kind='unwound',radius=1,mesh=n,exit=shape.exit,throughput=flux,
         max_normal_mach=normal.max(),area_fraction_normal_above_sound=area[normal>1].sum()/area.sum(),
         flux_fraction_normal_above_sound=(area*j)[normal>1].sum()/flux,
         throughput_disagreement=abs(flux-state.info['terminal_current']),n_range=[density.min(),density.max()])
    if np.min(j)<=0 or abs(flux-state.info['terminal_current'])>1e-6*flux:
        raise ArithmeticError('Invalid flux at the held zero-winding control')
    return flux


def run():
    for radius in (1,2):
        shape,order,digest,load=preparation(radius)
        b=shape.exit
        emit(kind='input',radius=radius,exit=b,length=shape.length,
             material_sha256=digest,inward_load=load)
        if radius==1:
            for n in (20,32,48):
                throughput=unwound_one(shape,order,n)
        else:
            record=json.loads((ROOT/'research/pde_ledger_v3/_measurements/mouth_curved_exit_gate.json').read_bytes())
            prior=record['runs']['mouth_curved_exit01/steady']
            raw=prior['literal_stdout'].encode()
            if prior['guard_exit'] or prior['child_exit'] or hashlib.sha256(raw).hexdigest()!=prior['stdout_sha256']:
                raise ArithmeticError('Unverified radius-two speed reference')
            baseline,=[row for row in map(json.loads,raw.splitlines()) if row['offset']==0 and row['field_mesh']==48]
            if baseline['material_sha256']!=digest or abs(b-baseline['disk_radius'])>1e-12:
                raise ArithmeticError('Mismatched disk geometry')
            throughput=baseline['flow']
            emit(kind='unwound_adopted',radius=2,mesh=48,throughput=throughput,
                 max_normal_mach=baseline['max_normal_mach'],
                 area_fraction_normal_above_sound=baseline['area_fraction_local_mach_over_one'],
                 flux_fraction_normal_above_sound=baseline['throughflow_fraction_local_mach_over_one'],
                 prior_sha256=prior['stdout_sha256'])
        for m in (0,1,-1,2):
            area,cap,second=capacity(b,m)
            emit(kind='kinematic_bound',radius=radius,winding=m,phase_circulation=2*np.pi*m,
                 rim_equator_speed_floor=abs(m)/b,
                 volume_fraction_total_mach_at_least_one_lower_bound=1-area,
                 max_subsonic_mass_throughput=cap,independent_capacity=second,
                 capacity_over_unwound_throughput=cap/throughput)
            if m==0 and abs(cap-4*np.pi*b**3/3)>2e-11:
                raise ArithmeticError('Zero-winding capacity control')


if __name__=='__main__':
    run()
