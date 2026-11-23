"""Track the lowest radius-two breathing resonance with joint side and disk radiation.

The uniform-current front is a comparison law. It deliberately does not assign a
subsonic planar response to the fast local rim. Outgoing stiffness is continued
to the free resonance; real-frequency energy flux is measured separately.
"""
import argparse
from functools import lru_cache

import numpy as np
from scipy.linalg import eig, eigh

from mouth_breathing_matched import displacement
from mouth_curved_cap_ports import port
from mouth_curved_sound import CurvedSound
from mouth_curved_steady import SteadyMouth, reference
from mouth_poly_shallow_release import emit


def as_pairs(x):
    x=np.asarray(x)
    return np.stack([x.real,x.imag],axis=-1).tolist()


def unit_mode(stiffness,mass,reference_vector):
    values, vectors=eig(stiffness,mass)
    norms=np.sqrt(np.einsum('ik,ij,jk->k',vectors.conj(),mass,vectors).real)
    vectors=vectors/norms
    overlap=reference_vector.conj()@mass@vectors
    index=np.argmax(abs(overlap))
    vector=vectors[:,index]*np.exp(-1j*np.angle(overlap[index]))
    return values[index],vector,float(abs(overlap[index]))


def follow_pole(H,M,response):
    vals, vec=eigh(H,M)
    z=np.sqrt(vals[0])+0j
    ref=vec[:,0]
    trace=[]
    for fraction in (.125,.25,.5,.75,1.):
        def mismatch(freq):
            val,v,ov=unit_mode(H+fraction*response(freq),M,ref)
            return val-freq**2,v,ov
        old=z*(1+2e-4)
        prev=mismatch(old)[0]
        for iteration in range(18):
            value,v,ov=mismatch(z)
            if abs(value)<2e-8*(1+abs(z)**2):
                break
            if value==prev or abs(z-old)==0:
                raise ArithmeticError('No slope in coupled pole solve')
            step=value*(z-old)/(value-prev)
            if abs(step)>.4*abs(z):
                step*=.4*abs(z)/abs(step)
            old,prev,z=z,value,z-step
        else:
            raise ArithmeticError('Unresolved continuation of breathing pole')
        if ov<.8 or z.real<=0:
            raise ArithmeticError('Lost breathing mode along coupling continuation')
        trace.append([fraction, *as_pairs(z),ov,float(abs(value))])
        ref=v
    return z,ref,trace


def run(args):
    shape,chi,native,scale,ext,rows,digest=reference(args.offset)
    state=SteadyMouth(shape,chi,native)
    step=.002 if abs(args.offset)>.01 else .001
    row,=[r for r in rows if r.get('relative_exit')==args.offset and r['shape_step']==step]
    H,M=np.array(row['hessian']),np.array(row['mass'])
    mean=state.info['terminal_current']/(4*np.pi*shape.exit**3/3)

    def sheet_at_mean(w,r):
        p=port(w,mean)
        return np.broadcast_to(p,(len(r),5,5))
    args_dict=vars(args)
    emit(kind='input',args=args_dict,material_sha256=digest,mean_current=mean,
         exterior_side='sealed',terminal_front=args.front,
         unforced_frequency=row['frequency'][0],physical_exit=shape.exit, mass=M.tolist())
    op=CurvedSound(shape,scale,n=args.mesh,box=args.box,pml=args.pml,stretch=args.stretch,
                   base=state,sheet=sheet_at_mean if args.front=='mean' else None)

    @lru_cache(maxsize=32)
    def response(w):
        return op.stiffness(w)
    z,v,trace=follow_pole(H,M,response)
    w=z.real
    report=op.evaluate(w)
    v*=np.exp(-1j*np.angle(v[abs(v).argmax()]))
    def projection(A):
        return float(np.vdot(v,A@v).real)
    def derivative(h):
        return (response(w+h)-response(w-h))/(2*h)-2*w*M
    energy=-w*projection(derivative(.001*w))/4
    energy_coarse=-w*projection(derivative(.002*w))/4
    if energy<=0 or z.imag>=0:
        raise ArithmeticError('No decaying positive-energy branch')
    fraction=-np.expm1(4*np.pi*z.imag/w)
    powers={name:projection(report[name]) for name in ('feed_power','bulk_power','feed_cut','bulk_cut')}
    work=-.5*w*np.vdot(v,report['W']@v).imag
    quantities={key:as_pairs(report[key]@v) for key in ('relative_current','absolute_current','front_force',
                                                         'pressure_minus','pressure_plus')}
    # Actual four-dimensional volume, signed relative to the outward normal of
    # the cavity. Axial motion also contributes a side term on the sloping wall.
    t,weight=np.polynomial.legendre.leggauss(32)
    t=(t+1)/2
    b,slope,_=shape.profile(t)
    radial=scale*t*v[1]+4*shape.neck*t*(1-t)*v[2]
    side=-1j*w*shape.length*np.sum(weight/2*4*np.pi*b*b*(radial-slope*shape.length*t*v[0]))
    geometric_disk=-1j*w*shape.length*v[0]*report['disk_weights'].sum()
    anchor=2.96  # Common real driving frequency across signed geometries.
    anchor_stiffness=H+response(anchor)-anchor**2*M
    anchor_slope=(response(anchor+.002)-response(anchor-.002))/.004-2*anchor*M
    emit(kind='mode',anchor_omega=anchor,anchor_stiffness=as_pairs(anchor_stiffness),
         anchor_slope=as_pairs(anchor_slope),offset=args.offset,pole=as_pairs(z),tracking=trace,mass_unit_vector=as_pairs(v),
         outgoing_Q=-w/(2*z.imag),pole_cycle_loss=fraction,peak_unit_stored_energy=energy,
         stored_energy_double_step=energy_coarse,real_frequency_external_work=float(work),
         real_frequency_power=powers,
         power_per_period_over_stored_energy=float(2*np.pi*(powers['feed_power']+powers['bulk_power'])/(w*energy)),
         relative_cut_discrepancy=float(abs(powers['feed_cut']+powers['bulk_cut']-work)/abs(work)),
         source_stress_work_discrepancy=float(abs(powers['feed_power']+powers['bulk_power']-work)/abs(work)),
         side_sweep=as_pairs(side),disk_sweep=as_pairs(geometric_disk),
         disk_relative_current=as_pairs(report['disk_weights']@(report['relative_current']@v)),
         disk_absolute_current=as_pairs(report['disk_weights']@(report['absolute_current']@v)),
         disk_front_reaction=as_pairs(report['disk_weights']@(report['front_force']@v)),
         modal_residual=float(np.linalg.norm((H+response(z)-z*z*M)@v)/np.linalg.norm(H@v)))
    emit(kind='traces',radius=report['radii'].tolist(),weights=report['disk_weights'].tolist(),
         radii_normalized=(report['radii']/shape.exit).tolist(),**quantities)


def readout(folder):
    import json
    from pathlib import Path
    folder=Path(folder)
    def result(path):
        guard=json.loads((path/'outcome.json').read_text())
        if guard['exitCode'] or guard['childOutcome']['exitCode'] or not guard['limitsVerified']:
            raise ValueError(f'Invalid guarded result at {path}')
        rows=list(map(json.loads,(path/'stdout').read_text().splitlines()))
        if [r['kind'] for r in rows]!=['input','mode','traces']:
            raise ValueError(f'Incomplete result at {path}')
        return rows[:2]
    sets={}
    for path in sorted(folder.iterdir()):
        if (path/'outcome.json').is_file() and path.name.startswith('case_'):
            setup,measured=result(path)
            spec=setup['args']
            key=(spec['mesh'],spec['box'],spec['pml'],spec['stretch'],spec['front'])
            sets.setdefault(key,{})[spec['offset']]=measured
    shape,_,_,scale,ext,rows,_=reference()
    for case, by_offset in sets.items():
        if 0. not in by_offset:
            continue
        center=by_offset[0.]
        w=center['pole'][0]
        ev=np.array(center['mass_unit_vector'])
        v=ev[:,0]+1j*ev[:,1]
        for h in (.006,.012):
            if -h not in by_offset or h not in by_offset:
                continue
            grad=-(by_offset[h]['pole'][0]-by_offset[-h]['pole'][0])/(2*h*shape.exit)
            load,=[r['current_closing_load'] for r in rows if r['kind']=='signed' and r['index']==0 and r['exit_step']==h]
            def array(pairs):
                x=np.array(pairs)
                return x[...,0]+1j*x[...,1]
            kb=(array(by_offset[h]['anchor_stiffness'])-array(by_offset[-h]['anchor_stiffness']))/(2*h*shape.exit)
            kw=array(center['anchor_slope'])
            reactive_force=np.vdot(v,kb@v).real/np.vdot(v,kw@v).real
            record=dict(kind='signed_exit',reactive_force_at_anchor=reactive_force,case=case,signed_fraction=h,
                        real_pole_frequency_slope_force=grad,closing_load=load,
                        force_ratio_to_locally_loaded=grad/1.07730,loss_fraction=center['pole_cycle_loss'])
            if grad>0:
                J=load/grad
                amp=np.sqrt(J*w/center['peak_unit_stored_energy'])
                phase=np.linspace(0.,np.pi,129)
                frames=[displacement(shape,ext,scale,amp*(v*np.exp(1j*p)).real,w,w/2) for p in phase]
                record.update(formal_action=J,peak_coordinate=as_pairs(amp*v),
                              wall_displacement=max(f['wall']['peak_displacement'] for f in frames),
                              wall_strain=max(f['wall']['peak_principal_strain'] for f in frames),
                              collar_strain=max(f['outermost_registered']['peak_principal_strain'] for f in frames))
            emit(**record)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh',type=int,default=9)
    parser.add_argument('--box',type=float,default=5.)
    parser.add_argument('--pml',type=float,default=3.)
    parser.add_argument('--stretch',type=float,default=2.1)
    parser.add_argument('--offset',type=float,choices=(-.012,-.006,0.,.006,.012),default=0.)
    parser.add_argument('--front',choices=('mean','thin'),default='mean')
    parser.add_argument('--reduce')
    args=parser.parse_args()
    if args.reduce: readout(args.reduce)
    else: run(args)
