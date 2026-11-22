"""Finite observation of tracers on a half-order local surface.

Axial entrants follow the signed flux measure in a one-pitch band. The offset
of the atomic plane, radial velocity dispersion, and initially cold in-plane
host are controlled preparations. The impact beam instead conditions on an
incident wall-normal crossing. It has a different denominator.
"""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
from scipy.special import ndtri
from scipy.stats import qmc

from mouth_attachment_incidence import HalfOrderBand, incident_normal_speed, thermal_at_fixed_axial
from mouth_segment_contact_hold import PITCH
from neutral_patch_sites import TEMPERATURE, emit


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(args):
    band=HalfOrderBand.load(args.section)
    random=qmc.Sobol(5,scramble=True,seed=args.seed).random_base2(args.power)
    registry=random[:,:2]
    thermal=np.sqrt(TEMPERATURE)
    clearance,normal_draw,yaw=random[:,2:].T
    yaw=2*np.pi*yaw
    cs,sn=np.cos(yaw),np.sin(yaw)
    directions=[(0.,0.)]
    if args.signed:
        directions.extend((drive,cosine) for drive in (.002,-.002,.001,-.001)
                          for cosine in (-1/np.sqrt(3),1/np.sqrt(3)))
    meta=[]
    input_rows=[]
    for drive,cosine in directions:
        arrival=band.launch(drive,cosine,clearance)
        if args.preparation=='radial':
            arrival=thermal_at_fixed_axial(band,arrival,thermal*ndtri(normal_draw))
        tangent=np.array(arrival['tangent'])
        polar=np.array(arrival['polar'])
        normal=np.array(arrival['normal'])
        weight=np.array(arrival['weight'])
        if args.preparation=='beam':
            # Maxwellian flux onto an exposed face, conditioned on approaching.
            z=np.full_like(clearance,args.offset)
            normal=-incident_normal_speed(normal_draw,TEMPERATURE)
            weight=np.ones_like(clearance)
            gate=0.
        else:
            # A transverse kick has zero axial w component. It broadens the
            # normal incidence while retaining the measured axial flux weights.
            z=args.offset+arrival['z_offset']
            gate=args.gate_depth*np.hypot(1,band.wall_slope)/PITCH
        chi=band.scalar_order+drive*cosine*band.order_dipole
        for k in range(len(random)):
            uid=len(meta)
            vx=tangent[k]*cs[k]-polar[k]*sn[k]
            vy=tangent[k]*sn[k]+polar[k]*cs[k]
            row=(uid,*registry[k],z[k],vx,vy,normal[k],chi,cs[k],sn[k],gate)
            input_rows.append(' '.join(str(x) for x in row))
            meta.append(dict(id=uid,sample=k,U=drive,cosine=cosine,weight=float(weight[k]),
                             initial=[float(s) for s in row[1:]]))
    steps=round(args.time/args.dt)
    dt=args.time/steps
    command=[str(args.engine.resolve()),'--side',str(args.side),'--pitch',str(band.pitch),
             '--dt',str(dt),'--steps',str(steps),'--host',args.host,'--mixed',str(args.mixed),
             '--gate-slope',str(band.wall_slope),'--input','-']
    emit(test='preparation',section=str(args.section),section_sha256=sha(args.section),
         engine=str(args.engine),engine_sha256=sha(args.engine),engine_command=command,
         preparation=args.preparation,seed=args.seed,power=args.power,offset=args.offset,
         side=args.side,dt=dt,native_time=args.time,model_time=args.time*PITCH/np.sqrt(band.scale),
         scale=band.scale,pitch=band.pitch,gate_depth=args.gate_depth,thermal_native=TEMPERATURE,
         chi=band.scalar_order,chi_dipole=band.order_dipole,wall_slope=band.wall_slope,
         sample_count=len(random),trajectories=len(meta))
    if args.output.exists():
        raise FileExistsError(args.output)
    result=subprocess.run(command,input='\n'.join(input_rows)+'\n',text=True,stdout=subprocess.PIPE,check=True)
    readings=[json.loads(s) for s in result.stdout.splitlines()]
    if len(readings)!=len(meta):
        raise ArithmeticError('Incomplete trajectory readouts')
    with args.output.open('x') as stream:
        for m,diagnostics in zip(meta,readings,strict=True):
            if diagnostics['id']!=m['id']:
                raise ArithmeticError('Trajectory order changed')
            stream.write(json.dumps(dict(**m,result=diagnostics),allow_nan=False)+'\n')
    emit(test='trajectory_file',path=str(args.output),sha256=sha(args.output),traces=len(meta),
         initial_energy_range=[min(r['energy_tag_initial'] for r in readings),
                               max(r['energy_tag_initial'] for r in readings)],
         max_energy_error=max(r['max_energy_error'] for r in readings),
         max_momentum_error=max(r['max_momentum_error'] for r in readings))
    if min(r['energy_tag_initial'] for r in readings)<=0:
        raise ArithmeticError('An incoming tag was already energetically bound')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--engine',type=Path,required=True)
    p.add_argument('--section',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--power',type=int,required=True)
    p.add_argument('--seed',type=int,required=True)
    p.add_argument('--host',choices=('mobile','rigid'),default='mobile')
    p.add_argument('--mixed',choices=(0,1),type=int,default=1)
    p.add_argument('--side',type=int,default=5)
    p.add_argument('--offset',type=float,default=2.)
    p.add_argument('--gate-depth',type=float,default=.01)
    p.add_argument('--time',type=float,default=8.)
    p.add_argument('--dt',type=float,default=.002)
    p.add_argument('--signed',action='store_true')
    p.add_argument('--preparation',choices=('cold','radial','beam'),required=True)
    run(p.parse_args())
