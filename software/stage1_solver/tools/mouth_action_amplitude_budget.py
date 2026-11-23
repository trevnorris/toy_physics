"""Scale the physical amplitude readings to the matched work threshold and integrate where the kinetic energy sits."""
import argparse
import copy
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss,legval

from mouth_action_amplitude import DEGREE,EIGENVALUE
from mouth_exit_witness import inputs
from mouth_poly_shallow_release import Trial,emit
from mouth_strained_wave import core_factors
from mouth_trapping_angular_tail import recorded_seed
from mouth_trapping_modes import Sector


def verified(path):
    guard=json.loads((path/'outcome.json').read_bytes())
    raw=(path/'stdout').read_bytes()
    if guard['exitCode'] or guard['childOutcome']['exitCode'] or not guard['limitsVerified']:
        raise ValueError(path)
    return list(map(json.loads,raw.splitlines())),hashlib.sha256(raw).hexdigest()


def run(radius,mesh):
    case=f'r{radius}n72'
    src=Path('_scratch/s11c')
    rows,amp_sha=verified(src/'mouth_action_amplitude01'/f'r{radius}m{mesh}')
    source,amp=rows
    work,work_sha=verified(src/'mouth_gallery_action01'/case)
    work=min((r for r in work if r['kind']=='exit_balance'),key=lambda row:row['log_step'])
    action=work['action']
    scale=np.sqrt(action)
    scaled={}
    floor={area:{key:v*scale for key,v in values.items()}
           for area,values in amp['angular_independent_lower_bounds'].items()}
    for state,regions in amp['sampled_peaks'].items():
        scaled[state]={}
        for name,stats in regions.items():
            item={key:values[0]*scale for key,values in stats.items()}
            item['strain_site_order']=stats['strain'][1]['cell_order']
            if state=='standing_m0':
                s=item['strain']
                item['simple_shear_maximum_stretch']=.5*(np.hypot(s,2)+s)
            scaled[state][name]=item
    emit(kind='action_scale',radius=radius,wave_mesh=mesh,action=action,
         action_source_wave_mesh=20,amplitude_stdout_sha256=amp_sha,threshold_stdout_sha256=work_sha,
         phase_and_angle_rms_lower_bounds=floor,sampled_peaks=scaled)
    tr=Trial(radius,72,bonded_core=True)
    sh,c,_,digest,status,*_=inputs(src/'mouth_exit_release01'/case/'endpoint.npz',radius,72,tr)
    if status!='balanced' or digest!=source['endpoint_sha256']:
        raise ValueError('Changed mouth')
    material=copy.copy(tr)
    material.mesh=mesh
    op=Sector(material,sh,c,DEGREE)
    row=dict(frequency=amp['omega'],loss=amp['pole_energy_loss'],cavity_fraction=amp['cavity_fraction'])
    (om,v),_=recorded_seed(op,row)
    full=np.zeros(len(op.live),complex)
    full[op.live]=v
    cell=full[op.ids]*(amp['unit_action_amplitude']*scale)
    nr=mesh;nt=len(cell)//nr
    ii,jj=np.indices((nr,nt)).reshape(2,-1)
    x0=ii/nr;x1=(ii+1)/nr
    chi=op.weights
    sums=dict(total=0.,collar=0.,collar_order_half=0.,order_below_one_percent=0.,
              order_below_one_tenth=0.,order_below_half=0.,meridian_rms_strain_over_one=0.,
              meridian_rms_speed_above_flat=0.,collar_half_rms_strain_over_one=0.)
    maxima={cut:[0.,0.] for cut in (.1,.5,.9)}
    nodes,w=leggauss(3)
    for b,wb in zip((nodes+1)/2,w/2):
        t=(jj+b)/nt
        width,slope,_=op.shape.profile(t)
        p=np.exp(legval(2*t-1,op.coeff))
        mid=np.clip(np.exp(-p/2),x0,x1)
        stretch,_,axial=core_factors(op.shape,op.coeff,op.reference,t)
        for lo,hi,collar in ((x0,mid,False),(mid,x1,True)):
            for eta,wa in zip((nodes+1)/2,w/2):
                x=lo+(hi-lo)*eta
                a=(x-x0)*nr
                r=x*width
                f=(1-b)*((1-a)*cell[:,0]+a*cell[:,1])+b*((1-a)*cell[:,2]+a*cell[:,3])
                fx=nr*((1-b)*(cell[:,1]-cell[:,0])+b*(cell[:,3]-cell[:,2]))
                ft=nt*((1-a)*(cell[:,2]-cell[:,0])+a*(cell[:,3]-cell[:,1]))
                R=x*fx; W=r*ft/op.shape.length-x*x*slope*fx
                rms_strain=np.sqrt((abs(R)**2+abs(W)**2+(EIGENVALUE-2)*abs(f)**2)/(8*np.pi))
                rms_speed=r*abs(om*f)/np.sqrt(8*np.pi)
                if collar:
                    tag=-np.log(np.where(hi>lo,x,1.))/p
                    invj=(op.shape.neck+tag)**2/(op.shape.length*p*np.where(hi>lo,r,1.)**3)
                else:
                    invj=1/(stretch**3*axial)
                amt=wa*wb*(hi-lo)*width*op.shape.length/nt*r**4*invj*chi*abs(f)**2
                sums['total']+=sum(amt)
                if collar:
                    sums['collar']+=sum(amt)
                    sums['collar_order_half']+=sum(amt[chi>=.5])
                    sums['collar_half_rms_strain_over_one']+=sum(amt[(chi>=.5)&(rms_strain>1)])
                    for cut,m in maxima.items():
                        take=(hi>lo)&(chi>=cut)
                        if take.any():
                            m[0]=max(m[0],max(rms_speed[take]))
                            m[1]=max(m[1],max(rms_strain[take]))
                for key,cut in [('order_below_one_percent',.01),('order_below_one_tenth',.1),('order_below_half',.5)]:
                    sums[key]+=sum(amt[chi<cut])
                sums['meridian_rms_strain_over_one']+=sum(amt[rms_strain>1])
                sums['meridian_rms_speed_above_flat']+=sum(amt[rms_speed>1])
    cavity=sums.pop('total')
    check=cavity/(amp['unit_action_amplitude']**2*action)
    if abs(check-amp['cavity_fraction'])>2e-10:
        raise ArithmeticError('Independent squared-amplitude integral lost cavity mass')
    emit(kind='kinetic_location',radius=radius,wave_mesh=mesh,cavity_fraction_check=check,
         fractions_of_cavity_kinetic_energy={key:val/cavity for key,val in sums.items()},
         interior_quadrature_collar_rms_floors_speed_strain=maxima)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--radius',type=int,choices=(1,2),required=True)
    p.add_argument('--mesh',type=int,choices=(20,48),required=True)
    a=p.parse_args();run(a.radius,a.mesh)
