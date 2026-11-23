"""Physical cycle peaks for specified angular preparations of the saved degree-12 pole."""
import argparse
import copy
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import Legendre, leggauss, legval
from scipy.special import beta

from mouth_exit_witness import inputs
from mouth_poly_shallow_release import Trial, emit
from mouth_shallow_puncture_fields import PORT
from mouth_trapping_angular_tail import recorded_seed
from mouth_trapping_modes import Sector


DEGREE = 12
EIGENVALUE = DEGREE*(DEGREE+1)


def angular(cos_theta, state):
    """T_theta, T_phi, 2D_theta T_theta, and D_theta T_phi + D_phi T_theta."""
    c = cos_theta
    s = np.sqrt(np.maximum(0.,1-c*c))
    if state == 'standing_m0':
        p = Legendre.basis(DEGREE)
        norm = np.sqrt((2*DEGREE+1)/(4*np.pi*EIGENVALUE))
        # Choose the overall phase to agree with the exact positive azimuthal twist.
        return (0j*c,norm*s*p.deriv()(c),0j*c,-norm*s*s*p.deriv(2)(c))
    if state != 'traveling_m12':
        raise ValueError(state)
    norm = DEGREE/np.sqrt(EIGENVALUE*2*np.pi*beta(.5,DEGREE+1))
    a = norm*s**(DEGREE-1)
    g = norm*(DEGREE-1)*s**(DEGREE-2)
    return -1j*a,a*c,-2j*g*c,g*(1+c*c)


def angular_groups(theta,state):
    a,b,c,d=angular(np.cos(theta),state)
    return abs(a)**2+abs(b)**2, a*a+b*b, abs(c)**2+abs(d)**2, c*c+d*d


def validate_angular(state):
    c,w=leggauss(48)
    a,b,d,e = angular(c,state)
    inertia=2*np.pi*(w@(abs(a)**2+abs(b)**2))
    strain=2*np.pi*(w@(abs(d)**2+abs(e)**2))
    if abs(inertia-1.)>2e-12 or abs(strain-(EIGENVALUE-2))>2e-10:
        raise ArithmeticError('Angular norm or symmetric strain disagrees with the material pencil')
    return inertia,strain


def observe(op,seed,row,radius,mesh):
    omega,v=seed
    # The physical action here means the energy inside the matching port divided
    # by Re omega. Its exterior is outgoing and no infinite-volume norm is used.
    kin=abs(omega)**2*np.vdot(v,op.M@v).real
    proj=op.port.project@v[op.port.ids]
    potential=np.vdot(v,op.K@v).real-2*PORT**3*np.vdot(proj,proj).real
    stored=.25*(kin+potential)
    A=np.sqrt(omega.real/stored)
    full=np.zeros(len(op.live),complex)
    full[op.live]=v
    cell=full[op.ids]*A
    nr=op.mesh
    nt=len(cell)//nr
    ii,jj=np.indices((nr,nt)).reshape(2,-1)
    sample={}
    bounds={}
    groups=('shear_bearing_cavity','registered_collar','collar_order_at_least_half')
    statuses={}
    for state in ('standing_m0','traveling_m12'):
        statuses[state]=validate_angular(state)
        theta=np.linspace(0,np.pi/2,385)
        tnorm,tcomplex,gnorm,gcomplex=angular_groups(theta,state)
        # The physical speed of a real vector harmonic has this exact phase maximum.
        angular_speed=np.sqrt(.5*(tnorm+abs(tcomplex)))
        vel_theta=np.argmax(angular_speed)
        groups_data={g:{'speed':(0.,None),'strain':(0.,None)} for g in groups}
        lower={g:{'speed':0.,'strain':0.} for g in groups}
        for a in (0.,.5,1.):
            x=(ii+a)/nr
            for b in (0.,.5,1.):
                t=(jj+b)/nt
                width,slope,_=op.shape.profile(t)
                r=x*width
                f=(1-b)*((1-a)*cell[:,0]+a*cell[:,1])+b*((1-a)*cell[:,2]+a*cell[:,3])
                fx=nr*((1-b)*(cell[:,1]-cell[:,0])+b*(cell[:,3]-cell[:,2]))
                ft=nt*((1-a)*(cell[:,2]-cell[:,0])+a*(cell[:,3]-cell[:,1]))
                # Components r∂_r f and r∂_w f evaluated at fixed physical r and w.
                R=x*fx
                W=r*ft/op.shape.length-x*x*slope*fx
                uv=r*abs(omega*f)
                mabs=abs(R)**2+abs(W)**2
                msq=R*R+W*W
                fa=abs(f)**2
                fs=f*f
                shear_upper=np.sqrt(tnorm.max()*mabs+gnorm.max()*fa)
                # Use the actual wave element's averaged order in the FE pencil.
                active=op.weights>0
                collar=x>=np.exp(-.5*np.exp(legval(2*t-1,op.coeff)))
                areas=(active,active & collar,active & collar & (op.weights>=.5))
                def where(k,angle):
                    return dict(relative_depth=float(t[k]),relative_radius=float(x[k]),
                                physical_radius=float(r[k]),theta=float(angle),cell_order=float(op.weights[k]))
                for region,take in zip(groups,areas):
                    if not take.any():
                        continue
                    record=groups_data[region]
                    k=np.argmax(np.where(take,uv,-np.inf))
                    speed=uv[k]*angular_speed[vel_theta]
                    if speed>record['speed'][0]:
                        record['speed']=speed,where(k,theta[vel_theta])
                    lo=lower[region]
                    lo['speed']=max(lo['speed'],float(np.max(uv[take])/np.sqrt(8*np.pi)))
                    lo['strain']=max(lo['strain'],float(np.sqrt(np.max((mabs+(EIGENVALUE-2)*fa)[take])/(8*np.pi))))
                    # At a given point maximize over time analytically, then sample
                    # polar angle. Skip only cells whose exact angular upper bound at
                    # the sampled angles cannot beat an already measured lower bound.
                    targets=np.flatnonzero(take & (shear_upper>record['strain'][0]))
                    for start in range(0,len(targets),256):
                        ids=targets[start:start+256]
                        gmax2=.5*(mabs[ids,None]*tnorm+fa[ids,None]*gnorm +
                                  abs(msq[ids,None]*tcomplex+fs[ids,None]*gcomplex))
                        k,j=np.unravel_index(np.argmax(gmax2),gmax2.shape)
                        strain=np.sqrt(gmax2[k,j])
                        if strain>record['strain'][0]:
                            record['strain']=strain,where(ids[k],theta[j])
        sample[state]=groups_data
        bounds=lower
    emit(kind='amplitude',radius=radius,wave_mesh=mesh,material_mesh=72,degree=DEGREE,
         omega=omega.real,unit_action=1.,pole_energy_loss=row['loss'],unit_action_amplitude=A,
         unscaled_mean_kinetic=.25*kin,unscaled_mean_potential=.25*potential,
         action_check=A*A*stored/omega.real,cavity_fraction=row['cavity_fraction'],
         angular_integrals=statuses,angular_independent_lower_bounds=bounds,sampled_peaks=sample,
         samples='3x3 one-sided in each Q1 cavity cell, 385 polar angles over one reflected half sphere; cycle maxima analytic. Entries are lower bounds on spatial suprema.')


def main(radius,mesh):
    n=72
    tr=Trial(radius,n,bonded_core=True)
    case=f'r{radius}n{n}'
    base=Path('_scratch/s11c/mouth_exit_release01')/case/'endpoint.npz'
    sh,c,_,digest,status,*_=inputs(base,radius,n,tr)
    archive=Path('_scratch/s11c/mouth_trapping_tail_mesh01')/f'r{radius}'
    guard=json.loads((archive/'outcome.json').read_bytes())
    raw=(archive/'stdout').read_bytes()
    records=list(map(json.loads,raw.splitlines()))
    if (not guard['limitsVerified'] or guard['exitCode'] or guard['childOutcome']['exitCode']
            or status!='balanced' or records[0]['endpoint_sha256']!=digest):
        raise ValueError('Nonmatching or invalid source')
    if mesh==20:
        parent=Path('_scratch/s11c/mouth_trapping_tail01')/case
        raw=(parent/'stdout').read_bytes()
        if hashlib.sha256(raw).hexdigest()!=records[0]['tail_stdout_sha256']:
            raise ValueError('Changed angular source')
        target,=[r for r in map(json.loads,raw.splitlines()) if r['kind']=='angular_mode' and r['degree']==DEGREE]
    else:
        target,=[r for r in records if r['kind']=='wave_refinement' and r['wave_mesh']==mesh]
    material=copy.copy(tr)
    material.mesh=mesh
    op=Sector(material,sh,c,DEGREE)
    seed,_=recorded_seed(op,target)
    emit(kind='inputs',radius=radius,wave_mesh=mesh,endpoint_sha256=digest,source_sha256=hashlib.sha256(raw).hexdigest())
    observe(op,seed,target,radius,mesh)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--radius',type=int,choices=(1,2),required=True)
    p.add_argument('--mesh',type=int,choices=(20,48),required=True)
    a=p.parse_args()
    main(a.radius,a.mesh)
