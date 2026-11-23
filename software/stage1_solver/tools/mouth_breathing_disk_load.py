"""Load the wave-empty breathing trial with the measured normal traction of a local translating planar front."""
import argparse
from functools import cache
import hashlib
import json

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import newton

from mouth_breathing_matched import displacement
from mouth_breathing_mobile_patch import preparation, response
from mouth_exit_witness import inputs
from mouth_live_waist_fields import Carrier
from mouth_poly_shallow_release import BIAS, ROOT, Trial, emit
from mouth_sealed_shape_relax import Meridian


def records(path):
    g = json.loads((path/'outcome.json').read_bytes())
    if not g['limitsVerified'] or g['exitCode'] or g['childOutcome']['exitCode']:
        raise ValueError('Unqualified source '+str(path))
    raw = (path/'stdout').read_bytes()
    return list(map(json.loads,raw.splitlines())),hashlib.sha256(raw).hexdigest()


def mode_with_disk(H, M, disk, j, archived):
    """Real-axis conservative root and outgoing complex pole via a local analytic continuation."""
    vacuum, _ = eigh(H,M)
    bare = np.sqrt(vacuum[0])
    @cache
    def sample(w):
        if archived and abs(w-archived['frequency'])<1e-12 and abs(j-archived['normal_current'])<1e-11:
            r = archived
        else:
            r = response(j,w,21.,6)
        if max(r['acoustic_work_relative_error'],0)>2e-5 or min(r['toward_feed'],r['toward_receiver'])<=0:
            raise ArithmeticError('Patch power balance failed')
        return complex(r['force_re'],r['force_im']),r

    def equation(w):
        f, _ = sample(float(w))
        ev = eigh(H+f.real*disk,M,eigvals_only=True)
        return ev[0]-w**2

    # The patch has zero static translation stiffness in this homogeneous planar law.
    w = float(newton(equation,bare,x1=bare*.998,tol=3e-12,maxiter=20))
    f,r = sample(w)
    vals,vecs = eigh(H+f.real*disk,M)
    vec = vecs[:,0]
    weight = vec@disk@vec
    controls = []
    for step in (.008,.004):
        fminus = sample(w*(1-step))[0]
        fplus = sample(w*(1+step))[0]
        df = (fplus-fminus)/(2*w*step)
        quadratic = (fplus-2*f+fminus)/(2*(w*step)**2)
        effective_mass = vec@M@vec-weight*df.real/(2*w)
        width = -weight*f.imag/(w*effective_mass)
        def det(z):
            dw = z-w
            return np.linalg.det(H-z*z*M+(f+df*dw+quadratic*dw**2)*disk)
        z = newton(det,w-1j*.5*width,x1=w-1j*.49*width,tol=2e-13,maxiter=20)
        energy=.5*w*w*effective_mass
        P = weight*r['power_per_unit_displacement_squared']
        Q = w*energy/P
        if (effective_mass<=0 or z.imag>=0 or abs((z.real/(-2*z.imag))/Q-1)>1e-3
                or abs(w**2-vals[0])>1e-7 or abs(width)/(w*step)>.1):
            raise ArithmeticError('Unresolved reactive dispersion or continuation')
        controls.append(dict(fractional_step=step,effective_mass=effective_mass,
                             pole=[z.real,z.imag],weak_loss_width=width,
                             Q=Q,stored_energy_per_amplitude_squared=energy,
                             mean_power_per_amplitude_squared=P,disk_area_motion_weight=weight,
                             feed_power=weight*r['toward_feed'],receiver_power=weight*r['toward_receiver'],
                             analytic_force_frequency_derivative=[df.real,df.imag],
                             translation_order_peak_per_displacement=r['peak_linear_order_per_unit_displacement']))
    a,b=controls
    if abs(a['pole'][0]-b['pole'][0])>1e-8 or abs(a['pole'][1]/b['pole'][1]-1)>1e-4:
        raise ArithmeticError('Pole derivative control did not converge')
    return bare,w,vec,controls


def run(radius):
    n=72
    scratch=ROOT/'_scratch/s11c'
    source='r1n72_direction' if radius==1 else 'r2n72'
    rows,digest=records(scratch/'mouth_breathing_matched01'/source)
    archived,patch_hash=records(scratch/'mouth_breathing_mobile_patch01'/f'matched_r{radius}')
    target=next(r for r in archived if r['kind']=='two_sided_patch' and r['normal_current']>0 and r['mesh_refinement']==6)
    case,_,_=preparation(radius)
    sha=rows[0]['material_sha256']
    if sha!=case['endpoint_sha256'] or archived[0]['source_sha']!=digest:
        raise ValueError('Patch and shape background mismatch')
    trial=Trial(radius,n,bonded_core=True)
    loc=scratch/f'mouth_exit_release01/r{radius}n{n}/endpoint.npz'
    sh,chi,_,actual,status,*_=inputs(loc,radius,n,trial)
    if status!='balanced' or sha!=actual:
        raise ValueError('Failed material fingerprint')
    b0=case['physical_exit']
    init=rows[0]['original_exit_scale']
    cr=trial.cr.copy()
    cr[0]+=sh[5]
    emit(kind='setup',radius=radius,material_sha256=sha,shape_source_sha256=digest,
         patch_source_sha256=patch_hash,full_three_coordinate=radius==2)
    central_rows=[r for r in rows if r['kind']=='spectrum' and r['relative_exit']==0]
    central=min(central_rows,key=lambda row:row['shape_step'])
    # R1 has signed Rayleigh data for one central mass-normalized shape vector.
    r1v=np.array(central['mass_normalized_vectors'][0])
    @cache
    def data(offset):
        config=trial.x+sh[:4]
        config[1]+=np.log1p(offset)
        shape=Meridian(radius,config)
        area=4*np.pi*shape.exit**3/3
        if offset==0:
            current=case['current']
        else:
            _,_,d=Carrier(shape,n,2*n,trial.mesh).response(chi,BIAS)
            current=d['current']
            if abs(current-d['terminal_current'])>1e-6*abs(current):
                raise ArithmeticError('Lateral leakage in saved carrier')
        return shape,current/area,area

    for h,step in ((.012,.002),(.006,.001)):
        mode=next(r for r in rows if r['kind']=='signed' and r['index']==0 and r['exit_step']==h)
        poles,carrier=[],[]
        for offset in (-h,0.,h):
            shape,j,area=data(offset)
            if radius==2:
                row=next(r for r in rows if r['kind']=='spectrum' and r['relative_exit']==offset and r['shape_step']==step)
                H=np.array(row['hessian']); M=np.array(row['mass'])
                transform=np.eye(3)
            else:
                transform=r1v[:,None]
                if offset==0:
                    M=transform.T@np.array(central['mass'])@transform
                    H=transform.T@np.array(central['hessian'])@transform
                else:
                    read=next(r for r in rows if r['kind']=='work' and r['index']==0 and
                              r['relative_exit']==offset and r['step']==step)
                    # Directional work normalized its central eigenvector to max |component|=1.
                    scale=np.max(abs(r1v))**2
                    M=np.array([[read['mass']*scale]]); H=np.array([[read['curvature']*scale]])
            disk=area*(shape.length**2)*np.outer(transform[0],transform[0])
            bare,w,vector,checks=mode_with_disk(H,M,disk,j,target if offset==0 else None)
            physical=transform@vector
            norm=checks[-1]['effective_mass']
            poles.append(complex(*checks[-1]['pole']))
            carrier.append(bare)
            emit(kind='local_loading',exit_step=h,offset=offset,area=area,mean_normal_current=j,
                 total_current=j*area,axial_depth=shape.length,bare=bare,real_axis_root=w,
                 unit_lattice_mass_mode=physical.tolist(),checks=checks)
            if offset==0:
                center=shape,w,physical,norm,checks[-1]
        shape,w,v,eff,check=center
        force=-(poles[-1].real-poles[0].real)/(2*b0*h)
        bare_force=-(carrier[-1]-carrier[0])/(2*b0*h)
        if abs(bare_force-mode['outward_force_per_action'])>1e-7:
            raise ArithmeticError('Removed-load derivative failed to recover the saved force')
        if force<=0:
            emit(kind='result',step=h,outward_force=force,bare_force=bare_force)
            continue
        action=mode['current_closing_load']/force
        ampl=np.sqrt(2*action/(w*eff))
        depth_motion=abs(ampl*shape.length*v[0])
        Q=poles[1].real/(-2*poles[1].imag)
        emit(kind='result',step=h,outward_force=force,bare_force=bare_force,closing=mode['current_closing_load'],
             action=action,Q=Q,energy_fraction_lost_per_cycle=-np.expm1(-2*np.pi/Q),
             total_decay_time=1/(-2*poles[1].imag),
             displacement=displacement(shape,cr,init,v,w,action/eff),
             eulerian_patch_order_linear_at_load=check['translation_order_peak_per_displacement']*depth_motion)


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('radius',type=int,choices=(1,2))
    run(cli.parse_args().radius)
