"""Wave-empty tangent breathing modes with dead holding loads on a measured draining shape."""
import argparse
from functools import cache
import hashlib
import json

import numpy as np
from numpy.polynomial.legendre import legval
from scipy.linalg import eigh

from mouth_action_breathing import volume_derivatives
from mouth_exit_witness import inputs
from mouth_live_waist_fields import Carrier
from mouth_poly_bonded_core import core_source
from mouth_poly_shallow_release import BIAS, CUT, ROOT, Order, Trial, emit, exterior, source
from mouth_sealed_shape_relax import Meridian


def displacement(shape, exterior_registration, original_exit, vector, omega, action):
    # Ordered lattice labels r/b(t) are fixed, including the bonded core. For
    # q^T M q=1 the peak real mode coordinate is sqrt(2 E)/omega, E=J omega.
    depth, terminal, middle = np.sqrt(2*action/omega)*vector
    t = np.linspace(0, 1, 2049)
    b, slope, _ = shape.profile(t)
    dr = t*original_exit*terminal + 4*shape.neck*t*(1-t)*middle
    beta = dr/b
    dr_dt = original_exit*terminal + 4*shape.neck*(1-2*t)*middle
    b_beta_prime = dr_dt/shape.length - beta*slope
    dz = t*shape.length*depth
    result = dict(peak_exit_radial_displacement=abs(dr[-1]), peak_depth_displacement=abs(dz[-1]),
                  max_radial_fraction=float(max(abs(beta))), max_axial_fraction=abs(depth))
    # Cylindrical diagonal has three brane components beta and fourth component
    # depth. The off-diagonal r,w component of infinitesimal strain is r beta'/2.
    for name, x in [('wall',1.),('outermost_registered',np.exp(.5*np.exp(legval(2*t-1,exterior_registration))))]:
        grad = x*b_beta_prime
        mean=(beta+depth)/2
        half=np.hypot(beta-depth,grad)/2
        strains=np.maximum(abs(mean-half),abs(mean+half))
        result[name] = dict(peak_displacement=float(np.hypot(x*dr,dz).max()),
                            peak_principal_strain=float(strains.max()),
                            peak_engineering_strain_norm=float(np.sqrt(2*(3*beta**2+depth**2)+grad**2).max()),
                            peak_shear_depth_fraction=float(t[strains.argmax()]))
    return result


def run(radius, mesh):
    folder=ROOT/'_scratch/s11c/mouth_action_breathing01/run'
    evidence=ROOT/'research/pde_ledger_v3/_measurements/mouth_gallery_action.json'
    old=json.loads(evidence.read_bytes())['runs']['mouth_action_breathing01/run']
    recorded=(folder/'stdout').read_bytes()
    if hashlib.sha256(recorded).hexdigest()!=old['stdout_sha256'] or old['guard_exit'] or old['child_exit']:
        raise ValueError('Previously measured lattice inertia failed validation')
    rows=[d for d in map(json.loads,recorded.splitlines()) if d['radius']==radius and d['n']==mesh]
    setup,=[r for r in rows if r['kind']=='setup']
    inertia={r['relative_exit']:np.array(r['mass']) for r in rows if r['kind']=='basis_state'}
    trial=Trial(radius,mesh,bonded_core=True)
    path=ROOT/f'_scratch/s11c/mouth_exit_release01/r{radius}n{mesh}/endpoint.npz'
    sh,chi,_,digest,status,*_=inputs(path,radius,mesh,trial)
    records=list(map(json.loads,(path.parent/'stdout').read_bytes().splitlines()))
    if status!='balanced' or digest!=setup['endpoint_sha256']:
        raise ValueError('Prepared geometry and measured inertia differ')
    init=records[0]['exit_radius']
    shape=Meridian(radius,trial.x+sh[:4])
    b=shape.exit
    co=trial.co.copy(); co[0]+=sh[4]
    cr=trial.cr.copy(); cr[0]+=sh[5]
    depth=np.linspace(0,1,2*mesh+1)
    radial=np.linspace(0,1,mesh+1)

    @cache
    def energy(dL,db,dm):
        # All three coordinates change material positions affinely. Otherwise a
        # nonlinear chart and a nonzero held-shape force change the tangent.
        config=sh[:4].copy()
        config[0]+=np.log1p(dL)
        config[1]+=np.log1p(init*db/b)
        config[2]+=dm
        moved=Meridian(radius,trial.x+config)
        inner=source(moved,co,depth,radial,False)
        inner+=core_source(moved,co,*trial.core_reference,depth,radial)
        outer=exterior(moved,mesh,CUT,side_order_work=lambda w,x:
                       source(moved,cr,w/moved.length,x,True))
        carrier=Carrier(moved,mesh,2*mesh,trial.mesh)
        eflow,_,data=carrier.response(chi,BIAS)
        if abs(data['terminal_current']-data['current'])>1e-6*abs(data['current']):
            raise ArithmeticError('Stationary carrier currents disagree')
        return (outer['energy']+outer['static_work']+np.sum(inner*chi)
                +Order(moved,mesh).energy_gradient(chi)[0]+eflow)

    reference=sum(value for key,value in records[-1]['energy'].items() if key!='strained_wave')
    ebase=energy(0.,0.,0.)
    if abs(ebase-reference)>1e-8:
        raise ArithmeticError('Independent wave-empty energy reconstruction failed')
    emit(kind='input',radius=radius,mesh=mesh,material_sha256=digest,mass_sha256=old['stdout_sha256'],
         b=b,base_depth=shape.length,original_exit_scale=init,center_material_energy=ebase,
         energy_disagreement=ebase-reference)

    def system(exit_fraction,step):
        center=np.array([0.,exit_fraction*b/init,0.])
        e0=energy(*center)
        basis=step*np.eye(3)
        H=np.zeros((3,3))
        g=np.empty(3)
        for i,u in enumerate(basis):
            ep,em=energy(*(center+u)),energy(*(center-u))
            H[i,i]=(ep-e0+em-e0)/step**2
            g[i]=(ep-em)/(2*step)
            for j in range(i):
                v=basis[j]
                H[i,j]=H[j,i]=((energy(*(center+u+v))-energy(*(center+u-v)))
                                -(energy(*(center-u+v))-energy(*(center-u-v))))/(4*step**2)
        M=inertia[exit_fraction]
        eigen,vec=eigh(H,M)
        emit(kind='spectrum',relative_exit=exit_fraction,shape_step=step,hessian=H.tolist(),mass=M.tolist(),
             generalized_gradient=g.tolist(),omega_squared=eigen.tolist(),mass_normalized_vectors=vec.T.tolist(),
             frequency=np.sqrt(np.maximum(eigen,0)).tolist())
        return eigen,vec,g

    for exit_step,h in ((.012,.002),(.006,.001)):
        mid,states,g=system(0.,h)
        left,vleft,_=system(-exit_step,h)
        right,vright,_=system(exit_step,h)
        for index,omega2 in enumerate(mid):
            if min(omega2,left[index],right[index])<=0:
                emit(kind='nonpositive_stiffness',index=index,exit_step=exit_step)
                continue
            vector=states[:,index]
            overlaps=[]
            for offset,side in ((-exit_step,vleft),(exit_step,vright)):
                weight=inertia[offset]
                ov=np.abs(vector@weight@side)/np.sqrt(vector@weight@vector)
                if np.argmax(ov)!=index or ov[index]<.94:
                    raise ArithmeticError('The signed mode cannot be tracked')
                overlaps.append(ov[index])
            force=-(np.sqrt(right[index])-np.sqrt(left[index]))/(2*b*exit_step)
            load=(energy(0.,exit_step*b/init,0.)-energy(0.,-exit_step*b/init,0.))/(2*b*exit_step)
            omega=np.sqrt(omega2)
            report=dict(kind='signed',index=index,exit_step=exit_step,omega=omega,outward_force_per_action=force,
                        minimum_mass_overlap=min(overlaps),current_closing_load=load,
                        volume_change_at_unit_modal_coordinate=float(volume_derivatives(shape,init)@vector))
            if force>0:
                action=load/force
                report.update(conditional_action=action,
                              displacement=displacement(shape,cr,init,vector,omega,action))
            emit(**report)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('radius',choices=(1,2),type=int)
    parser.add_argument('mesh',choices=(48,72),type=int)
    arguments=parser.parse_args()
    run(arguments.radius,arguments.mesh)
