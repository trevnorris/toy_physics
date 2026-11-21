"""Does identical support on the prescribed first bend fix the free entrance response?"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

from mouth_carrier_lateral_response import prepare, weighted_gram
from mouth_free_bend_drive import Response, real_pair
from mouth_free_bend_held_flow import volume_reaction
from mouth_free_bend_incident import boundary_momentum
from mouth_free_bend_inertia import entrance_force
from mouth_free_bend_release import pencil
from mouth_order_wave_pulse import digest, emit

OLD=Path('research/pde_ledger_v3/_measurements/mouth_free_bend_response.json')


def calculation(args):
    dom, field = prepare(args)
    response=Response(dom,field,args.extension)
    energy,time=pencil(response,0.)
    free=response.dipole.free
    mass=weighted_gram(dom.N,dom.weight)
    field_skew=sparse.bmat([[None,-mass],[mass,None]],format='csr')
    boundary_field=1j*dom.points[:,0]*field*dom.feed
    boundary=real_pair(boundary_field)
    flux=volume_reaction(dom,field,boundary_field,args.extension)
    momentum=boundary_momentum(dom,field,boundary_field,args.extension)
    bulk=-(2/3)*(response.dipole.operator@boundary)[free]
    bulk_time=(2/3)*(field_skew@boundary)[free]
    solve=response.dipole.factor.solve
    tangent=solve(response.bend)
    velocity=response.motion+response.skew@tangent
    inertia=(2/3)*velocity.T@solve(velocity)
    cross=response.motion.T@tangent
    gyro=(2/3)*(cross-cross.T-tangent.T@(response.skew@tangent))
    static=response.force_matrix(0.)
    np.testing.assert_allclose(gyro,-gyro.T,atol=1e-8)
    omega=.005
    z=response.force_matrix(1j*omega)
    applied=np.zeros(energy.shape[0])
    applied[-3]=1.
    emit(test='input',archive_sha256=digest(OLD),omega=omega,area_per_tension=response.area.diagonal().tolist(),
         active_coordinates=['entrance translation','first bend','deeper bend'])
    outputs={}
    for first,deeper in ((.02,.02),(.02,.04),(.04,.02),(.04,.04)):
        # Supports do not overlap in depth. Every static first-bend contact
        # experiment holds both prescribed coordinates and feeds the same sites.
        support=np.array([0.,first,deeper])*response.area.diagonal()
        augmented=energy+sparse.diags(np.r_[np.zeros(energy.shape[0]-3),support],format='csc')
        active=np.ones(energy.shape[0],bool)
        active[-3]=False
        rhs0=np.r_[bulk,flux]
        state=np.zeros_like(rhs0)
        state[active]=splu(augmented[active][:,active].tocsc()).solve(rhs0[active])
        held=float((rhs0-augmented@state)[-3])
        k=static+np.diag(support)
        v=np.r_[1.,-np.linalg.solve(k[1:,1:],k[1:,0])]
        y=(gyro@v)[1:]
        effective_mass=v@inertia@v+y@np.linalg.solve(k[1:,1:],y)
        op=augmented-1j*omega*time
        loads=np.column_stack((np.r_[bulk+1j*omega*bulk_time,flux-1j*omega*momentum],applied))
        dynamic=splu(op.tocsc()).solve(loads)
        residual=float(np.max(abs(op@dynamic-loads)))
        local_z=z+np.diag(support)
        np.testing.assert_allclose(dynamic[-3:,1],np.linalg.solve(local_z,[1.,0.,0.]),rtol=2e-7,atol=1e-8)
        if residual>3e-8 or min(np.linalg.eigvalsh(k))<=0:raise ArithmeticError('Free-response validation failed')
        observed=-omega**2*dynamic[-3,:]
        value=float(observed[0].real)
        outputs[first,deeper]=value
        emit(test='release',first=first,deeper=deeper,positive_step_radius=5e-6,
             first_bend_area_holder=float(support[1]*5e-6),second_bend_area_holder=float(support[2]*5e-6),
             static_pull_per_incident_speed=held,effective_mass=effective_mass,bath_spring=entrance_force(k),
             rigid_mass=float(inertia[0,0]),bend_flow_mass=float(y@np.linalg.solve(k[1:,1:],y)),
             pull_over_mass=held/effective_mass,free_entrance_per_speed=[value,float(observed[0].imag)],
             entrance_per_external_push=[float(observed[1].real),float(observed[1].imag)],balance=residual)
    for a in (.02,.04):
        # Independent, earlier homogeneous-tension builds solve the full free system.
        archive=json.loads(OLD.read_text())
        candidates=[]
        def visit(item):
            if isinstance(item,dict):
                if item.get('test')=='free_incident' and item.get('tension')==a and item.get('omega')==omega:
                    candidates.append(item['entrance_and_flare_real'][0])
                for v in item.values(): visit(v)
            elif isinstance(item,list):
                for v in item: visit(v)
            elif isinstance(item,str) and item.startswith('{'):
                try: parsed=json.loads(item)
                except ValueError:return
                if isinstance(parsed,dict):visit(parsed)
        visit(archive)
        if not candidates or min(abs(v-outputs[a,a]) for v in candidates)>1e-8:
            raise ArithmeticError('Homogeneous result failed archived comparison')
    delta=outputs[.02,.04]/outputs[.02,.02]-1
    emit(test='same_first_bend',fractional_change_free_acceleration=delta)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=Path('_scratch/s11c/mouth_flared_cubic_n34'))
    p.add_argument('--extension',choices=['wide','compact'],default='wide')
    calculation(p.parse_args())
