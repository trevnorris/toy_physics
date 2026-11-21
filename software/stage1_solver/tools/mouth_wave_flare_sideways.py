"""Quasistatic stiffness and low-frequency inertia of a trapped real scalar.

Translate the centers of the entrance and the distal flare independently at
R=H=1. Read carrier loads from the two existing material-shape tests. No bare
wall mass or imposed sideways spring. Wave action is an independent input.
"""
import argparse
import json
import numpy as np
import yaml
from scipy import sparse
from scipy.sparse.linalg import splu

from mouth_wave_flare_work import ARCHIVE, RECORD, Wave, emit, taper


def run(step, padding):
    params = yaml.safe_load(ARCHIVE.read_text())['parameters']
    wave = Wave(1.,1.,step,padding,params)
    omega, u, op, _, radial_gate = wave.ground(1.)
    gate = wave.gate(1.)[0]
    br = params['delta_r']/2
    vr = wave.contrast*gate*(1-radial_gate)/br
    source = vr.ravel()*u/np.sqrt(3)
    radius = np.broadcast_to(wave.r,gate.shape).ravel()
    dipole = op + sparse.diags(2*params['c_w']**2/radius**2-omega**2,format='csc')
    inverse = splu(dipole).solve
    full = inverse(source)
    # A uniformly translated, prescribed potential has a separate solution
    # by a phase boost and a reduced kinetic coefficient along its velocity.
    radial_u = np.pad(u.reshape(gate.shape), ((1,1),(0,0)))
    rigid_mass = omega/params['c_w']**2 + np.sum(np.diff(radial_u,axis=0)**2)/(3*omega*step**2)
    reciprocal = np.array([[1.,-1.],[-1.,1.]])
    run_input = json.loads(RECORD.read_text())['runs']['carrier_attachment_40_01']
    if run_input['exit_code'] or run_input['worker_exit'] or not run_input['limits_verified']:
        raise ValueError('Invalid independent carrier attachment data')
    action = params['I_wave']
    emit(test='settings', step=step, padding=padding, action=action, omega=omega)
    for line in run_input['stdout_lines']:
        row = json.loads(line)
        near = np.broadcast_to(taper(wave.w,*row['switch'])[0],gate.shape).ravel()
        q = source[:,None]*np.column_stack((near,1-near))
        first = inverse(q)
        k = ((near**2*source)@full-q[:,0]@first[:,0])/omega
        vxx = (vr*((2*radial_gate-1)/br+2/wave.r)/3).ravel()
        k_direct = ((u*near)**2@vxx/2-q[:,0]@first[:,0])/omega
        mass = (first.T@first+4*omega**2*(first.T@inverse(first)))/omega
        mc = np.array(row['mass_matrix'])
        held = np.array(row['held_per_intake'])
        total = mc+action*mass
        center = total.sum(axis=0)/total.sum()
        emit(test='response', switch=row['switch'], k_per_action=k, direct_k_per_action=k_direct,
             mass_matrix_per_action=mass.tolist(), wave_rigid_mass_per_action=mass.sum(),
             boosted_rigid_mass_per_action=rigid_mass,
             carrier_mass=mc.sum(), common_acceleration=held.sum()/total.sum(),
             limit_distal_force_share=(total.sum(axis=0)*held.sum()/total.sum()-held)[1]/held.sum())
        if k<=0 or np.linalg.eigvalsh(mass)[0]<=0 or abs(mass.sum()/rigid_mass-1)>.01:
            raise ValueError('Stiffness, inertia, or independent rigid boost check failed')
        # Evaluate the sideband resolvents themselves, then subtract the static
        # term before assembling the dynamic wave force. Positive kinetic sign
        # for the carrier gives the usual minus nu^2 M in displacement units.
        empty = np.linalg.solve(mc,held)
        empty_center = mc.sum(axis=0)@empty/mc.sum()
        for nu in (.02,.01,.005):
            response = []
            for sign in (-1,1):
                shift = nu**2+sign*2*omega*nu
                displaced = splu(dipole - shift*sparse.eye(dipole.shape[0],format='csc')).solve(q)
                response.append(q.T@(displaced-first))
            dynamic = k*reciprocal-sum(response)/(2*omega)
            x = np.linalg.solve(action*dynamic-nu**2*mc,held)
            slow = np.linalg.solve(action*k*reciprocal-nu**2*total,held)
            xc = center@x
            emit(test='driven_following', switch=row['switch'], frequency=nu,
                 signed_excursion=(x[1]-x[0])/xc,
                 low_frequency_signed_excursion=(slow[1]-slow[0])/(center@slow),
                 acceleration_per_unit_held_sum=(-nu**2*xc/held.sum()),
                 empty_signed_excursion=(empty[1]-empty[0])/empty_center)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--step',required=True,type=float)
    parser.add_argument('--pad',required=True,type=float)
    args=parser.parse_args()
    run(args.step,args.pad)
