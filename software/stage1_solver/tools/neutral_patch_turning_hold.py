"""Join an exterior material's angular displacement to that of the held patch.

The angular component is the sole attachment. The neighboring sheet and patch
move together at that angular coordinate. Inner material force, constraint
reaction, and the torque transmitted through this connection are read separately.
"""
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import block_diag, cho_factor, cho_solve, eigh, solve

from neutral_patch_resting_order_hold import LEVELS, coordinates, digest, emit, hessian
from neutral_patch_sites import LAWS, Patch
from neutral_patch_turning_neighbor import annulus_matrices


RECORD = Path('research/pde_ledger_v3/_measurements/neutral_patch_strain_hold.json')
DIRECTORY = Path('_scratch/s11c/patch_resting_order_hold_v1')


def read_patch(law, order, old):
    path = DIRECTORY/'traces'/f'{law}_{order:.2f}.npz'
    if digest(path) != old['trace_sha256']:
        raise ValueError('Prepared patch differs from the prior angle comparison')
    with np.load(path) as data:
        xy = data['rested_xy']
        times = data['time']
        refs = data['kernel']
    potential = Patch(law)
    potential.depth[potential.cross] *= order/.95
    hp = hessian(potential, xy.ravel())
    center = xy-xy.mean(axis=0)
    turn = np.stack((-center[:,1],center[:,0]),axis=1).ravel()
    angle = turn/(turn@turn)
    if not np.isclose(angle@turn, 1):
        raise ArithmeticError('Rigid-angle scale changed')
    return hp, angle, times[-1], refs[-1]


def coupled_response(hp, angle, end_time, order, ratio, law, cells, check_ode=False):
    radii, stiffness, damping = annulus_matrices(order, ratio, law, cells)
    n = len(hp)
    count = len(radii)
    if stiffness.shape != damping.shape or stiffness.shape != (count,count):
        raise ValueError('Annulus matrix dimension mismatch')
    # x contains all microscopic coordinates followed by the interior angles.
    # The inner sheet angle follows the patch and the remote edge is held still.
    trace = np.zeros((count,n+count-1))
    trace[0,:n] = angle
    trace[1:,n:] = np.eye(count-1)
    system = trace.T@stiffness@trace
    system[:n,:n] += hp
    metric = trace.T@damping@trace
    metric[:n,:n] += np.eye(n)
    fixed, inside, translation = coordinates()
    free = block_diag(inside,np.eye(count-1))
    drive = np.r_[translation,np.zeros(count-1)]
    mr = free.T@metric@free
    kr = free.T@system@free
    # A displacement step leaves the unconstrained coordinates free of externally
    # applied impulses. This is the limit of a short ramp in the drag metric.
    first = drive-free@solve(mr, free.T@metric@drive,assume_a='pos')
    rates, modes = eigh(kr,mr)
    if min(rates)<=0 or not np.isfinite(rates).all():
        raise ArithmeticError('Free material modes must decay')
    modes = free@modes
    projection = modes.T@system@first
    sample = np.array([0.,end_time/4,end_time,np.inf])
    decay = np.exp(-sample[:,None]*rates)
    x = first+np.expm1(-sample[:,None]*rates)*(projection/rates)@modes.T
    dx = -(decay*projection)@modes.T
    angles = x@trace.T
    torques = angles@stiffness+dx@trace.T@damping
    force = x[:,:n]@hp@translation
    total = (x@system+dx@metric)@drive
    # The fixture that controls the two layer centers also does work against
    # the moment applied at the patch's angle, and no direct drag after the step.
    virtual_work = total-force-(angle@translation)*torques[:,0]
    means = x[:,:n]@fixed.T
    balance = abs((x@system+dx@metric)@free).max()
    # Condense the sheet at zero frequency and solve the patch by itself.
    k_static = stiffness[0,0]-stiffness[0,1:]@solve(stiffness[1:,1:],stiffness[1:,0],assume_a='pos')
    static = hp+k_static*np.outer(angle,angle)
    end = translation-inside@solve(inside.T@static@inside,inside.T@static@translation,assume_a='pos')
    static_residual = abs(end@hp@translation-force[-1])
    ode_error = None
    if check_ode:
        mass_factor = cho_factor(mr)
        drift = -cho_solve(mass_factor,kr)
        source = -cho_solve(mass_factor,free.T@system@first)
        response = solve_ivp(lambda t,y: drift@y+source,(0,end_time),np.zeros(len(kr)),
                             method='BDF',jac=drift,rtol=2e-9,atol=1e-13)
        if not response.success:
            raise ArithmeticError(response.message)
        direct = first+free@response.y[:,-1]
        ode_error = float(abs(direct[:n]@hp@translation-force[2]))
    errors = [np.max(abs(virtual_work)),balance,np.max(abs(means-fixed@translation)),static_residual,ode_error or 0.]
    if max(errors)>1e-6:
        raise ArithmeticError(f'Coupled force, work, integration or static check failed: {errors}')
    return dict(ratio=ratio,cells=cells,sheet=law,material_force=force.tolist(),
                center_fixture_force=total.tolist(),patch_angle=x[:,:n].dot(angle).tolist(),
                opposing_torque=torques[:,0].tolist(),static_sheet_stiffness=float(k_static),
                bond_force_retained_at_T=float(force[2]/force[0]),static_force=float(force[-1]),
                max_constraint_or_balance_error=float(max(errors)),direct_ode_error=ode_error)


def main():
    prior = json.loads(RECORD.read_text())['guarded_measurements'][str(DIRECTORY/'angle_guarded')]['literal_stdout']
    old = {(r['law'],r['chi']):r for r in map(json.loads,prior) if r['test']=='stationary_angle_control'}
    for patch_law in LAWS:
        for order in LEVELS:
            ref = old[patch_law,order]
            hp, theta, end_time, free_force = read_patch(patch_law,order,ref)
            pinned_force = ref['material_gripped'][-2]
            np.testing.assert_allclose(free_force,ref['material_free'][-2],rtol=1e-9)
            emit('patch_reference', patch=patch_law,order=order,T=end_time,trace_sha256=ref['trace_sha256'],
                 detached_force=free_force,angle_fixed_force=pinned_force,free_initial_force=ref['material_free'][0],
                 fixed_initial_force=ref['material_gripped'][0],detached_angle=ref['free_angle_per_q'][-2])
            for ratio in (2.,4.,8.,16.):
                for law in ('curl','symmetric'):
                    results = []
                    for cells in (48,96):
                        ode = patch_law=='central' and order==.5 and ratio==4 and cells==48
                        results.append(coupled_response(hp,theta,end_time,order,ratio,law,cells,ode))
                    coarse, fine = results
                    err = abs(coarse['material_force'][2]-fine['material_force'][2])/abs(fine['material_force'][2])
                    emit('bonded_material',patch=patch_law,order=order,relative_mesh_change_at_T=err,
                         coarse_ode_error=coarse['direct_ode_error'],**fine)
                    if err>.004:
                        raise ArithmeticError('Coupled material force is spatially underresolved')


if __name__ == '__main__':
    main()
