"""Fix the linearized whole-patch Procrustes angle during a cold held step.

The rotation generator r=J(X-X_com) includes both layer centroids and all
within-layer offsets. Its linear angle is r.du/(r.r), in radians. At fixed
angle the extra grip does zero angular work. The material force on the lower
layer and the force conjugate to q at fixed angle are separately recorded.
"""
import argparse
from pathlib import Path

import numpy as np
from scipy.linalg import null_space

from neutral_patch_resting_order_hold import LEVELS, coordinates, digest, emit, hessian
from neutral_patch_sites import COUNT, LAWS, Patch


def held_displacement(full, free, drive, times):
    rates, vectors = np.linalg.eigh(free.T@full@free)
    if min(rates) <= 0.:
        raise ValueError('Nonpositive curvature in the internal coordinate space')
    modes = free@vectors
    amplitudes = (modes.T@full@drive)/rates
    fraction = -np.expm1(-times[:, None]*rates)
    return drive-(fraction*amplitudes)@modes.T


def read_control(path, law, level):
    with np.load(path) as archive:
        xy = archive['rested_xy']
        t = archive['time']
        kernel = archive['kernel']
    patch = Patch(law)
    patch.depth[patch.cross] *= level/.95
    full = hessian(patch, xy.ravel())
    center = xy-xy.mean(axis=0)
    rot = np.stack((-center[:, 1], center[:, 0]), axis=1).ravel()
    inertia = rot@rot
    fixed, free, drive = coordinates()
    rot_free = free@(free.T@rot)
    # Independent angular actuator: one radian, the two x means and common
    # y origin fixed. The pure q motion at fixed angle compensates its own
    # small Procrustes rotation. Both columns are orthogonal to the remaining
    # internal motions, preserving unit mobility under constrained drift.
    angle_actuator = rot_free*(inertia/(rot_free@rot_free))
    affine_angle = (rot@drive)/inertia
    grip_drive = drive-affine_angle*angle_actuator
    grip_free = free@null_space((rot@free)[None, :])
    idx = [0, len(t)//64, len(t)//16, len(t)//4, len(t)-1]
    times = np.r_[t[idx], np.inf]
    loose = held_displacement(full, free, drive, times)
    gripped = held_displacement(full, grip_free, grip_drive, times)
    loose_grad = loose@full
    grip_grad = gripped@full
    loose_bond = loose_grad@drive
    grip_bond = grip_grad@drive
    torque = grip_grad@angle_actuator
    conjugate = grip_grad@grip_drive
    work_residual = max(abs(conjugate-(grip_bond-affine_angle*torque)))
    angle_gripped = gripped@rot/inertia
    angle_loose = loose@rot/inertia
    loose_final = loose[-1]
    rigid_final = angle_loose[-1]*rot
    rigid_residual = np.linalg.norm(rigid_final-loose_final)/np.linalg.norm(loose_final)
    kernel_error = max(abs(loose_bond[:-1]-kernel[idx]))
    mean_error = np.max(abs(gripped@fixed.T-drive@fixed.T))
    if max(work_residual, kernel_error, mean_error, max(abs(angle_gripped))) > 1e-8:
        raise ArithmeticError('Original kernel, constrained displacements or virtual work mismatch')
    def fraction_removed(index):
        # Compare absolute losses from each protocol's own instantaneous load.
        loose_loss = loose_bond[0]-loose_bond[index]
        gripped_loss = grip_bond[0]-grip_bond[index]
        return (loose_loss-gripped_loss)/loose_loss
    emit('stationary_angle_control', law=law, chi=level, trace_path=str(path),
         trace_sha256=digest(path), finite_times=times[:-1].tolist(), last_row_is_infinite_time=True,
         archived_rerest=False, layer_mean_y_separation=xy[COUNT:, 1].mean()-xy[:COUNT, 1].mean(),
         free_angle_per_q=angle_loose.tolist(), fixed_angle_max_per_q=max(abs(angle_gripped)),
         relative_norm_nonrigid_part_of_free_static_motion=rigid_residual,
         fixed_means_max_error=mean_error, old_kernel_max_abs_error=kernel_error,
         force_virtual_work_max_abs_error=work_residual,
         material_free=loose_bond.tolist(), material_gripped=grip_bond.tolist(),
         layer_actuator_conjugate_gripped=conjugate.tolist(),
         angle_actuator_torque_per_q=torque.tolist(), affine_angle_per_q=affine_angle,
         fixed_angle_velocity_drag_coefficient=grip_drive@grip_drive,
         gripped_material_T_retention=grip_bond[-2]/grip_bond[0],
         gripped_material_infinite_retention=grip_bond[-1]/grip_bond[0],
         fraction_T_absolute_material_loss_removed=fraction_removed(-2),
         fraction_infinite_absolute_material_loss_removed=fraction_removed(-1))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trace-dir', type=Path, required=True)
    args = parser.parse_args()
    for law in LAWS:
        for level in LEVELS:
            read_control(args.trace_dir/f'{law}_{level:.2f}.npz', law, level)
