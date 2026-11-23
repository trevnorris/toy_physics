"""Remeasure the exit work at progressively finer wave resolution on one fixed material field."""
import argparse
import copy
import hashlib
import json

import numpy as np

from mouth_gallery_action import archived
from mouth_poly_shallow_release import ROOT, emit
from mouth_trapping_angular_tail import recorded_seed
from mouth_trapping_modes import Sector, pole


def checked_rows(path):
    guard = json.loads((path / 'outcome.json').read_bytes())
    if not guard['limitsVerified'] or guard['exitCode'] or guard['childOutcome']['exitCode']:
        raise ValueError(f'Unqualified source {path}')
    raw = (path / 'stdout').read_bytes()
    return list(map(json.loads, raw.splitlines())), hashlib.sha256(raw).hexdigest()


def measure(radius):
    trial, shape_coords, chi, _, digest, _ = archived(radius, 72)
    scratch = ROOT / '_scratch/s11c'
    modes, modal_hash = checked_rows(scratch / f'mouth_trapping_tail_mesh01/r{radius}')
    loads, load_hash = checked_rows(scratch / f'mouth_gallery_action01/r{radius}n72')
    if modes[0]['endpoint_sha256'] != digest or loads[0]['digest'] != digest:
        raise ValueError('The two calculations use different material states')
    emit(kind='inputs', radius=radius, endpoint_sha256=digest, mode_hash=modal_hash, load_hash=load_hash)
    for mesh in (28, 36, 48):
        target, = [row for row in modes if row.get('wave_mesh') == mesh]
        wave_trial = copy.copy(trial)
        wave_trial.mesh = mesh
        center = Sector(wave_trial, shape_coords, chi, 12)
        seed, _ = recorded_seed(center, target)
        for step in (.00125, .000625):
            load, = [row for row in loads if row.get('log_step') == step]
            signed = []
            for sign in (-1, 1):
                coords = shape_coords.copy()
                coords[1] += sign * step
                shifted = Sector(wave_trial, coords, chi, 12)
                if not np.array_equal(shifted.live, center.live):
                    raise ValueError('The exit displacement changed wave support')
                root, read = pole(shifted, seed)
                if read['seed_overlap'] < .98 or read['residual'] > 1e-7:
                    raise ArithmeticError('Unresolved signed pole continuation')
                signed.append((shifted.shape.exit, root[0].real, read['seed_overlap']))
            left, right = signed
            force = -(right[1]-left[1])/(right[0]-left[0])
            if force <= 0:
                raise ArithmeticError('The refined wave no longer supplies outward force')
            action = load['physical_closing_load']/force
            emit(kind='balance', mesh=mesh, log_step=step, outward_force_per_action=force,
                 closing=load['physical_closing_load'], action=action, cycle_energy=action*seed[0].real,
                 signed_exit_frequency_overlap=signed)
            if mesh == 48:
                peaks, peak_hash = checked_rows(scratch / f'mouth_action_amplitude_budget01/r{radius}m48')
                native = peaks[0]
                if native['threshold_stdout_sha256'] != load_hash:
                    raise ValueError('Amplitude and load preparation mismatch')
                scale = np.sqrt(action/native['action'])
                fields = {name: {k: state['shear_bearing_cavity'][k]*scale for k in ('speed','strain')}
                          for name, state in native['sampled_peaks'].items()}
                floor = native['phase_and_angle_rms_lower_bounds']['shear_bearing_cavity']['strain']*scale
                emit(kind='amplitude', mesh=mesh, log_step=step, peak_hash=peak_hash, action=action,
                     sampled=fields, angular_independent_strain_floor=floor,
                     max_fraction_of_needed_force_if_shear_at_most_one=1/floor**2)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--radius', required=True, type=int, choices=(1, 2))
    measure(parser.parse_args().radius)
