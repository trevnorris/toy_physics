"""Read the occupied pole and outgoing wave power at the two archived exit states."""
import argparse
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from mouth_exit_release import Release
from mouth_exit_witness import inputs
from mouth_live_waist_active_wave import ActiveWave
from mouth_live_waist_wave import OccupiedWave
from mouth_sealed_shape_relax import Meridian
from mouth_shallow_puncture_fields import PORT
from mouth_strained_wave import install


def read_mode(trial, shifts, chi, seed, cavity_fraction):
    shape = Meridian(trial.radius, trial.x+shifts[:4])
    coefficients = trial.co.copy()
    coefficients[0] += shifts[4]
    wave = OccupiedWave(shape, trial.n, 2*trial.n, trial.mesh, PORT, seed)
    install(wave, shape, coefficients, core_reference=trial.core_reference)
    wave = ActiveWave.install(SimpleNamespace(wave=wave)).wave
    weights = wave.project @ chi.T.ravel()
    cavity_mass = wave.matrix(weights[:, None, None]*wave.mc)
    mass = wave.cut_m+cavity_mass
    stiffness = wave.cut_k+wave.matrix(weights[:, None, None]*wave.kc)
    omega, v = seed
    port = wave.port
    edge = v[port.ids]
    mv, kv = mass@v, stiffness@v
    norm = np.vdot(v, mv).real
    fraction = np.vdot(v, cavity_mass@v).real/norm
    error = (kv-omega**2*mv)
    error[port.ids] -= port.block(omega)@edge
    relative_error = np.linalg.norm(error)/np.linalg.norm(kv)
    # Cycle-averaged port work is computed as a single boundary quadratic form.
    # The stored energy includes the cut feeding slab as well as the cavity.
    flux = .5*omega.real*np.vdot(edge, port.block(omega.real)@edge).imag
    stored = .25*(abs(omega)**2*norm+np.vdot(v, kv).real)
    rate = -2*omega.imag
    result = dict(radial_exit=shape.exit, cavity_kinetic_fraction=fraction,
                  pole=[omega.real, omega.imag], pole_energy_decay_rate=rate,
                  real_frequency_port_power=flux, stored_energy_to_port=stored,
                  power_over_stored=flux/stored, energy_efold_periods=omega.real/(2*np.pi*rate),
                  pole_residual=relative_error, cavity_fraction_error=abs(cavity_fraction-fraction))
    print(json.dumps(result, allow_nan=False), flush=True)
    if relative_error>1e-7 or abs(fraction-cavity_fraction)>1e-8 or min(rate, stored, flux)<=0:
        raise ArithmeticError('Stored fields or the outgoing work are unresolved')


def run(radius, mesh, folder):
    problem = Release(radius, mesh)
    start = problem.start()
    path = folder/f'r{radius}n{mesh}'/'endpoint.npz'
    shifts, chi, seed, digest, status, _, _ = inputs(path, radius, mesh, problem.trial)
    if status != 'balanced':
        raise ValueError('Invalid endpoint origin')
    raw = (path.parent/'stdout').read_bytes()
    final, = [p for p in map(json.loads, raw.splitlines()) if p['kind']=='endpoint']
    print(json.dumps(dict(radius=radius, mesh=mesh, input_sha256=problem.trial.fingerprints,
                          endpoint_sha256=digest, parent_stdout_sha256=hashlib.sha256(raw).hexdigest(),
                          row_order=['prior_three_coordinate','released_four_coordinate'])), flush=True)
    read_mode(problem.trial, start.shifts, start.chi, start.seed, start.info['cavity_fraction'])
    read_mode(problem.trial, shifts, chi, seed, final['cavity_fraction'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--origin-folder', type=Path, required=True)
    args = parser.parse_args()
    for radius in (1,2):
        for mesh in (48,72):
            run(radius, mesh, args.origin_folder)
