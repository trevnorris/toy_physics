"""Release the fixed shape in the same steady incident stream, at linear order.

Incremental coordinates are measured from the clamped, stream-loaded state.
The carrier and wave evolve after the opposing holder is smoothly removed.
"""
import argparse
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

from mouth_carrier_lateral_response import prepare
from mouth_flared_receiver import emit
from mouth_free_bend_drive import Response
from mouth_free_bend_held_flow import volume_reaction, wall_reaction
from mouth_occupied_companion import OccupiedMode

DEEP_TENSIONS = (.02, .04)
FIRST_TENSION = .02
RELEASE_TIME = 4.
END = 64.
AVERAGING_RATE = .35
SNAPSHOTS = (1., 2., 4., 8., 16., 32., 48., 64.)


def released_fraction(t):
    return (1 - np.cos(np.pi * min(t / RELEASE_TIME, 1.))) / 2


def check_diagnostics(state, field, wave_position, wave_velocity, wall, wave,
                      local_stiffness, work):
    # Rotating-frame quadratic energy, including the carrier/wall cross term.
    field_energy = (np.sum(field * (wall.operator @ field), axis=0)
                    - 2 * np.sum((wall.bend.T @ field) * state, axis=0)) / 3
    wall_energy = .5 * np.einsum('ik,kij,jk->k', state, local_stiffness, state)
    wave_energy = (np.sum(abs(wave_velocity)**2, axis=0)
                   + np.sum((wave_position.conj() * (wave.dipole @ wave_position)).real, axis=0)
                   - 2 * np.sum((wave.q.T @ wave_position.real) * state, axis=0)
                   ) * wave.action / (2 * wave.omega)
    return field_energy + wall_energy + wave_energy - work


def run(args):
    if args.bath:
        from mouth_occupied_push_baths import load_bath
        dom, background = load_bath(args.bath)
    else:
        dom, background = prepare(args)
    wall = Response(dom, background, 'wide')
    incident, residual = wall.dipole.steady(1., 0.)
    surface = wall_reaction(dom, background, incident)
    load = volume_reaction(dom, background, incident, 'wide')
    wave = OccupiedMode(args.wave_step, args.pad)
    support = np.array([wall.area * [0., FIRST_TENSION, deep] for deep in DEEP_TENSIONS])
    static = wall.force_matrix(0.) + wave.action * wave.static + support
    emit(test='held', generalized_force_per_speed=load.tolist(), direct_pressure=surface.tolist(),
         field_residual=residual, first_tension=FIRST_TENSION, deep_tensions=DEEP_TENSIONS,
         wave_action=wave.action, min_static_eigenvalues=np.linalg.eigvalsh(static)[:,0].tolist(),
         dt=args.dt, far_stream=args.speed, release_duration=RELEASE_TIME,
         feed_radius=dom.far, receiver_roof=dom.roof)
    if max(abs(load-surface)) > .0003 or np.min(np.linalg.eigvalsh(static)) <= 0:
        raise ArithmeticError('Held pressure or bend stability failed')
    # Condense each implicit-midpoint step at its center. Keeping the imaginary
    # amplitude evolves both sidebands of the real occupied mode.
    dt = args.dt
    s = 2 / dt
    carrier_solve = splu((wall.operator - s*wall.skew).tocsc()).solve
    carrier_shape = carrier_solve(wall.bend + s*wall.motion)
    carrier_force = (2/3) * (wall.bend - s*wall.motion).T
    wave_shift = s*s - 2j*wave.omega*s
    wave_solve = splu((wave.dipole + wave_shift*sparse.eye(wave.dipole.shape[0],format='csc')).tocsc()).solve
    wave_shape = wave_solve(wave.q.astype(complex))
    wave_scale = wave.action / wave.omega
    local = (wall.frozen + support + wave.action*wave.static + wave_scale*(wave.q.T @ wave.first))
    matrices = local - carrier_force @ carrier_shape - wave_scale*(wave.q.T @ wave_shape.real)
    # Columns here describe distinct supports under the identical physical load.
    state = np.zeros((3, len(DEEP_TENSIONS)))
    field = np.zeros((wall.operator.shape[0], len(DEEP_TENSIONS)))
    pos = np.zeros((wave.dipole.shape[0], len(DEEP_TENSIONS)), complex)
    vel = np.zeros_like(pos)
    work = np.zeros(len(DEEP_TENSIONS))
    absolute_work = work.copy()
    weighted_stations = state.copy()
    test_rate = AVERAGING_RATE
    for i in range(round(END/dt)):
        t = (i+1) * dt
        released = released_fraction(t-dt/2)
        force = args.speed*released*load[:,None]
        base_field = carrier_solve(-s*(wall.skew @ field + wall.motion @ state))
        base_wave = wave_solve(s*vel + wave_shift*pos)
        rhs = (force + (2/3)*s*(wall.motion.T @ field) + carrier_force @ base_field
               + wave_scale*(wave.q.T @ base_wave.real))
        center = np.linalg.solve(matrices, rhs.T[...,None])[...,0].T
        mid_wave = base_wave + wave_shape @ center
        increment = 2*(center-state)
        delta_work = np.sum(force*increment,axis=0)
        work += delta_work
        absolute_work += abs(delta_work)
        vel = 2*s*(mid_wave-pos)-vel
        pos = 2*mid_wave-pos
        field = 2*(base_field+carrier_shape@center)-field
        state += increment
        stations = np.cumsum(state, axis=0)
        weight = np.exp(-test_rate*t) * dt * (.5 if t == END else 1)
        weighted_stations += weight*stations
        if t in SNAPSHOTS:
            error = check_diagnostics(state, field, pos, vel, wall, wave, local, work)
            relative_error = abs(error) / np.maximum(absolute_work, 1.e-40)
            emit(test='release', time=t, removed=released_fraction(t),
                 stations_per_speed=(stations/args.speed).T.tolist(),
                 work=work.tolist(), energy_minus_work=error.tolist(),
                 relative_energy_error=relative_error.tolist())
            if not np.isfinite(stations).all() or max(relative_error) > 1.e-6:
                raise ArithmeticError('Incremental release work failed')
    # Independent continuous-frequency operator. Compare the finite-duration
    # quadrature to the infinite-time transform (the small omitted tail and
    # discrete time errors are both included in the discrepancy).
    rate = test_rate
    drive_frequency = np.pi/RELEASE_TIME
    transform = drive_frequency**2*(1+np.exp(-rate*RELEASE_TIME))/(2*rate*(rate**2+drive_frequency**2))
    impedance = wall.force_matrix(rate) + wave.action*wave.force(rate) + support
    transformed = np.linalg.solve(impedance, np.broadcast_to(load,(len(support),3))[...,None])[...,0]
    predicted = transform*np.cumsum(transformed,axis=1)
    measured = weighted_stations.T / args.speed
    relative = abs(measured-predicted)/np.maximum(abs(predicted),1.e-15)
    emit(test='laplace', rate=rate, window=END, predicted=predicted.tolist(),
         integrated=measured.tolist(), discrepancy=relative.tolist())
    if max(relative.ravel()) > .025:
        raise ArithmeticError('Release trace and independent causal transform disagree')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sources = parser.add_mutually_exclusive_group(required=True)
    sources.add_argument('--source', type=Path)
    sources.add_argument('--bath', type=Path)
    parser.add_argument('--wave-step', type=float, required=True)
    parser.add_argument('--pad', type=float, required=True)
    parser.add_argument('--dt', type=float, choices=(.25,.125,.0625), required=True)
    parser.add_argument('--speed', type=float, default=1.e-6)
    args = parser.parse_args()
    if args.speed == 0:
        parser.error('A nonzero signed far-brane speed is required for normalization')
    run(args)
