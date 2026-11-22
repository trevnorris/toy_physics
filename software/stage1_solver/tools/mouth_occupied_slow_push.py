"""Slow entrance push on an occupied draining mouth with two flexible bends.

The far feed and bulk dipole data are zero. The applied load is conjugate to
the common translation, with no force applied to either relative bend.
"""
import argparse
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

from mouth_carrier_lateral_response import prepare
from mouth_free_bend_drive import Response
from mouth_free_bend_release import pencil
from mouth_flared_receiver import emit
from mouth_occupied_companion import OccupiedMode, parts

TENSIONS = (.02, .04)
FREQUENCIES = (.005, .0025, .00125, .000625, .0003125, .00015625)
LOAD = np.array([1., 0., 0.])
PUSH = 1.e-12


def static_expansion(wall, wave):
    solve = wall.dipole.factor.solve
    tangent = solve(wall.bend)
    velocity = wall.motion + wall.skew @ tangent
    mass = (2/3) * velocity.T @ solve(velocity) + wave.action * wave.mass
    cross = wall.motion.T @ tangent
    gyro = (2/3) * (cross - cross.T - tangent.T @ (wall.skew @ tangent))
    np.testing.assert_allclose(gyro, -gyro.T, atol=1e-10)
    return mass, gyro


def exact_coefficients(stiffness, mass, gyro):
    bend = stiffness[1:, 1:]
    shape = np.r_[1., -np.linalg.solve(bend, stiffness[1:, 0])]
    load = (gyro @ shape)[1:]
    induced = load @ np.linalg.solve(bend, load)
    spring = stiffness[0] @ shape
    inertia = shape @ mass @ shape + induced
    if min(np.linalg.eigvalsh(stiffness)) <= 0 or inertia <= 0:
        raise ArithmeticError('Stationary state lacks positive slow support')
    return spring, inertia, induced, shape


def signed_push_check(wall, rate, extra, predicted):
    energy, time = pencil(wall, 0.)
    system = energy - rate*time + sparse.block_diag(
        (sparse.csc_matrix(wall.operator.shape), sparse.csc_matrix(extra)), format='csc')
    amplitudes = PUSH * np.array([1., -1., .5, -.5])
    rhs = np.zeros((system.shape[0], len(amplitudes)), complex)
    rhs[-3] = amplitudes
    solved = splu(system.tocsc()).solve(rhs)
    per_force = solved[-3:] / amplitudes
    error = np.max(abs(per_force - predicted[:, None])) / max(abs(predicted))
    residual = np.max(abs(system @ solved - rhs) / abs(amplitudes))
    if error > 5e-6 or residual > 5e-6:
        raise ArithmeticError(f'Full signed push failed: {error}, {residual}')
    return dict(max_response_disagreement=error, residual_per_force=residual,
                amplitudes=amplitudes.tolist(), max_station_excursion=np.max(abs(np.cumsum(solved[-3:],axis=0))))


def run(args):
    if args.bath:
        from mouth_occupied_push_baths import load_bath
        dom, field = load_bath(args.bath)
    else:
        dom, field = prepare(args)
    wave = OccupiedMode(args.step, args.pad)
    if wave.action != 1:
        raise ArithmeticError('This comparison requires wave action one')
    wall = Response(dom, field, 'wide')
    static = wall.force_matrix(0.) + wave.action*wave.static
    mass, gyro = static_expansion(wall, wave)
    coefficients = {}
    emit(test='setup', far=dom.far, roof=dom.roof, wave_action=wave.action,
         wall_area=wall.area.diagonal().tolist(), rigid_inertia=mass[0, 0],
         rigid_spring=static[0, 0], force_amplitude=PUSH)
    for tension in TENSIONS:
        spring, inertia, induced, shape = exact_coefficients(static+tension*wall.area, mass, gyro)
        coefficients[tension] = (spring, inertia)
        emit(test='slow', tension=tension, bath_spring=spring, entrance_inertia=inertia,
             inverse_inertia=1/inertia, bend_velocity_inertia=induced,
             static_stations=np.cumsum(shape).tolist(), pin_frequency=np.sqrt(spring/inertia))
    estimates = {t: [] for t in TENSIONS}
    for omega in FREQUENCIES:
        rate = 1j*omega
        carrier = wall.force_matrix(rate)
        wave_force = wave.action*wave.force(rate)
        for tension in TENSIONS:
            extra = wave_force + tension*wall.area
            matrix = carrier + extra
            np.testing.assert_allclose(matrix, matrix.T.conj(), atol=3e-9)
            q = np.linalg.solve(matrix, LOAD)
            stations = np.cumsum(q)
            z = 1 / q[0]
            spring, inertia = coefficients[tension]
            measured = (spring-z.real)/omega**2
            estimates[tension].append(measured)
            emit(test='push', omega=omega, tension=tension,
                 acceleration_per_force=parts(-omega**2*stations),
                 station_over_entrance=parts(stations/q[0]), rigid_acceleration_per_force=parts(-omega**2/matrix[0,0]),
                 dynamic_inertia=measured, spring_over_inertial_force=spring/(inertia*omega**2),
                 coordinate_residual=np.max(abs(matrix@q-LOAD)))
            if omega == .0025:
                emit(test='coupled_check', omega=omega, tension=tension, **signed_push_check(wall,rate,extra,q))
    for tension, values in estimates.items():
        estimate = (4*values[-1]-values[-2])/3
        error = abs(estimate/coefficients[tension][1]-1)
        emit(test='derivative_check', tension=tension, richardson=estimate, relative_error=error)
        if error > 1.e-5:
            raise ArithmeticError('Low frequency coefficient failed direct frequency comparison')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sources = parser.add_mutually_exclusive_group(required=True)
    sources.add_argument('--source', type=Path)
    sources.add_argument('--bath', type=Path)
    parser.add_argument('--step', type=float, default=.03)
    parser.add_argument('--pad', type=float, default=5.)
    run(parser.parse_args())
