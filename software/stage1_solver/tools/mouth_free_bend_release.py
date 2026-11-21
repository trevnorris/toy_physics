"""Unforced modes and short-push releases of the finite-bath carrier and free wall coordinates.

Use the full field-coordinate pencil, retaining the carrier as a dynamic variable.
All forces are conjugate to physical entrance or distal displacements.
"""
import argparse
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import LinearOperator, eigs, splu

from mouth_carrier_lateral_response import prepare
from mouth_flared_receiver import emit
from mouth_free_bend_drive import Response


def pencil(response, tension):
    factor = 2/3
    stiffness = response.frozen+tension*response.area
    b = sparse.csc_matrix(-factor*response.bend)
    p = sparse.csc_matrix(factor*response.motion)
    energy = sparse.bmat([[factor*response.operator, b], [b.T, stiffness]], format='csc')
    time = sparse.bmat([[factor*response.skew, p], [-p.T, sparse.csc_matrix((3,3))]],
                       format='csc')
    return energy, time


def modes(energy, time, tension, seed):
    shift = .08+.10j
    factor = splu((energy-shift*time).tocsc())
    inverse = LinearOperator(energy.shape, matvec=lambda x: factor.solve(time@x), dtype=complex)
    vals, vecs = eigs(inverse, k=8, which='LM', v0=seed, tol=1e-10)
    for index in np.argsort(abs(1/vals)):
        rate = shift+1/vals[index]
        v = vecs[:, index]
        av, bv = energy@v, time@v
        err = np.linalg.norm(av-rate*bv)/(np.linalg.norm(av)+abs(rate)*np.linalg.norm(bv))
        shape = v[-3:]/v[np.argmax(abs(v[-3:]))+len(v)-3]
        emit(test='mode', area_tension=tension, real=float(rate.real), imag=float(rate.imag),
             full_relative_residual=float(err), coordinates_real=shape.real.tolist(),
             coordinates_imag=shape.imag.tolist())
        if err > 3e-6:
            raise ArithmeticError('Mode failed the field and wall equations')


def release(energy, time, tension, dt):
    # (energy - time*d/dt)y = applied force. Implicit midpoint after release
    # preserves the quadratic energy of these linear, conservative equations.
    advance = splu((energy/2-time/dt).tocsc())
    past = -energy/2-time/dt
    ports = np.zeros((energy.shape[0], 2))
    ports[-3:, :] = [[1,1], [0,1], [0,1]]
    amplitude = 1.e-9
    state = np.zeros_like(ports)
    work = np.zeros(2)
    work_scale = np.zeros(2)
    error_max = np.zeros(2)
    for i in range(round(120/dt)):
        mid = (i+.5)*dt
        force = amplitude*np.sin(np.pi*mid/8)**2 if mid <= 8 else 0.
        future = advance.solve(past@state+force*ports)
        dw = force*np.sum(ports*(future-state), axis=0)
        work += dw
        work_scale += abs(dw)
        state = future
        t = (i+1)*dt
        if t in (8., 16., 32., 48., 64., 80., 96., 112., 120.):
            stored = .5*np.sum(state*(energy@state), axis=0)
            error = abs(stored-work)/np.maximum(work_scale,1.e-40)
            error_max = np.maximum(error_max,error)
            emit(test='release', area_tension=tension, t=t, dt=dt, peak_force=amplitude,
                 drives=['entrance','flare'], entrance=(state[-3]/amplitude).tolist(),
                 flare=(state[-3:].sum(axis=0)/amplitude).tolist(),
                 relative_excursion=((state[-2]+state[-1])/amplitude).tolist(),
                 energy=stored.tolist(), work=work.tolist(), work_balance_relative=error.tolist())
    if np.max(error_max) > 2e-5:
        raise ArithmeticError('Pulse work and stored field plus wall energy disagree')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--extension', choices=['compact','wide'], default='compact')
    p.add_argument('--dt', type=float, choices=[.25,.125], default=.25)
    p.add_argument('--modes-only', action='store_true')
    args = p.parse_args()
    dom, field = prepare(args)
    response = Response(dom, field, args.extension)
    seed = np.random.default_rng(7491).normal(size=response.operator.shape[0]+3)
    for tension in (0., .02):
        energy, time = pencil(response, tension)
        modes(energy, time, tension, seed)
        if not args.modes_only:
            release(energy, time, tension, args.dt)
