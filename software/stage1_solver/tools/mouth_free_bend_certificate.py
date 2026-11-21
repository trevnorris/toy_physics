"""Certify the finite-bath modes and check signed pulses over four pulse durations."""
import argparse
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

from mouth_carrier_lateral_response import prepare
from mouth_flared_receiver import emit
from mouth_free_bend_drive import Response
from mouth_free_bend_release import modes, pencil


def pulses(energy, time, tension, dt):
    # Two physically distinct push points. Use equal and opposite full pushes
    # at the entrance and opposite half pushes at both points as linearity controls.
    ports = np.zeros((energy.shape[0], 6))
    strengths = np.array([1., -1., .5, -.5, .5, -.5])*1e-9
    ports[-3, :] = strengths
    ports[-2:, 4:] = strengths[4:]
    # Solve for the state increment, independent of the earlier next-state update.
    factor = splu((energy*.5-time/dt).tocsc())
    state = np.zeros_like(ports)
    work = np.zeros(len(strengths))
    absolute_work = np.zeros_like(work)
    length = 8.
    for step in range(round(4*length/dt)):
        midpoint = (step+.5)*dt
        load = np.sin(np.pi*midpoint/length)**2 if midpoint < length else 0.
        change = factor.solve(load*ports-energy@state)
        state += change
        increment = load*np.sum(ports*change, axis=0)
        work += increment
        absolute_work += abs(increment)
        t = (step+1)*dt
        if t not in (length, 2*length, 4*length):
            continue
        stored = .5*np.sum(state*(energy@state), axis=0)
        error = abs(stored-work)/absolute_work
        gap = (state[-2]+state[-1])/strengths
        emit(test='short_pulse', area_tension=tension, t=t, dt=dt,
             strengths=strengths.tolist(), gap_per_peak=gap.tolist(),
             entrance_per_peak=(state[-3]/strengths).tolist(),
             flare_per_peak=(state[-3:].sum(axis=0)/strengths).tolist(),
             work=work.tolist(), stored=stored.tolist(), relative_error=error.tolist())
        if np.max(error) > 2e-5 or not np.all(np.isfinite(error)):
            raise ArithmeticError('Short-pulse work check failed')
        np.testing.assert_allclose(gap[:4], gap[0], rtol=1e-8, atol=1e-9)
        np.testing.assert_allclose(gap[4], gap[5], rtol=1e-8, atol=1e-9)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--dt', type=float, choices=[.25,.125], required=True)
    p.add_argument('--extension', choices=['compact','wide'], required=True)
    args = p.parse_args()
    dom, field = prepare(args)
    response = Response(dom, field, args.extension)
    sampled = dom.N@field
    # L has the gradient Gram matrix, this local multiplier, and a PSD density Hessian.
    field_floor = np.min((abs(sampled)**8-1)/4+1/response.dipole.rq**2)
    static = response.force_matrix(0.)
    seed = np.random.default_rng(60341).normal(size=response.operator.shape[0]+3)
    for tension in (0., .02):
        schur = static+tension*response.area
        values = np.linalg.eigvalsh((schur+schur.T)/2)
        emit(test='energy_curvature', area_tension=tension, field_quadrature_floor=field_floor,
             relaxed_curvatures=values.tolist())
        energy, time = pencil(response, tension)
        modes(energy, time, tension, seed)
        pulses(energy, time, tension, args.dt)
