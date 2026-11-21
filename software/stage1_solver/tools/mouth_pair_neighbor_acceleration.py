"""Slow prescribed acceleration of the renewed packet population in resting brane.

Each occupied age has its own birth time. A drive proportional to time since
birth minus present age has zero present velocity and unit acceleration.
Two harmonic experiments check that low-frequency derivative independently.
"""
import argparse
from hashlib import file_digest
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from mouth_pair_neighbor_impulse import (Carrier, PRIOR, constitutive_data, emit,
                                         feeding_packet, grid, samples, transit)


FREQUENCIES = np.array([.08, .04])
DRIVE_RIM = np.tile([1., 0.], 4)


def age_response(path, law):
    clocks, mass, scale = path
    spring, drag, contact = law
    coeff = spring[:, None]
    # Rows: velocity, three extensions, rim contact extension, two force integrals,
    # and the occupied population's momentum. Columns: rim/support for each drive.
    state = np.zeros((8, 8), dtype=complex)

    def drivers(t):
        waves = np.exp(1j * FREQUENCIES * t)
        return np.repeat([1., t, *waves], 2)

    def observe(signal, t, phase):
        result = signal * phase
        result[2:4] -= t * signal[:2]
        return result

    for start, stop, attached in zip(clocks[:-1], clocks[1:], (True, False, True)):
        if not attached:
            state[4] = 0.

        def rates(t, z):
            z = z.reshape(state.shape)
            stiffness, vol = scale(t)
            drive = drivers(t)
            phase = np.ones(8, dtype=complex)
            phase[4:] = drive[4:].conj()
            visc = vol * drag[:, None]
            neighbor = stiffness * coeff * z[1:4]
            far = (z[0] - drive * (1.-DRIVE_RIM) + (neighbor/visc).sum(axis=0)) / (1./visc).sum()
            rim = stiffness * contact * z[4] if attached else np.zeros_like(far)
            out = np.empty_like(z)
            out[0] = -(far + rim) / mass
            out[1:4] = (far - neighbor) / visc
            out[4] = z[0] - drive * DRIVE_RIM if attached else 0.
            out[5] = observe(rim, t, phase) / mass
            out[6] = observe(far, t, phase) / mass
            out[7] = observe(z[0], t, phase)
            return out.ravel()

        sol = solve_ivp(rates, (start, stop), state.ravel(), rtol=1e-9,
                        atol=3e-12, max_step=clocks[-1]/192, method='DOP853')
        if not sol.success:
            raise ArithmeticError(sol.message)
        state = sol.y[:, -1].reshape(state.shape).copy()

    phase = np.ones(8, dtype=complex)
    phase[4:] = drivers(clocks[-1])[4:].conj()
    export = observe(state[0], clocks[-1], phase)
    budget = state[5] + state[6] + export
    budget[2:4] += state[7, :2]
    budget[4:] += 1j * np.repeat(FREQUENCIES, 2) * state[7, 4:]
    return state[5:7], float(abs(budget).max())


def opening_weights(opening):
    target, = [Path(p) for p in opening['input_pins'] if 'alone' in p or 'baseline' in p]
    geo = grid(*opening['feed_grid'])
    data, _ = feeding_packet(target, {target:opening['input_pins'][str(target)]}, geo, [0])
    radius, mu, volume, fields = samples(geo, data, 'cell')
    inward = volume * fields[:, 0]
    _, groups = np.unique(radius.round(12), return_inverse=True)
    grouped = np.bincount(groups, weights=inward)
    rank = ((grouped.cumsum() - .5*grouped) / grouped.sum())[groups]
    return rank, inward, inward*(1.-mu**2)


def run(folder):
    for name in ('child-outcome', 'outcome'):
        if json.loads((folder/f'{name}.json').read_text())['exitCode']:
            raise ValueError('Failed static reference')
    if not json.loads((folder/'limit-validation.json').read_text())['verified']:
        raise ValueError('Unverified static reference')
    rows = [json.loads(line) for line in (folder/'stdout').read_text().splitlines()]
    rays = [r for r in rows if r['test']=='ray' and r['policy']=='reconnect']
    opening, = [r for r in rows if r['test']=='opening' and r['policy']=='reconnect'
                and r['d']==8 and r['method']=='cell' and r['feed_grid'][:2]==[96,6]]
    state, mesh = opening['state'], opening['clock_mesh']
    archive = json.loads(PRIOR.read_text())['runs']['mouth_mixed_order_ripple02']
    records = map(json.loads, archive['stdout_lines'])
    ref, = [r for r in records if r['test']=='archive_patch'
            and (r['mesh'],r['profile'],r['state'])==(mesh,'cubic',state)]
    carrier = Carrier(Path(ref['carrier_file']).parent)
    with Path(ref['order_file']).open('rb') as stream:
        if file_digest(stream,'sha256').hexdigest()!=ref['order_sha256']:
            raise ValueError('Wrong order preparation')
    with np.load(ref['order_file']) as data:
        order = data['chi']
    law = constitutive_data()
    responses, ages = [], []
    for ray in rays:
        path = transit(carrier, order, ray['eta'], .501 if state=='near' else 1.15)
        np.testing.assert_allclose(path[0], ray['clocks'], atol=1e-10,rtol=0)
        gain, balance = age_response(path, law)
        c, b, _, minus_g, bv, _ = ray['responses']
        expected = [[minus_g,-(c+minus_g)], [bv,-(b+bv)]]
        dc_error = float(abs(gain[:, :2]-expected).max())
        emit(test='acceleration_ray', state=state, eta=ray['eta'],
             momentum_balance_error=balance, static_reference_error=dc_error)
        if balance>2e-8 or dc_error>3e-7:
            raise ArithmeticError('Age balance or independent static response disagrees')
        responses.append(gain)
        ages.append(path[0][-1])

    rank, inward, tangent = opening_weights(opening)
    labels = [r['rate_label'] for r in rays]
    table = np.asarray(responses)
    # Accumulate per-label weights once, including endpoint extension used in the static run.
    interp = np.column_stack([np.interp(rank, labels, col) for col in np.eye(len(labels))])
    gains = np.tensordot(tangent @ interp, table, axes=1)
    dc, derivative, slow, slower = np.split(gains, 4, axis=1)
    dc, derivative = dc.real, derivative.real
    for mode, impedance in [('fixed',None), ('free',0.), ('matched',-dc[1,1])]:
        def effective(matrix):
            if impedance is None:
                return matrix[0,0]
            return matrix[0,0] + matrix[0,1]*matrix[1,0]/(impedance-matrix[1,1])
        w0 = 0. if impedance is None else dc[1,0]/(impedance-dc[1,1])
        w1 = 0. if impedance is None else (derivative[1,0]+derivative[1,1]*w0)/(impedance-dc[1,1])
        material = derivative[0,0] + derivative[0,1]*w0 + dc[0,1]*w1
        harmonic = [effective(matrix) for matrix in (slow,slower)]
        slopes = np.imag(harmonic)/FREQUENCIES
        errors = abs(slopes-material)
        emit(test='opening_acceleration',state=state,mesh=mesh,rays=len(rays),mode=mode,
             impedance=impedance,dc=dc.tolist(),derivative=derivative.tolist(),
             steady_drag=-effective(dc),signed_actuator_mass=-material,
             omega=FREQUENCIES.tolist(),harmonic_force_over_velocity=[[a.real,a.imag] for a in harmonic],
             harmonic_actuator_mass=(-slopes).tolist(),absolute_errors=errors.tolist(),
             active_mass=inward@interp@ages,projected_active_mass=tangent@interp@ages,
             support_velocity_per_V=w0,support_velocity_per_A=w1)
        if errors[1] > max(2e-8, .32*errors[0]) or errors[1]>max(2e-8,.005*abs(material)):
            raise ArithmeticError('Finite-frequency extrapolation disagrees with acceleration response')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference', type=Path)
    run(parser.parse_args().reference)
