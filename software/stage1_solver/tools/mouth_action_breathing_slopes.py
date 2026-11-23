"""Exit-frequency slopes of three coupled axisymmetric lattice trial modes."""
import argparse
import hashlib
import json

import numpy as np

from mouth_action_breathing import coordinates
from mouth_exit_witness import inputs
from mouth_poly_shallow_release import BIAS, CUT, ROOT, Order, Trial, emit, exterior, source
from mouth_live_waist_fields import Carrier
from mouth_poly_bonded_core import core_source
from mouth_sealed_shape_relax import Meridian
from mouth_trapping_modes import Sector, pole


def projected_energy(trial, displaced, order, archived_wave):
    """Native static terms and the occupied pole using its full-domain eigenvector."""
    shape = Meridian(trial.radius, trial.x + displaced[:4])
    inner = trial.co.copy()
    outer = trial.cr.copy()
    inner[0] += displaced[4]
    outer[0] += displaced[5]
    depths = np.linspace(0, 1, 2*trial.n+1)
    sites = np.linspace(0, 1, trial.n+1)
    inner_elastic = source(shape, inner, depths, sites, False)
    inner_elastic += core_source(shape, inner, *trial.core_reference, depths, sites)
    receiver = exterior(shape, trial.n, CUT, side_order_work=lambda w, r:
                        source(shape, outer, w/shape.length, r, True))
    order_energy = Order(shape, trial.n).energy_gradient(order)[0]
    carrier = Carrier(shape, trial.n, 2*trial.n, trial.mesh)
    carrier_energy = carrier.response(order, BIAS)[0]
    operator = Sector(trial, displaced, order, 4)
    full = np.zeros(len(operator.live), complex)
    full[operator.allpoints[:, 0] > 0] = archived_wave[1]
    wave, data = pole(operator, (archived_wave[0], full[operator.live]))
    if data['seed_overlap'] < .95 or data['residual'] > 1e-7:
        raise ArithmeticError('Unresolved outgoing background wave')
    return (order_energy + receiver['energy'] + receiver['static_work']
            + np.sum(inner_elastic*order) + carrier_energy + wave[0].real)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('radius', type=int, choices=(1, 2))
    parser.add_argument('--complete-polar', action='store_true',
                        help='Adopt a partial signed readout and use the full-domain pole for missing directions.')
    args = parser.parse_args()
    radius = args.radius
    mesh = 72
    folder = ROOT / '_scratch/s11c/mouth_action_breathing01/run'
    guard = json.loads((folder / 'outcome.json').read_text())
    if guard['exitCode'] or guard['childOutcome']['exitCode'] or not guard['limitsVerified']:
        raise ValueError('Original guarded load and lattice inertia unavailable')
    raw = (folder / 'stdout').read_bytes()
    rows = [json.loads(row) for row in raw.splitlines()]
    rows = [r for r in rows if r['radius'] == radius and r['n'] == mesh]
    setup, = [r for r in rows if r['kind'] == 'setup']
    states = {r['relative_exit']: r for r in rows if r['kind'] == 'basis_state'}
    modes = [r for r in rows if r['kind'] == 'three_shape_trial']
    path = ROOT / f'_scratch/s11c/mouth_exit_release01/r{radius}n{mesh}/endpoint.npz'
    trial = Trial(radius, mesh, bonded_core=True)
    shifts, order, seed, digest, status, *_ = inputs(path, radius, mesh, trial)
    if digest != setup['endpoint_sha256'] or status != 'balanced':
        raise ValueError('Preparation mismatch')
    original = json.loads((path.parent/'stdout').read_text().splitlines()[0])['exit_radius']
    b = setup['physical_exit']
    emit(kind='input', radius=radius, n=mesh, parent_stdout_sha256=hashlib.sha256(raw).hexdigest(),
         endpoint_sha256=digest)
    completed = set()
    if args.complete_polar:
        old = ROOT / f'_scratch/s11c/mouth_action_breathing01/slopes{radius}'
        guard = json.loads((old/'outcome.json').read_text())
        if not guard['limitsVerified'] or guard['exitCode'] != 1:
            raise ValueError('Expected the explicitly recorded incomplete directional probe')
        old_raw = (old/'stdout').read_bytes()
        previous = list(map(json.loads, old_raw.splitlines()))
        if previous[0]['endpoint_sha256'] != digest:
            raise ValueError('Different background in retained readings')
        completed = {(r['index'], r['relative_exit_step']) for r in previous[1:]}
        center_error = abs(projected_energy(trial, shifts, order, seed) - states[0.]['energy'])
        if center_error > 1e-8:
            raise ValueError('Full-domain pole failed to reproduce archived center energy')
        emit(kind='retained', prior_sha256=hashlib.sha256(old_raw).hexdigest(),
             samples=sorted(completed), center_energy_error=center_error)
    # A real mode's first eigenfrequency derivative needs only its unperturbed
    # eigenvector in the generalized Rayleigh quotient (stationarity principle).
    for mode in modes:
        v = np.array(mode['basis'])
        v /= max(abs(v))
        for exit_step, mode_step in ((.012, .002), (.006, .001)):
            if (mode['index'], exit_step) in completed:
                continue
            frequencies = []
            for offset in (-exit_step, 0., exit_step):
                center, shape = coordinates(trial, shifts, b*(1 + offset))
                state = states[offset]
                direction_energy = []
                for sign in (-1, 1):
                    dq = sign*mode_step*v
                    displaced = center.copy()
                    displaced[0] += dq[0]
                    displaced[1] += np.log1p(original*dq[1]/shape.exit)
                    displaced[2] += dq[2]
                    if args.complete_polar:
                        energy = projected_energy(trial, displaced, order, seed)
                    else:
                        energy, _, _, _, info, _ = trial.evaluate(displaced, order, seed)
                        if info['mode_overlap'] < .95 or info['pole_residual'] > 2e-7:
                            raise ValueError('Degree-four preparation lost in signed energy change')
                    direction_energy.append(energy-state['energy'])
                stiffness = sum(direction_energy)/mode_step**2
                mass = v @ np.array(state['mass']) @ v
                if min(stiffness, mass) <= 0:
                    raise ValueError('Unstable trial variation')
                frequencies.append(float(np.sqrt(stiffness/mass)))
            mismatch = abs(frequencies[1]/mode['freq'] - 1)
            if mismatch > .002:
                raise ValueError('Independent directional curvature did not recover the archived mode')
            emit(kind='three_shape_derivative', radius=radius, n=mesh, index=mode['index'],
                 relative_exit_step=exit_step, shape_step=mode_step,
                 rayleigh_frequencies=frequencies, center_relative_error=mismatch,
                 domega_dexit=(frequencies[2]-frequencies[0])/(2*b*exit_step))


if __name__ == '__main__':
    main()
