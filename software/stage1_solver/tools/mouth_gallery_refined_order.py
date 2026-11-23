"""Quasistatic crystalline-order response in the measured coordinate at the refined balancing action."""
import argparse
import copy
from functools import cache

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import expit, logit

from mouth_exchange_carrier import BIAS
from mouth_exit_release import Release
from mouth_gallery_action import archived
from mouth_gallery_refinement import checked_rows
from mouth_poly_shallow_release import ROOT, emit
from mouth_trapping_angular_tail import recorded_seed
from mouth_trapping_modes import Sector, pole


def measure(radius):
    n = 72
    trial, shifts, chi, old_mode, digest, _ = archived(radius, n)
    scratch = ROOT / '_scratch/s11c'
    fine, action_hash = checked_rows(scratch / f'mouth_gallery_refinement01/r{radius}')
    original, order_hash = checked_rows(scratch / f'mouth_gallery_order01/r{radius}n{n}')
    spectrum, pole_hash = checked_rows(scratch / f'mouth_trapping_tail_mesh01/r{radius}')
    order_origin, = [r for r in original if r['kind'] == 'origin']
    if {digest, fine[0]['endpoint_sha256'], order_origin['order_field_sha'], spectrum[0]['endpoint_sha256']} != {digest}:
        raise ValueError('Mismatched material preparations')
    target, = [r for r in spectrum if r.get('wave_mesh') == 48]
    balance, = [r for r in fine if r['kind'] == 'balance' and r['mesh'] == 48 and r['log_step'] == .000625]
    action = balance['action']
    _, _, owner, material, _, elastic = trial.evaluate(shifts, chi, old_mode)
    wave_trial = copy.copy(trial)
    wave_trial.mesh = 48
    operator = Sector(wave_trial, shifts, chi, 12)
    seed, _ = recorded_seed(operator, target)
    family = Release(radius, n)
    coordinate = logit(chi[family.active])

    def field(q):
        c = chi.copy()
        c[family.active] = expit(coordinate + q*family.rate)
        return c

    base_order = owner.energy_gradient(chi)[0]
    base_carrier = material.carrier.response(chi, BIAS)[0]

    @cache
    def work(q):
        c = field(q)
        op = Sector(wave_trial, shifts, c, 12)
        if not np.array_equal(op.live, operator.live) or not np.array_equal(op.labels, operator.labels):
            raise ArithmeticError('Change of shear support along the order coordinate')
        mode, read = pole(op, seed)
        if read['seed_overlap'] < .98 or read['residual'] > 1e-7:
            raise ArithmeticError('Lost the prepared torsional wave')
        return (owner.energy_gradient(c)[0]-base_order + material.carrier.response(c, BIAS)[0]-base_carrier
                + np.sum(elastic*(c-chi)) + action*(mode[0].real-seed[0].real))

    minima = {r['label']: r for r in original if r['kind'] == 'order_minimum'}
    q0 = minima['threshold']['q']
    result = minimize_scalar(work, bracket=(q0, q0+.002), method='brent', options={'xtol': 1e-8})
    if not result.success or result.fun > 0:
        raise ArithmeticError('Unresolved stationary order')
    q = result.x
    sides = [[work(q-h)-result.fun, work(q+h)-result.fun] for h in (.0002, .0001)]
    if min(min(pair) for pair in sides) <= 0:
        raise ArithmeticError('Signed order changes did not resolve a minimum')
    off = field(minima['no_wave']['q'])
    on = field(q)
    diff = on-off
    ij = np.unravel_index(abs(diff).argmax(), diff.shape)
    emit(kind='refined_order', radius=radius, input_sha=digest, action_hash=action_hash, order_hash=order_hash,
         pole_hash=pole_hash, action=action, wave_mesh=48, material_mesh=n, q_off=minima['no_wave']['q'],
         q_on=q, signed_costs=sides, work_from_reference=result.fun, peak_signed=diff[ij], location=list(map(int,ij)),
         chi_off=off[ij], chi_on=on[ij])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--radius', type=int, choices=(1, 2), required=True)
    measure(parser.parse_args().radius)
