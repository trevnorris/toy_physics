"""Compare unit-energy releases with two independent localized-torque experiments."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

from mouth_order_wave_pulse import digest, emit
from mouth_shear_ringdown_compare import read as read_burst


def harmonic_prediction(t, data, train, future):
    """Fit one constant-frequency sinusoid per port on train and freeze all coefficients."""
    y = data/np.sqrt(np.mean(data[train]**2, axis=0))
    def evaluate(omega, full=False):
        angle = omega*(t-t[train][0])
        basis = np.column_stack((np.cos(angle), np.sin(angle)))
        coef = np.linalg.lstsq(basis[train], y[train], rcond=None)[0]
        diff = basis@coef-y
        if full:
            return dict(frequency=omega, train_error=(np.linalg.norm(diff[train],axis=0)/np.linalg.norm(y[train],axis=0)).tolist(),
                        prediction_error=(np.linalg.norm(diff[future],axis=0)/np.linalg.norm(y[future],axis=0)).tolist())
        return np.sum(diff[train]**2)
    grid = np.arange(1., 16.01, .03)
    best = grid[np.argmin([evaluate(omega) for omega in grid])]
    fit = minimize_scalar(evaluate, bounds=(best-.03, best+.03), method='bounded', options={'xatol':1e-12})
    return evaluate(fit.x, True)


def main(path):
    load = lambda f: json.loads((path/f).read_text())
    inv = load('invocation.json')
    outcome = load('outcome.json')
    limits = load('limit-validation.json')
    if not (limits['verified'] and outcome['limitsVerified'] and outcome['exitCode'] == 0
            and outcome['childOutcome']['exitCode'] == 0 and inv['swapMax'] == 0 and not inv['seconds']):
        raise ValueError(f'Guard failed for {path}')
    s, end = [json.loads(x) for x in (path/'stdout').read_text().splitlines()]
    if s['test'] == 'burst_preparation':
        s, d = read_burst(path)
        t = d['time']-s['drive_duration']
        work = d['work'][-1]
        data = np.column_stack((d['near_velocity'], d['far_velocity']))
        record = dict(active=(d['ordered']+d['upper']+d['lower'])/work, feed=(d['feed']+d['loss'])/work,
                      cut=d['cut_export']/work, weak=d['tail']/work)
        err = end['max_balance_over_final_work']
    else:
        src = path/('trace.npz' if s['test'] == 'preparation' else 'probe.npz')
        if digest(src) != end['trace_sha256']:
            raise ValueError('Time history digest differs')
        with np.load(src) as f:
            d = dict(zip(f['columns'], f['trace'].T, strict=True))
        if s['test'] == 'preparation':
            t = d['time']
            data = None
            record = dict(active=d['active_throat'], feed=d['feed_net'], cut=d['cut_export'], weak=d['weak_throat'])
            ledger = sum(record.values())+d['receiver']-1.
            err = np.max(abs(ledger))
            np.testing.assert_allclose(ledger, d['closure'], atol=2e-14, rtol=0)
        else:
            t = d['time']-s['drive_duration']
            data = np.column_stack((d['projection_velocity'], d['patch_velocity']))
            work = d['work'][-1]
            record = dict(feed=d['feed_net']/work, cut=d['cut_loss']/work)
            err = np.max(abs(d['closure']))/work
        np.testing.assert_allclose(err,end['max_energy_error'],atol=2e-14)
    if err > 1e-8:
        raise ArithmeticError('Energy ledger failed')
    def between(lo, hi): return (t>=lo)&(t<hi)
    windows = {}
    for lo, hi in ((4,8),(14,18),(44,48),(94,98),(144,148),(194,198)):
        if t[-1] < hi-.12: continue
        keep = between(lo,hi)
        row = {k:v[keep].mean() for k,v in record.items()}
        if s['test'] == 'local_source':
            power = np.mean(d['projection_velocity'][keep]**2)
            row['kinetic_fraction_in_prepared_shape'] = power/d['twice_kinetic'][keep].mean()
            row['kinetic_fraction_of_waist_interior_in_prepared_shape'] = power/d['behind_twice_kinetic'][keep].mean()
            row['twice_mean_projected_kinetic_per_input'] = power/work
        windows[f'{lo}-{hi}'] = row
    slopes = {}
    if '194-198' in windows:
        for key in ('feed','cut'):
            slopes[key] = (windows['194-198'][key]-windows['44-48'][key])/150
    fits = {}
    if data is not None and t[-1] >= 48:
        fits['train14-30_predict32-48'] = harmonic_prediction(t,data,between(14,30),between(32,48))
        if t[-1] >= 198:
            fits['train60-100_predict170-198'] = harmonic_prediction(t,data,between(60,100),between(170,198))
    emit(run=path.name, energy_error=err, windows=windows, net_power_secants_per_initial=slopes,
         frozen_two_port_predictions=fits)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('paths', type=Path, nargs='+')
    for path in p.parse_args().paths:
        main(path)
