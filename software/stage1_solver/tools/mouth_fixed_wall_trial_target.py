"""Recover the first previously rejected fine trial wall from its saved descent rule."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_joined_front_audit import Front
from mouth_joint_polar_shape import Trial
from mouth_polar_exit_readout import accepted
from mouth_supercore_joined_front import emit, file_sha


def recover(state, snapshot):
    prep, record, x, d, u = accepted(state, snapshot)
    directory = snapshot.parent
    with (directory/'outcome.json').open() as file:
        outcome = json.load(file)
    if (record['iteration'] != 0 or prep['material_cells'] != [33,32] or
            not outcome['limitsVerified'] or outcome['exitCode'] != 1
            or outcome['childOutcome']['exitCode'] != 1
            or '0.0008128832206617774' not in (directory/'stderr').read_text()):
        raise ValueError('Expected the recorded fine transfer and its material stall')
    p = directory.parent/'src/mouth_joint_polar_shape.py'
    if file_sha(p) != 'bcdcd7639d4eb13eacafca6d20bef8d33f690dbff4259b6993d0af9fef449eb6':
        raise ValueError('Unknown prior trajectory rule')
    base = Front(state, tuple(prep['field_grid']), 1)
    def energy(y):
        return Trial(base,state,y,11,32,False).energy(d,u)
    e = record['energy']
    h = 1e-4
    pairs = np.array([[energy(x+sign*h*v) for sign in (-1.,1.)] for v in np.eye(6)])
    g = (pairs[:,1]-pairs[:,0])/(2*h)
    diag = (pairs.sum(axis=1)-2*e)/h**2
    force_error = float(max(abs(g+record['global_forces'])))
    if force_error>5e-7 or abs(energy(x)-e)>2e-9:
        raise ArithmeticError('Cannot reproduce the starting energy and signed shape force')
    dx = -g/np.maximum(abs(diag),1.)
    dx *= min(1., .05/max(abs(dx)))
    for half in range(22):
        step = 2.**-half*dx
        fixed = energy(x+step)
        if fixed>e and half<20:
            continue
        armijo = e+1e-4*(g@step)
        if fixed>=armijo-1e-8:
            raise ArithmeticError('First evaluated trial cannot be uniquely identified by descent')
        return dict(kind='reconstructed_fixed_wall', accepted_sha256=file_sha(snapshot),
                    source_sha256=file_sha(p), x=(x+step).tolist(), original_x=x.tolist(),
                    signed_curvatures=diag.tolist(), start_energy=e, fixed_energy=fixed,
                    old_armijo_ceiling=armijo, trial_halvings=half, gradient_disagreement=force_error,
                    prior_attempt_residual=.0008128832206617774)
    raise ArithmeticError('Could not reconstruct previous trial geometry')


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('state',type=Path)
    parser.add_argument('snapshot',type=Path)
    args=parser.parse_args()
    emit(**recover(args.state,args.snapshot))
