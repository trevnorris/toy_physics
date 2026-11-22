"""Independent depth, waist and order displacements at fixed slow core attachment."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import expit, logit

from mouth_poly_shallow_release import Trial, emit
from mouth_sealed_shape_relax import Meridian


def run(radius, mesh):
    folder = Path(f'_scratch/s11c/poly_bonded01/r{radius}n{mesh}')
    guard = json.loads((folder/'outcome.json').read_bytes())
    if not guard['limitsVerified'] or guard['exitCode'] or guard['childOutcome']['exitCode']:
        raise ValueError('The preceding one-direction relaxation is incomplete')
    rows = [json.loads(line) for line in (folder/'stdout').read_text().splitlines()]
    initial, direction, end = rows[0], rows[1], rows[-1]
    if end['kind'] != 'directional_minimum':
        raise ValueError('Independent release requires a resolved directional minimum')
    trial = Trial(radius,mesh,bonded_core=True)
    if initial['inputs'] != trial.fingerprints:
        raise ValueError('Mismatched frozen reference or seed')
    zero = np.zeros(6)
    e0,g0,order,fields,_,_ = trial.evaluate(zero,trial.chi)
    active = order.free & (trial.chi > 0.) & (trial.chi < 1.)
    chi0 = trial.chi[active]
    rate = -g0[active]/(2*order.stiff.diagonal()+1.5*order.lump).reshape(g0.shape)[active]
    rate /= np.max(abs(rate*chi0*(1-chi0)))
    logits = logit(chi0)
    frozen = end['step']*np.array(direction['displacement'])
    prepared = np.array([frozen[0], frozen[2], end['step']])

    def measure(x, seed, extra=None):
        shifts = frozen.copy()
        shifts[[0,2]] = x[:2]
        if extra is not None:
            shifts += extra
        shape = Meridian(radius,trial.x+shifts[:4])
        if shape.wall()['minimum_radius'] <= 0.:
            raise ValueError('Positive-radius chart ended')
        chi = trial.chi.copy()
        chi[active] = expit(logits+x[2]*rate)
        energy,g,_,material,info,_ = trial.evaluate(shifts,chi,seed)
        if (abs(info['current']-info['outflow']) > 1e-6*abs(info['current'])
                or info['mode_overlap'] < .94 or info['pole_residual'] > 1e-7):
            raise ArithmeticError('Current or wave continuation unresolved')
        dorder = rate*chi[active]*(1-chi[active])
        return energy, g[active]@dorder, material.wave.seed, info

    def slopes(x, seed, step):
        ds = []
        for axis in (0,1):
            dx = np.zeros(3)
            dx[axis] = step
            plus = measure(x+dx,seed)[0]
            minus = measure(x-dx,seed)[0]
            ds.append((plus-minus)/(2*step))
        return ds

    # Continue the pole along the already measured family to its saved minimum.
    # This reconstructs its vector only; no relaxation is repeated.
    seed = fields.wave.seed
    for fraction in (.25,.5,.75,1.):
        x = prepared*fraction
        extra = (fraction-1.)*frozen
        extra[[0,2]] = 0.
        seed = measure(x,seed,extra=extra)[2]
    x = prepared.copy()
    energy,ds,seed,info = measure(x,seed)
    if abs(energy-e0-end['change']) > 2e-6:
        raise ArithmeticError('Did not recover the starting minimum')
    emit(kind='origin', radius=radius, mesh=mesh, previous_stdout_sha256=hashlib.sha256(
        (folder/'stdout').read_bytes()).hexdigest(),inputs=trial.fingerprints)
    # A positive inverse metric supplies a length scale only to the numerical
    # search. Each accepted independent motion must lower the measured work.
    scale = max(1.,abs(direction['slopes'][2]))
    metric = np.eye(3)/scale
    gradient = np.array([*slopes(x,seed,.0004),ds])
    iteration = 0
    while True:
        emit(kind='independent_release',iteration=iteration,coordinates=x.tolist(),
             gradient=gradient.tolist(),change=energy-e0,**info)
        if np.max(abs(gradient)) < .001:
            break
        step = -metric@gradient
        if step@gradient >= 0.:
            raise ArithmeticError('Lost positive descent metric')
        step /= max(1.,np.max(abs(step))/.15)
        amplitude = 1.
        reason = 'residual exceeds resolved work'
        while True:
            nxt = x+amplitude*step
            try:
                candidate,ds,candidate_seed,new_info = measure(nxt,seed)
                if candidate <= energy+1e-4*amplitude*(step@gradient):
                    break
            except (ArithmeticError,ValueError,RuntimeError) as err:
                reason = str(err)
            amplitude *= .5
            if np.max(abs(amplitude*step)) < 1e-6:
                emit(kind='unresolved_release', reason=reason, last_gradient=gradient.tolist(),
                     last_coordinates=x.tolist())
                return
        next_gradient = np.array([*slopes(nxt,candidate_seed,.0004),ds])
        dx = nxt-x
        dg = next_gradient-gradient
        curvature = dx@dg
        if curvature > 1e-10*np.linalg.norm(dx)*np.linalg.norm(dg):
            left = np.eye(3)-np.outer(dx,dg)/curvature
            metric = left@metric@left.T+np.outer(dx,dx)/curvature
        x,energy,seed,info,gradient = nxt,candidate,candidate_seed,new_info,next_gradient
        iteration += 1

    # A pair of independent signed differences, then a Hessian from slopes.
    # The reserved shape direction holds all three solved coordinates fixed.
    controls = []
    hessian = np.zeros((3,3))
    for i in range(3):
        delta = np.zeros(3)
        delta[i] = .001
        plus, _, plus_seed, _ = measure(x+delta,seed)
        minus, _, minus_seed, _ = measure(x-delta,seed)
        p2, order_plus, _, _ = measure(x+.5*delta,seed)
        m2, order_minus, _, _ = measure(x-.5*delta,seed)
        controls.append([(plus-minus)/.002,(p2-m2)/.001])
        pg = [*slopes(x+.5*delta,plus_seed,.0004),order_plus]
        mg = [*slopes(x-.5*delta,minus_seed,.0004),order_minus]
        hessian[:,i] = (np.array(pg)-mg)/.001
    remaining = []
    for h in (.001,.0005):
        extra = np.zeros(6)
        extra[1] = h
        remaining.append((measure(x,seed,extra)[0]-measure(x,seed,-extra)[0])/(2*h))
    emit(kind='balance', coordinates=x.tolist(),change=energy-e0,signed_slope_pairs=controls,
         projected_hessian=hessian.tolist(),
         projected_eigenvalues=np.linalg.eigvalsh(.5*(hessian+hessian.T)).tolist(),
         independent_log_exit_slopes=remaining,**info)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--radius',type=int,choices=(1,2),required=True)
    parser.add_argument('--mesh',type=int,choices=(48,72),required=True)
    args=parser.parse_args()
    run(args.radius,args.mesh)
