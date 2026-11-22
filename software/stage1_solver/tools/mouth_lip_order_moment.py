"""Neutral Route A wave-off curvature energy and moment of a pinned coherent order layer."""
import argparse
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_bvp
from scipy.special import expit

from mouth_curved_order_patch import well


def emit(**row):
    print(json.dumps(row, allow_nan=False), flush=True)


def patch(kind, curvatures, resolution, extent):
    """Length unit delta, energy density a_B. Hold central physical three-area fixed."""
    k = np.asarray(curvatures)
    width = 25.0
    if kind == 'front':
        edges = [-extent, 0., extent]
        pins = [1., .5, 0.]
    elif kind == 'strip':
        edges = [-width/2-extent, -width/2, width/2, width/2+extent]
        pins = [0., .5, .5, 0.]
    else:
        raise ValueError(kind)
    if max(abs(k))*max(abs(x) for x in edges) >= .95:
        raise ValueError('Normal coordinates approach a focal surface')
    energy, second_moment, flat_profile_energy = 0., 0., 0.
    moments = np.zeros(3)
    max_residual = 0.
    midorder = None
    gauss, weights = leggauss(8)
    for index, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        s = np.linspace(lo, hi, int(np.ceil((hi-lo)*resolution))+1)
        def guess(s):
            if kind == 'front':
                f = expit(-s)
                return f, -f*(1-f)
            l, r = expit(s+width/2), expit(-s+width/2)
            return l*r, l*r*(r-l)
        def equation(t, y):
            c, deriv = y
            mean = np.sum(k[:, None]/(1.+k[:, None]*t), axis=0)
            return np.stack((deriv, c*(1-c)*(1-2*c)-mean*deriv))
        def ends(left, right):
            return [left[0]-pins[index], right[0]-pins[index+1]]
        sol = solve_bvp(equation, ends, s, np.stack(guess(s)), tol=2e-9,
                        max_nodes=20000)
        if not sol.success:
            raise ArithmeticError(sol.message)
        if lo <= 0. <= hi:
            midorder = float(sol.sol(0.)[0])
        max_residual = max(max_residual, max(sol.rms_residuals))
        half = np.diff(sol.x)[:, None]/2
        t = (sol.x[:-1,None]+half*(1+gauss)).ravel()
        ds = (half*weights).ravel()
        c, deriv = sol.sol(t)
        if min(c) < -1.e-8 or max(c) > 1.+1.e-8:
            raise ArithmeticError('Excursion outside order wells')
        area_factors = 1.+k[:,None]*t
        jac = np.prod(area_factors, axis=0)
        bare = (well(c)+deriv**2)*ds
        energy += bare@jac
        moments += ((jac*t)[None,:]/area_factors)@bare
        second_moment += bare@(t*t)
        c0, d0 = guess(t)
        flat_profile_energy += ((well(c0)+d0**2)*ds)@jac
    return dict(energy=energy, moments=moments.tolist(), second_moment=second_moment,
                rigid_profile_energy=flat_profile_energy, max_residual=max_residual,
                midorder=midorder)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--resolution', type=int, required=True)
    ap.add_argument('--extent', type=float, required=True)
    args = ap.parse_args()
    for kind in ('front', 'strip'):
        def evaluate(k):
            return patch(kind, k, args.resolution, args.extent)
        zero = evaluate([0.,0.,0.])
        emit(kind=kind, curvature=[0.,0.,0.], **zero)
        # Train at |k delta|=0.016, reserve half, quarter and the mixed case.
        moments = {}
        for step in (.016, .008, .004):
            for direction in (-1., 1.):
                q = step*direction
                center = evaluate([q,0.,0.])
                de = 2.e-5
                work = (evaluate([q+de,0.,0.])['energy']-
                        evaluate([q-de,0.,0.])['energy'])/(2*de)
                emit(kind=kind, curvature=[q,0.,0.], energy_difference_moment=work, **center)
                if abs(center['moments'][0]-work)>5.e-9:
                    raise ArithmeticError('Signed virtual moment/energy disagreement')
                moments[q] = center['moments'][0]
        def secant(q):
            return (moments[q]-moments[-q])/(2*q)
        B = (4*secant(.004)-secant(.008))/3
        m2 = zero['second_moment']
        emit(kind=kind, result='uniaxial', B_over_a_delta3=B,
             B_over_one_sigma_delta2=3*B, quadratic_jacobian_over_one_sigma_delta2=3*m2,
             signed_training_slope=secant(.016),
             predicted_half=secant(.016)*.008, observed_half=moments[.008],
             flat_applied_moment=zero['moments'][0], scale_delta_over_Rref=.04)
        if kind == 'front':
            mixed = np.array([-.006,.009,.009])
            held = evaluate(mixed)
            predicted = B*sum(mixed)+m2*(sum(mixed)-mixed)
            emit(kind=kind, curvature=mixed.tolist(), predicted=predicted.tolist(), **held)
            if max(abs(predicted-held['moments'])) > 2.e-5:
                raise ArithmeticError('Mixed-curvature validation failed')


if __name__ == '__main__':
    main()
