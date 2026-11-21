"""Stretch one material radius with a balanced load and all other coordinates free."""
import argparse
from pathlib import Path

import numpy as np
from scipy import linalg
from scipy.integrate import solve_ivp

from mouth_mobile_neighbor_probe import COUNT, TIMES, curvatures, incremental_gradient, relative_error, rigid_frame
from mouth_order_wave_pulse import digest, emit


def run(coupon):
    with np.load(coupon) as data:
        xy = data['upper_equilibrium']
        labels = data['labels'][:, 1]
    flat = xy.ravel()
    hessian = curvatures(flat)[0]
    change = incremental_gradient(flat)
    rigid = rigid_frame(flat)
    for label in labels:
        tag = label - COUNT
        radius = xy[tag] - xy.mean(axis=0)
        direction = radius / np.linalg.norm(radius)
        drive = (np.eye(COUNT)[tag, :, None] - 1 / COUNT) * direction
        g = drive.ravel()
        d = g / (g @ g)
        free = linalg.null_space(g[None, :])
        np.testing.assert_allclose(rigid.T @ g, 0., atol=1e-12)
        eig, v = linalg.eigh(free.T @ hessian @ free)
        modes = eig > 1e-8
        if np.count_nonzero(~modes) != 3 or eig.min() < -1e-8:
            raise ArithmeticError('Unresolved rigid modes under relative hold')
        cross = v[:, modes].T @ free.T @ hessian @ d
        initial = d @ hessian @ d
        gain = -np.expm1(-TIMES[:, None] * eig[modes]) / eig[modes]
        reading = initial - gain @ (cross**2)
        static = 1 / (g @ linalg.pinvh(hessian, atol=1e-8) @ g)
        np.testing.assert_allclose(reading[-1], static, rtol=1e-9)
        emit(test='internal_hold', tag=int(label), initial=initial, early=reading[2],
             transit=reading[-1], fraction=static/initial)
        if label != 12:
            continue
        def direct(amplitude):
            displaced = amplitude * d
            def force(q):
                return change(displaced + free @ q)
            sol = solve_ivp(lambda t, q: -free.T @ force(q), (0., TIMES[-1]),
                            np.zeros(free.shape[1]), method='BDF', atol=2e-13, rtol=2e-9,
                            jac=lambda t,q: -free.T @ curvatures(flat + displaced + free @ q)[0] @ free,
                            t_eval=TIMES)
            if not sol.success:
                raise RuntimeError(sol.message)
            return np.array([d @ force(q) for q in sol.y.T])
        for amp in (1e-5, 5e-6):
            actual = (direct(amp) - direct(-amp)) / (2*amp)
            err = relative_error(reading, actual)
            if err > 1e-6:
                raise ArithmeticError('Nonlinear relative hold outside tolerance')
            emit(test='direct_relative_hold', tag=int(label), amplitude=amp, peak_normalized_error=err)
    emit(test='coupon_used', sha256=digest(coupon))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--coupon', type=Path, required=True)
    run(p.parse_args().coupon)
