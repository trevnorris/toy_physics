"""Two-sided stationary carrier with a prescribed transition in its equation of state.

The feed coefficient is 1/4 at w<=0. A C-infinity step of width 1/8 connects
it to the bulk coefficient while preserving the homogeneous rest state n=1.
The potential term in the weak equation and its real Newton derivative both
scale by 4*c(w) relative to the constant-coefficient reference solver.
"""
import argparse
import json
from copy import copy

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_curved_phase import Domain, RATE, cap, integrated_radial_current, readings


WIDTH = 1/8
FEED_COEFFICIENT = 1/4


def transition(w):
    """Normalized compact smooth step and its physical derivative."""
    x = w/WIDTH
    step = np.asarray(x >= 1, dtype=float)
    gradient = np.zeros_like(step)
    active = (x > 0) & (x < 1)
    z = x[active]
    left = np.exp(-1/z)
    right = np.exp(-1/(1-z))
    step[active] = left/(left+right)
    gradient[active] = (left*right/(left+right)**2)*(1/z**2+1/(1-z)**2)/WIDTH
    return step, gradient


def quadrature_coordinates(domain):
    """Evaluate physical w in the element/quadrature order used for assembly."""
    roots, _ = leggauss(4)
    sub_p, sub_v = np.meshgrid((roots+1)/2, (roots+1)/2, indexing='ij')
    depths = []
    for block, (pg, vg, _) in domain.blocks.items():
        i, j = np.indices((len(pg)-1, len(vg)-1))
        i, j = i.ravel()[:, None], j.ravel()[:, None]
        p = pg[i] + np.diff(pg)[i]*sub_p.ravel()
        v = vg[j] + np.diff(vg)[j]*sub_v.ravel()
        _, w, *_ = domain.chart(block, p, v)
        depths.append(w.ravel())
    return np.concatenate(depths)


def density_energy(n):
    # (n**5-1)/5 - (n-1), factored to preserve precision near rest.
    return (n-1)**2*(n**3+2*n**2+3*n+4)/5


def run(n, far):
    domain = Domain(n, far)
    depths = quadrature_coordinates(domain)
    step, gradient = transition(depths)
    r = np.linspace(0, 1, 2001)
    for bulk_coefficient in (1/8, 1/2):
        c = FEED_COEFFICIENT + (bulk_coefficient-FEED_COEFFICIENT)*step
        solver = copy(domain)
        # Domain.solve uses its weight solely for the local nonlinear operator.
        # Keep the physical integration weight on the original domain.
        solver.weight = domain.weight*c/FEED_COEFFICIENT
        field, station = solver.solve(RATE)
        h, hp = cap(domain, field, r, 1)
        crossing = readings(domain, field, r, h, hp)
        radii = np.array([.99, 1.])
        heights, _ = cap(domain, field, radii, 1)
        cylinder = -4*np.pi*radii**2*integrated_radial_current(
            domain, field, radii, heights, 64)
        shape_curvature = (1+domain.slope**2)**1.5/domain.rounding
        # q'=j_r*(h_wall''-h_cut'') on a tangent smooth trace, in r coordinates.
        qprime = crossing['lip_jr']*(shape_curvature+domain.slope)
        # Carrier force from explicitly spatial energy c(w) G(n). Integrate
        # over the computational domain; this depends on the frozen material.
        density = abs(domain.N@field)**2
        force = -np.dot(domain.weight, density_energy(density)*gradient
                        *(bulk_coefficient-FEED_COEFFICIENT))
        row = dict(n=n, far=far, coefficient_feed=FEED_COEFFICIENT,
                   coefficient_bulk=bulk_coefficient, width=WIDTH, **station,
                   center_q=crossing['center_q'], q_099=crossing['q_at_099'],
                   min_to_099=crossing['min_to_099'], cap_rate=crossing['total'],
                   raw_tangent_lip_q=crossing['raw_lip_q'],
                   raw_wall_error=crossing['lip_wall_error'], lip_jr=crossing['lip_jr'],
                   tangent_continuum_q=0, tangent_qprime_from_slip=qprime,
                   cylinder_rate=cylinder[-1], shell_099_to_1=cylinder[1]-cylinder[0],
                   carrier_force_from_frozen_gradient=force)
        print(json.dumps(row), flush=True)
        if not (abs(crossing['total']-RATE) < 4e-4 and abs(cylinder[-1]-RATE) < 4e-4):
            raise ArithmeticError('Integrated flux failed the reservoir check')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, required=True)
    parser.add_argument('--far', type=float, required=True)
    args = parser.parse_args()
    run(args.n, args.far)
