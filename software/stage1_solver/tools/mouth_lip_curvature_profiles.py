"""Exact normal-order profiles on saved conditional sealed-duct geometries."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize_scalar

from mouth_lip_order_moment import patch
from mouth_sealed_shape_relax import DELTA, Meridian


def normal_work(shape, t, chart_fraction):
    radius, slope, second = shape.profile(t)
    arc = np.hypot(1, slope)
    hoop = 1/(radius*arc)
    k = DELTA*np.array([-second/arc**3, hoop, hoop])
    sigma = 1./3
    bend = -5./6*sigma
    m2 = (np.pi**2/3-2)*sigma
    trace = sum(k)
    s2 = 2*k[0]*k[1]+k[1]**2
    approx = sigma + bend*trace**2/2 + m2*s2
    mom = bend*trace + m2*(trace-k)
    extent = min(24., chart_fraction/max(abs(k)))
    exact = patch('front', k, 5, extent)
    return dict(fractional_depth=t, radius=radius, curvatures=k.tolist(), extent=extent,
                quadratic=approx, normal_energy=exact['energy'],
                quadratic_moments=mom.tolist(), normal_moments=exact['moments'],
                relative_surface_energy_difference=(exact['energy']-approx)/sigma,
                area_per_t=4*np.pi*radius**2*arc*shape.length)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--radius', choices=(1,2), type=int, required=True)
    parser.add_argument('--trajectory', type=Path, required=True)
    args = parser.parse_args()
    outcome = json.loads((args.trajectory/'outcome.json').read_text())
    if outcome['exitCode'] or outcome['childOutcome']['exitCode'] or not outcome['limitsVerified']:
        raise ValueError('Training trajectory incomplete')
    raw = (args.trajectory/'stdout').read_bytes()
    rows = [json.loads(line) for line in raw.splitlines()]
    if rows[0]['radius'] != args.radius or rows[-1].get('reason') != 'thin_front':
        raise ValueError('Input is not a curvature-limited trajectory of this radius')
    states = [r for r in rows if r['test']=='state']
    for state in states[-2:]:
        shape = Meridian(args.radius, state['coordinates'])
        at_neck = minimize_scalar(lambda t: shape.profile(t)[0], bounds=(0.,1.),
                                  method='bounded', options={'xatol':1e-11})
        t_neck = min((0.,at_neck.x,1.), key=lambda t: shape.profile(t)[0])
        print(json.dumps(dict(radius=args.radius, state=state['index'], test='narrowest',
                              source=hashlib.sha256(raw).hexdigest(),
                              readings=[normal_work(shape,t_neck,f) for f in (.83,.93)])),flush=True)
        for N in (12,20):
            x,w = leggauss(N)
            readings=[normal_work(shape,t,.93) for t in (1+x)/2]
            wall_area = .5*sum(v*r['area_per_t'] for r,v in zip(readings,w))
            surface_error = .5*sum(v*r['area_per_t']*r['relative_surface_energy_difference']
                                   for r,v in zip(readings,w))
            print(json.dumps(dict(radius=args.radius, state=state['index'], test='integral', quadrature=N,
                                  area=wall_area, E_normal_minus_quadratic_sigma1=surface_error,
                                  total_with_normal=state['total']+surface_error,
                                  maximum_abs_fitted_vs_normal_density=max(abs(r['relative_surface_energy_difference'])
                                                                           for r in readings))),flush=True)
        if abs(wall_area-state['area']) > 1.e-5:
            raise ArithmeticError('Meridian changed relative to the saved trajectory')


if __name__ == '__main__':
    main()
