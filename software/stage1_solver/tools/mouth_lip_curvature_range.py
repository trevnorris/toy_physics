"""Check the local surface expansion along representative meridional and hoop curvatures."""
import json

import numpy as np

from mouth_lip_order_moment import patch


def main():
    sigma = 1./3
    bend = -5./6 * sigma
    second = (np.pi**2/3-2)*sigma
    for hoop in (.02, .04, .07, .10):
        for meridian_ratio in (-1., 0., 1.):
            k = np.array([meridian_ratio*hoop, hoop, hoop])
            H = sum(k)
            s2 = k[0]*k[1]+k[0]*k[2]+k[1]*k[2]
            quadratic_energy = sigma + bend*H**2/2 + second*s2
            quadratic_moments = bend*H + second*(H-k)
            rows = []
            for fraction in (.83,.93):
                extent = min(24. if fraction>.9 else 20., fraction/hoop)
                result = patch('front', k, 5, extent)
                rows.append(dict(extent=extent, energy=result['energy'],
                                 moments=result['moments']))
            print(json.dumps(dict(curvatures=k.tolist(), cases=rows,
                                  quadratic_energy=quadratic_energy,
                                  fractional_energy_error=(quadratic_energy-rows[-1]['energy'])/sigma,
                                  quadratic_moments=quadratic_moments.tolist(),
                                  maximum_absolute_moment_error=max(abs(quadratic_moments-rows[-1]['moments'])),
                                  maximum_moment_extent_difference=max(abs(np.array(rows[0]['moments'])-rows[1]['moments'])))), flush=True)


if __name__ == '__main__':
    main()
