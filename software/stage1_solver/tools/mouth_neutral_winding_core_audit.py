"""Measure the bare vortex by direct annular finite differences in its physical amplitude."""
import json

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.linalg import solve_banded
from scipy.optimize import brentq


def neutral_state(outer_radius, spacing):
    """Fix a unit-density bath at the remote boundary and enforce f(0)=0."""
    count = round(outer_radius/spacing)
    r = np.linspace(0., outer_radius, count+1)
    x = r[1:-1]
    h = r[1]
    left = (x-h/2)/(x*h*h)
    right = (x+h/2)/(x*h*h)
    diagonal = left + right + 1/x**2

    def residual(f):
        u=f[1:-1]
        return (left*(u-f[:-2]) + right*(u-f[2:]) + u/x**2
                + .5*(u**8-1)*u)

    f = r/np.sqrt(4+r*r)
    f[-1] = 1.
    for _ in range(40):
        res=residual(f)
        J = np.zeros((3, len(x)))
        J[0,1:] = -right[:-1]
        J[2,:-1] = -left[1:]
        J[1] = diagonal + .5*(9*f[1:-1]**8-1)
        update = solve_banded((1,1), J, -res)
        if max(abs(update)) < 1.e-11:
            break
        norm = res@res
        alpha=1.
        while alpha >= 1.e-7:
            candidate = f.copy()
            candidate[1:-1] += alpha*update
            candidate_residual = residual(candidate)
            if candidate.min()>=0 and candidate.max()<=1 and candidate_residual@candidate_residual < norm:
                f=candidate
                break
            alpha*=.5
        else:
            raise ArithmeticError('No physical descending Newton step')
    else:
        raise ArithmeticError('Unconverged amplitude iteration')
    if f.min() < 0 or np.any(np.diff(f)<=0):
        raise ArithmeticError('Unphysical stationary branch')
    density=PchipInterpolator(r, f*f)
    radii=[brentq(lambda s: density(s)-fraction, 0., outer_radius)
           for fraction in (.5,.90,.95,.99)]
    sonic=brentq(lambda s:s*density(s)**2-1,0.,outer_radius)
    result=dict(radius=outer_radius,h=h,half_density=radii[0],density_90=radii[1],
                density_95=radii[2],density_99=radii[3],sonic=sonic,
                density_at_sonic=float(density(sonic)),net_circulation=2*np.pi,
                max_residual=max(abs(residual(f))),max_update=max(abs(update)))
    print(json.dumps(result,allow_nan=False),flush=True)
    if result['max_residual']*h*h>3e-14 or radii[0]>=sonic:
        raise ArithmeticError('Bad neutral force residual or cutoff')
    return result


if __name__ == '__main__':
    medium = neutral_state(20.,1/160)
    levels = [neutral_state(40.,h) for h in (1/40,1/80,1/160,1/320)]
    fine=neutral_state(64.,1/320)
    for key in ('half_density','sonic'):
        prev,mid,end = [v[key] for v in levels[-3:]]
        order = np.log2(abs((prev-mid)/(mid-end)))
        estimate=end+(end-mid)/3
        variation=abs(fine[key]-end)
        print(json.dumps(dict(quantity=key,order=float(order),extrapolated=estimate,
                              finest_finite_box_change=variation,
                              intermediate_finite_box_change=abs(medium[key]-levels[-2][key]))),
              flush=True)
        if not 1.65<order<2.35 or variation>2e-9 or abs(estimate-end)>2e-6:
            raise ArithmeticError('Radial refinement or bath-distance control failed')
