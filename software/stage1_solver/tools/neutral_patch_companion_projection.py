"""Angular average of a tangential bath velocity projected on the companion axis.

At fixed time each angular patch has the same causal linear scalar response in
its local meridian, and patches have uniform density per solid angle. Under
these assumptions projection and angular averaging commute with that response.
The resulting drive gives the mean force per active patch at first order.
"""
import argparse
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import Legendre, leggauss

from mouth_order_wave_pulse import digest, emit


def main(source, output):
    mu, weights = leggauss(24)
    with np.load(source,allow_pickle=False) as field:
        degrees = field['degrees']
        # The supplied fields contain equatorial velocities. One sine projects
        # the phase gradient to the meridian, another the rim load to the axis.
        coefficients = []
        for degree in degrees:
            derivative = Legendre.basis(int(degree)).deriv()
            coefficients.append(float(weights@((1-mu**2)*derivative(mu))/2/derivative(0)))
        result = {key:field[key] for key in ('w_nodes','wall_edges')}
        for d in (8,16):
            modes = field[f'modes{d}']
            np.testing.assert_allclose(modes.sum(axis=0),field[f'velocity{d}'],rtol=1e-14,atol=1e-17)
            result[f'velocity{d}'] = np.einsum('l,lwp->wp',coefficients,modes)
    emit(test='angular_response',source=str(source),sha256=digest(source),
         degrees=degrees.tolist(),weights=coefficients)
    # Integration by parts gives integral(mu P_l(mu)) after the boundary
    # terms vanish. Orthogonality leaves only l=1.
    np.testing.assert_allclose(coefficients,np.where(degrees==1,2/3,0),atol=2e-14,rtol=0)
    np.savez_compressed(output,**result)
    emit(test='projection_output',file=str(output),sha256=digest(output))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    main(args.source,args.output)
