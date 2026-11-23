"""Neutral, unit-wound aperture geometry and available subsonic current at fixed vortex density."""
import argparse
import json

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from mouth_circulation_material import profile


def apertures(exits):
    fields, residual = profile(1, 64., 1.e-9)

    def density(s):
        return float(fields(s)[0]**2)

    targets = (.1, .5, .9, .99)
    radii = {str(target): brentq(lambda s: density(s)-target, .001, 32.) for target in targets}
    sonic = brentq(lambda s: s*density(s)**2-1, .1, 20.)
    print(json.dumps(dict(kind='neutral',winding=1,bulk_density=1.,healing_definition='density fraction',
                          depletion_radii=radii,sonic=sonic,density_at_sonic=density(sonic),
                          collocation_residual=residual)),flush=True)
    for b in exits:
        if not np.isfinite(b) or not 0<b<=32:
            raise ValueError('Positive exit radius within the computed density profile required')
        # A three-ball slice at cylindrical distance s has measure 4*pi*s*sqrt(b²-s²) ds.
        # For this fixed neutral density, an axial test velocity must share the
        # sound-speed budget with the swirl, even if other meridional velocities vanish.
        annulus_start = max(sonic,radii['0.5'])
        area_fraction = max(0.,1-(annulus_start/b)**2)**1.5
        def cap(s):
            n=density(s)
            allowed = max(0., n**4 - 1/s**2)
            return 4*np.pi*s*np.sqrt(max(0., b*b-s*s))*n*np.sqrt(allowed)
        if b>annulus_start:
            current,err=quad(cap,annulus_start,b,epsabs=1e-8,epsrel=1e-9)
            measured,area_err=quad(lambda s:4*np.pi*s*np.sqrt(b*b-s*s),annulus_start,b)
        else:
            current,err,measured,area_err = (0.,)*4
        volume=4*np.pi*b**3/3
        result=dict(kind='neutral_cross_section',exit_radius=b,core_reference_density_at_rim=density(b),
                    azimuthal_mach_at_rim=1/(b*density(b)**2),
                    fraction_dense_and_subsonic_with_swirl_only=area_fraction,
                    integrated_fraction=measured/volume,volume_integration_error=area_err,
                    max_kinematic_dense_subsonic_intake_at_fixed_density=current,capacity_quadrature_error=err)
        print(json.dumps(result,allow_nan=False),flush=True)
        if abs(area_fraction-measured/volume)>1e-9 or err>1e-6*(1.+current):
            raise ArithmeticError('Disagreement in neutral aperture measure')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('exit_radii',type=float,nargs='+')
    apertures(p.parse_args().exit_radii)
