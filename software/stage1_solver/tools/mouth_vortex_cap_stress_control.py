"""Axial pressure on a three-ball cut by the neutral, straight coherent vortex.

The resting control distinguishes axial background traction and differential
flow-induced traction. The flat conversion sheet is in reversible coexistence.
"""
import argparse
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from mouth_circulation_material import profile, stress_check
from mouth_poly_shallow_release import emit


def main(radii):
    fields,_=profile(1,64.,1e-9)
    radial_balance=stress_check(fields,1)
    transition=brentq(lambda s:fields(s)[0]**2-.5,0.,5.)
    ones=brentq(lambda s:s*fields(s)[0]**4-1.,.001,5.)
    if radial_balance>1e-8:
        raise ValueError('Inconsistent circular vortex stress')
    for radius in radii:
        # w does not vary within the cut. Its ram and density-gradient normal
        # stresses vanish for the pure azimuthal flow. Isotropic pressure
        # relative to the undisturbed n=1 brane is (n^5-1)/5.
        def deficit(u):
            s=radius*np.sin(u)
            density=fields(s)[0]**2
            return 4*np.pi*radius**3*np.sin(u)*np.cos(u)**2*(1-density**5)/5
        end=np.pi/2
        angle=lambda s: np.arcsin(min(1.,s/radius))
        integrate=lambda lo,hi: quad(deficit,lo,hi,epsabs=1e-11,epsrel=1e-11)[0]
        total=integrate(0.,end)
        emit(kind='zero_drain_vortex_control', exit=radius, axial_pressure_deficit=total,
             dense_fraction_of_axial_deficit=integrate(angle(transition),end)/total,
             dense_total_subsonic_fraction_of_axial_deficit=integrate(angle(max(ones,transition)),end)/total,
             half_density_distance=transition,total_sonic_distance=ones,
             meridional_drain=0., radial_stress_residual=radial_balance)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('exit_radii',type=float,nargs='+')
    main(parser.parse_args().exit_radii)
