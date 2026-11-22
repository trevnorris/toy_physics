"""Can angular diffusion of bound sites change the sign of their order force?

Test a conditional material law on the two fixed spherical shells. A site
diffuses in its relaxed local adsorption energy and is replaced isotropically
with mean lifetime tau. The finite-volume calculation uses signed angular
potentials, conserves mean coverage, and resolves the steady distribution.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import solve_banded

from mouth_curved_site_front import sha
from mouth_curved_site_pair import material_calibration, validated
from neutral_patch_sites import TEMPERATURE, emit


def hemisphere_response(cells, speed, zeta, energy_dipole, temperature):
    # Uniform bins in x=cos(theta) have equal area on S². Potential is the
    # linear-in-flow site energy, divided by kBT. The rates on each face
    # have the exact Boltzmann ratio for the two adjacent bin centers.
    edges = np.linspace(-1.,1.,cells+1)
    x = (edges[:-1]+edges[1:])/2
    dx = 2./cells
    potential = speed*energy_dipole*x/temperature
    strength = zeta/2 * (1-edges[1:-1]**2)/dx**2
    forward = strength*np.exp(-np.diff(potential)/2)
    reverse = strength*np.exp(np.diff(potential)/2)
    band = np.zeros((3,cells))
    band[1] = 1.
    band[1,:-1] += forward
    band[1,1:] += reverse
    band[0,1:] = -reverse
    band[2,:-1] = -forward
    n = solve_banded((1,1),band,np.ones(cells))
    np.testing.assert_allclose(n.mean(),1.,rtol=0.,atol=3e-11)
    if min(n)<=0:
        raise ArithmeticError('Negative angular coverage')
    return n,3*np.mean(x*n)


def run(folder):
    rows = [r for r in validated(folder) if r['test']=='binding']
    if len(rows)!=4:
        raise ValueError('Expected four site-law and location combinations')
    _,scale = material_calibration(5)
    thermal = TEMPERATURE*scale
    emit(test='thermal_example', patch_kBT=TEMPERATURE, specific_energy_scale=scale,
         kBT_per_site_mass=thermal, input=str(folder/'stdout'), input_sha256=sha(folder/'stdout'))
    for row in rows:
        keys = {k:row[k] for k in ('mesh','location','width','law')}
        energy = row['energy_slope']*row['order_dipole_overlap_per_U']
        frozen = row['fixed_count_medium_force_per_mass_U']
        coefficient = row['count_dipole_force_coefficient']
        beta_equilibrium = -energy/thermal
        mobile = frozen+coefficient*beta_equilibrium
        fraction = -frozen/(coefficient*beta_equilibrium)
        critical = fraction/(1-fraction) if 0<fraction<1 else None
        emit(test='coverage_limit', **keys, polar_energy_per_mass_U=energy,
             fixed_count_force=frozen, full_equilibrium_force=mobile,
             equilibrium_fractional_count_dipole_per_U=beta_equilibrium,
             required_equilibrium_fraction=fraction, critical_2Ds_tau_over_radius2=critical)
        # These two mobility ratios straddle the transition in the sampled
        # material. The exact thresholds above are calculated without a scan.
        for zeta in (.01,.1):
            beta = beta_equilibrium*zeta/(1+zeta)
            predicted = frozen+coefficient*beta
            for ncell in (64,256):
                for speed in (.002,.001):
                    n1,p = hemisphere_response(ncell,speed,zeta,energy,thermal)
                    n2,m = hemisphere_response(ncell,-speed,zeta,energy,thermal)
                    measured = (p-m)/(2*speed)
                    err = np.max(abs(n1-n2[::-1]))
                    emit(test='signed_transport', **keys, cells=ncell, drive=speed,
                         zeta=zeta, pole_reflection_error=err,
                         predicted_count_dipole=beta, measured_count_dipole=measured,
                         predicted_force=predicted, finite_volume_force=frozen+coefficient*measured)
                    if err>1e-11 or abs(measured-beta)>4e-4*abs(beta):
                        raise ArithmeticError('Signed angular mobility check')


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder',type=Path)
    run(parser.parse_args().folder)
