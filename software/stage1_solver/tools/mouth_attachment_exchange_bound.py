"""Bound signed changes in dilute height occupancy after a finite order step.

Each angle is a separate periodic clamped patch, reflecting at z=0. The gas
starts in equilibrium at one common activity throughout z>=0. Its new passive
dynamics preserve the canonical measure at the stepped potential. B=[0,2]
is the bound layer. Use an upper reflecting box and then increase its size:
the bounds use only a local partition and a stationary flux through z=2.
"""
import argparse
from pathlib import Path

import numpy as np

from mouth_curved_attachment_exchange import Cell, Layer, LEVELS, Surface, archived_rows, sheet_mean
from mouth_curved_site_front import sha
from mouth_curved_site_pair import validated
from neutral_patch_sites import BETWEEN, TEMPERATURE, emit, morse


def linear_limit(exchange,turnover):
    rows=validated(exchange)
    geom,=[r for r in rows if r['test']=='geometry']
    clock,=[r for r in validated(turnover) if r['test']=='kinetic_inputs']
    if clock['exchange_sha256']!=sha(exchange/'stdout'):
        raise ValueError('Mismatched thermal and transit samples')
    scale=geom['specific_scale']
    thermal=geom['kBT']
    pref=np.sqrt(clock['thermal_specific']/(2*np.pi))/clock['microscopic_length_in_mouth_units']
    elapsed=clock['carrier_transit']
    for row in rows:
        if not (row['test']=='partition' and row['resolution']=='fine' and
                row['z_lo']==0 and row['z_hi']==2):
            continue
        z=row['z_hi']
        mixed=row['law']=='retain_mean_attraction'
        tail=geom['image_count']*abs(morse(z,*BETWEEN)[0])+mixed*abs(sheet_mean(z,geom['pitch'])[0])
        # For A=dV/dchi in microscopic units, F''=-scale Var(A)/T.
        second=(row['free_energy_slope_specific']/scale)**2-thermal*row['free_energy_curvature_specific']/scale
        crossings=elapsed*pref*np.exp(row['log_boundary_upper']-row['log_K'])
        # Invariance of the stepped equilibrium measure differentiates the
        # count into int A(0)[1_B(t)-1_B(0)]/T in the ORIGINAL ensemble.
        # On initially bound residents, use Cauchy and the crossing count.
        # On initial outsiders A is bounded throughout the exterior.
        loss=np.sqrt(second*crossings)
        gain=tail*crossings
        bound=abs(row['order_dipole_per_U'])*(loss+gain)/thermal
        emit(test='infinitesimal_count_bound',location=row['location'],law=row['law'],
             exchange_sha256=sha(exchange/'stdout'),clock_sha256=sha(turnover/'stdout'),
             elapsed=elapsed,baseline_crossings_per_resident=crossings,
             bound_second_moment_of_energy_derivative=second,exterior_slope_bound=tail,
             signed_count_dipole_derivative_ceiling=bound,
             required_magnitude=abs(row['zero_force_count_dipole']))


def run(exchange,turnover):
    sample,=[r for r in archived_rows('mouth_bound_site_transfer','primary') if r['test']=='sample']
    data=validated(exchange)
    clock,=[r for r in validated(turnover) if r['test']=='kinetic_inputs']
    if clock['exchange_sha256']!=sha(exchange/'stdout'):
        raise ValueError('Clock and adsorption sample disagree')
    t=clock['carrier_transit']
    surface=Surface(sample['side'],sample['shell'],sample['pitch'],True)
    cell=Cell(surface,LEVELS[-1])
    # One-sided Maxwellian normal speed in mouth units, per microscopic length.
    pref=np.sqrt(clock['thermal_specific']/(2*np.pi))/clock['microscopic_length_in_mouth_units']
    emit(test='inputs',exchange_hash=sha(exchange/'stdout'),clock_hash=sha(turnover/'stdout'),
         carrier_transit=t,one_sided_speed_per_native_length=pref)
    for row in data:
        if not (row['test']=='partition' and row['resolution']=='fine' and
                row['z_lo']==0 and row['z_hi']==2):
            continue
        chi=row['chi']
        mixed=row['law']=='retain_mean_attraction'
        z=row['z_hi']
        layer=Layer(cell,0.,z,mixed)
        log0=row['log_K']
        np.testing.assert_allclose(layer.partition(chi),log0,atol=1e-10,rtol=0)
        # At z>=2 all pair distances exceed the Morse minimum, so the
        # magnitude of each attractive pair energy decreases with height.
        assert z>BETWEEN[2]
        max_slope_free=len(surface.origins)*abs(morse(z,*BETWEEN)[0])+mixed*abs(sheet_mean(z,surface.pitch)[0])
        for speed in (.002,.001):
            amplitude=abs(speed*row['order_dipole_per_U'])
            gains=[]
            losses=[]
            for shift in (-amplitude,amplitude):
                wall=cell.boundary(z,chi+shift,mixed)
                mirror=layer.partition(chi-shift)
                # Outside B the ratio of the original and stepped densities
                # at fixed activity is at most exp(|shift| max_slope_free/T).
                # Any new resident from outside crosses the upper plane.
                log_gain=np.log(t*pref)+wall-log0+amplitude*max_slope_free/TEMPERATURE
                # Cauchy with the stepped path measure and the old ensemble
                # restricted to B. The partition of the squared density
                # ratio involves 2*chi - (chi+shift) by affinity in chi.
                log_loss=.5*(np.log(t*pref)+wall+mirror-2*log0)
                gains.append(np.exp(log_gain))
                losses.append(np.exp(log_loss))
            # log W and log K are convex in order. The endpoint maxima bound
            # all angles; the sum in log_loss is also convex.
            gain=max(gains)
            loss=min(1.,max(losses))
            if loss>=1:
                raise ArithmeticError('No positive lower bound on total occupancy')
            # Integrate x*(n/n0-1) on the two hemispheres and divide by mean n.
            dipole=3*(gain+loss)/(4*speed*(1-loss))
            emit(test='finite_count_bound',location=row['location'],law=row['law'],drive=speed,
                 baseline_logK=log0,uniform_free_energy_slope_bound=max_slope_free,
                 gain_fraction_ceiling=gain,loss_fraction_ceiling=loss,
                 magnitude_of_signed_count_dipole_per_U_ceiling=dipole,
                 magnitude_for_force_reversal=abs(row['zero_force_count_dipole']))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('exchange',type=Path)
    parser.add_argument('turnover',type=Path)
    parser.add_argument('--linear-limit',action='store_true')
    args=parser.parse_args()
    (linear_limit if args.linear_limit else run)(args.exchange,args.turnover)
