"""Two clocks for dilute attachment response on the same held mouth.

The first test bounds height departures after an instantaneous order change
on a thermal clamped patch. The second bounds lateral redistribution for a
separate inertial, isotropic shell with an effective polar energy. Compare to
an independently integrated carrier transit without assigning a site lifetime.
"""
import argparse
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import roots_laguerre

from mouth_carrier_lateral_response import trace
from mouth_curved_attachment_exchange import Cell, Layer, LEVELS, Surface, archived_rows
from mouth_curved_site_front import front_read, prepare, sha
from mouth_curved_site_pair import validated
from mouth_flared_receiver import emit
from mouth_segment_contact_hold import PITCH
from neutral_patch_sites import TEMPERATURE


def transit(dom,psi,lo,hi,degree):
    edges=dom.blocks[2][0]
    edges=np.r_[lo,edges[(lo<edges)&(edges<hi)],hi]
    q,weights=leggauss(degree)
    widths=np.diff(edges)[:,None]/2
    w=(edges[:-1,None]+widths*(q+1)).ravel()
    _,_,rp,_,wp,_,f,fr,fw=trace(dom,2,w,np.zeros_like(w),psi[:,None])
    speed=((fr[:,0]*rp+fw[:,0]*wp)/f[:,0]).imag*wp/(rp*rp+wp*wp)
    if min(speed)<=0:
        raise ValueError('Outward wall characteristic in reference transit')
    return float((widths*weights).ravel()@(1/speed))


def check_turnover(cell,row,t,specific_temperature):
    mixed=row['law']=='retain_mean_attraction'
    lo,hi=row['z_lo'],row['z_hi']
    layer=Layer(cell,lo,hi,mixed)
    chi=row['chi']
    log0=layer.partition(chi)
    np.testing.assert_allclose(log0,row['log_K'],rtol=0.,atol=2e-10)
    for speed in (.002,.001,-.001,-.002):
        shift=row['order_dipole_per_U']*speed
        log1=layer.partition(chi+shift)
        logmirror=layer.partition(chi-shift)
        surface=np.logaddexp(cell.boundary(lo,chi+shift,mixed),cell.boundary(hi,chi+shift,mixed))
        # Maxwellian positive or negative velocity through each face. K is
        # an average over cell area, in units of one microscopic length.
        rate=np.sqrt(specific_temperature/(2*np.pi))/PITCH*np.exp(surface-log1)
        # Density-ratio second moment on the bound interval for initial
        # vs post-step conditional equilibrium. The potential is affine.
        log_second=logmirror+log1-2*log0
        if log_second < -3e-12:
            raise ArithmeticError('Density-ratio second moment below one')
        bound=min(1.,np.exp(.5*(np.log(t*rate)+log_second)))
        emit(test='finite_step_departure',law=row['law'],location=row['location'],
             z_lo=lo,z_hi=hi,drive=speed,post_step_chi=chi+shift,
             logK_after=log1,logK_reflected=logmirror,log_squared_density_ratio=log_second,
             post_step_equilibrium_gross_exit_rate=rate,bound_on_old_residents_leaving=bound,
             equilibrium_fractional_count_change=np.expm1(log1-log0))


def run(folder):
    local=validated(folder)
    geom,=[r for r in local if r['test']=='geometry']
    sample,=[r for r in archived_rows('mouth_bound_site_transfer','primary') if r['test']=='sample']
    previous=archived_rows('mouth_curved_companion_order','mobile48')
    surface=Surface(sample['side'],sample['shell'],sample['pitch'],True)
    thermal=TEMPERATURE*sample['scale']
    np.testing.assert_allclose(thermal,geom['specific_kBT'])
    cell=Cell(surface,LEVELS[-1])
    dom,psi,chi=prepare(48)
    ends={'entry':1.175,'half':float(front_read(dom,chi,[],1.)[0])}
    elapsed=transit(dom,psi,*ends.values(),8)
    coarse=transit(dom,psi,*ends.values(),4)
    emit(test='kinetic_inputs',mesh=48,thermal_specific=thermal,
         microscopic_length_in_mouth_units=PITCH,depths=ends,carrier_transit=elapsed,
         coarse_carrier_transit=coarse,exchange_sha256=sha(folder/'stdout'))
    np.testing.assert_allclose(elapsed,coarse,rtol=1e-10)
    for r in local:
        if r['test']=='partition' and r['resolution']=='fine' and r['z_lo']==.25:
            check_turnover(cell,r,elapsed,thermal)
    nodes,weights=roots_laguerre(32)
    for r in previous:
        if r['test']!='coverage_limit':
            continue
        radius=float(dom.radius(ends[r['location']])[0])
        # For an isotropic canonical stationary shell the step susceptibility
        # to an effective ε cos(theta) is proportional to <|n(t)-n(0)|²>/2.
        # Kinetic equipartition and time Cauchy give the upper bound below.
        bound=thermal*(elapsed/radius)**2
        free=weights@(2*np.sin(.5*np.sqrt(2*nodes*bound))**2)
        needed=r['required_equilibrium_fraction']
        min_age=radius*np.sqrt(needed/thermal)
        emit(test='lateral_clock',law=r['law'],location=r['location'],radius=radius,
             required_equilibrium_fraction=needed,single_transit_response_fraction_ceiling=bound,
             free_geodesic_step_response=free,necessary_rms_age=min_age,
             necessary_rms_age_in_transits=min_age/elapsed)
        if not 0<free<=bound:
            raise ArithmeticError('Free geodesic cross-check violated the finite-speed bound')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('exchange',type=Path)
    run(parser.parse_args().exchange)
