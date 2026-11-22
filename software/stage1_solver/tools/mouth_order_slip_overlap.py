"""Apply a common site mass and Gaussian footprint to the early work budget.

Use order headroom at entry alone to bound the site mass. Give the site every
advantage allowed by this necessary constraint when maximizing delivered work.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_order_slip import NU, input_states, spectrum
from neutral_patch_sites import emit


def work_rows(path):
    outcome = json.loads((path/'outcome.json').read_text())
    if outcome['exitCode'] or outcome['childOutcome']['exitCode'] or not outcome['limitsVerified']:
        raise ValueError('Unverified entry-work result')
    raw = (path/'stdout').read_bytes()
    return list(map(json.loads,raw.splitlines())), hashlib.sha256(raw).hexdigest()


def run(args):
    rows, sha = work_rows(args.work)
    segment, = [r for r in rows if r['test']=='segment']
    endpoint, domain = input_states(segment['mesh'])
    c0, c1 = segment['order']
    np.testing.assert_allclose(c0, endpoint[0]['order'],rtol=0,atol=1e-12)
    kx, k2, count = spectrum(args.length, args.modes)
    shape = np.exp(-.5*args.width**2*k2)
    norm = domain['mobility']/args.length**4
    def rate(c):
        return NU*(k2+(1-6*c*(1-c))/domain['delta']**2)
    # Fourier amplitudes all have the same sign at the centered, static probe.
    # While c decreases to .8 the local curvature decreases monotonically.
    peak0 = norm*(count @ (shape/rate(c0)))
    overlap0 = norm*(count @ (shape**2/rate(c0)))
    overlap1 = norm*(count @ (shape**2/rate(c1)))
    if not (peak0 > overlap0 > 0 and overlap1 > overlap0 and rate(c1).min()>0):
        raise ArithmeticError('Invalid susceptibility or local stability')
    max_overlap = (1-c0)*overlap1/peak0
    emit(test='common_footprint', mesh=segment['mesh'], host_side=segment['side'], width=args.width,
         length=args.length, modes=args.modes, work_stdout_sha256=sha, entry_peak_per_abs_source=peak0,
         entry_overlap_per_abs_source=overlap0, destination_overlap_per_abs_source=overlap1,
         entry_order_headroom=1-c0, max_static_overlap_on_segment=max_overlap)
    for r in rows:
        if r['test']!='local_equilibrium_first_excursion':
            continue
        s0, s1 = (abs(r['branch_samples'][i]['slope']) for i in (0,-1))
        m_max = (1-c0)/(s0*peak0)
        # B=m*s^2*C/2 is specific energy. Gamma <= |s| A/nu and A <= m |s0| C1.
        dressing_difference = (s1**2*overlap1-s0**2*overlap0)/2
        slip_work = s0**2*overlap1*segment['carrier_phase_rise']/NU
        credit = m_max*max(0.,dressing_difference+slip_work)
        shortfall = r['bare_specific_energy_rise']-segment['launch_kinetic']-credit
        emit(test='best_first_excursion_credit', law=r['law'], mass_upper_from_entry=m_max,
             necessary_overlap=r['necessary_max_site_weighted_static_shift'],
             max_overlap=max_overlap, max_site_drag_work=m_max*slip_work,
             dressing_change_at_max_mass=m_max*dressing_difference,
             best_combined_credit=credit, remaining_specific_energy_shortfall=shortfall)
        np.testing.assert_allclose(m_max*s0*overlap1,max_overlap,rtol=2.e-14)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work',type=Path,required=True)
    parser.add_argument('--length',type=float,default=.12)
    parser.add_argument('--width',type=float,default=.01)
    parser.add_argument('--modes',type=int,default=18)
    args=parser.parse_args()
    if min(args.length,args.width,args.modes)<=0:
        parser.error('Use positive dimensions and cutoff')
    run(args)
