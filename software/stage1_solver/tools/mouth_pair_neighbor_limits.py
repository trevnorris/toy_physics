"""Locate any change of sign in the slow acceleration force for a resistive support."""
import argparse
from pathlib import Path
import json

import numpy as np

from mouth_contact_age import emit


def rows(folder, key):
    if (json.loads((folder/'outcome.json').read_text())['exitCode']
            or json.loads((folder/'child-outcome.json').read_text())['exitCode']
            or not json.loads((folder/'limit-validation.json').read_text())['verified']):
        raise ValueError('Unverified prerequisite')
    return [row for line in (folder/'stdout').read_text().splitlines()
            if (row := json.loads(line)).get('test')==key]


def classify(folder):
    trials = rows(folder,'opening_acceleration')
    fixed, = [r for r in trials if r['mode']=='fixed']
    dc, h = (np.array(fixed[name]) for name in ('dc','derivative'))
    loss = -dc[1,1]
    if loss<=0:
        raise ValueError('The resistive support denominator needs its own pole test')
    # x=L/(L+Z). The material force per acceleration is a quadratic in x.
    coupling = dc[1,0]/loss
    polynomial = np.polynomial.Polynomial([h[0,0], h[0,1]*coupling+dc[0,1]*h[1,0]/loss,
                                           dc[0,1]*h[1,1]*coupling/loss])
    roots = sorted(float(r.real) for r in polynomial.roots()
                   if abs(r.imag)<1e-12 and 0<r.real<1)
    boundaries = [0., *roots, 1.]
    for row in trials:
        impedance = row['impedance']
        x = 0. if impedance is None else loss/(loss+impedance)
        np.testing.assert_allclose(-polynomial(x), row['signed_actuator_mass'],rtol=1e-10)
    emit(test='support_range',state=fixed['state'],mesh=fixed['mesh'],rays=fixed['rays'],
         force_polynomial=polynomial.coef.tolist(),x_roots=roots,bulk_drag_scale=loss,
         critical_support_drag=[loss*(1/x-1) for x in roots],
         critical_drag_over_scale=[1/x-1 for x in roots],
         actuator_sign_on_increasing_x=[int(np.sign(-polynomial((a+b)/2)))
                                        for a,b in zip(boundaries,boundaries[1:])])


def partition(folder):
    trials = [r for r in rows(folder,'opening') if r['feed_grid'][:2]==[96,6]
              and r['policy']=='reconnect']
    def select(distance,quadrature):
        row, = [r for r in trials if r['d']==distance and r['method']==quadrature]
        return row
    small, large, spherical = select(8,'cell'), select(16,'cell'), select(8,'spherical')
    emit(test='partition',state=small['state'],mesh=small['clock_mesh'],rays=small['rays'],
         mouth_brane_export_percent=[100*small[c]/small['carrier_total'] for c in ('mouth','brane','export')],
         ratio_16_to_8=large['mouth']/small['mouth'],
         quadrature_relative_shift=spherical['mouth']/small['mouth']-1,
         cell_mass=small['integrated_mass'],spherical_mass=spherical['integrated_mass'])


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--acceleration',type=Path,nargs='+',required=True)
    parser.add_argument('--pair',type=Path,nargs='+',required=True)
    args = parser.parse_args()
    for folder in args.acceleration:
        classify(folder)
    for folder in args.pair:
        partition(folder)
