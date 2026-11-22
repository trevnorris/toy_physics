"""Probe displacements between saved samples for two weak trial occupations.

Only the axial order-profile coordinate is released. Carrier, geometry, order
shape family and the mode action remain prescribed. Positive generalized drive
advances the front for the passive projection with positive mobility.
"""

import json
from types import SimpleNamespace

import numpy as np

from mouth_flared_receiver import FlaredDomain
from mouth_front_empty_order_work import integrate
from mouth_order_front_rim_reaction import mode
from mouth_order_front_wave_force import Front
from mouth_order_wave_pulse import ARCHIVE, ROOT, digest, emit
from mouth_shear_pattern_release import prepared, truncate


def main():
    args = SimpleNamespace(radius=1., state=None, mesh=34, sub=2, feed_spacing=.09,
                           cut=2.15, degree=4, outer=12.)
    guide, order, intake, path = prepared(args)
    front = Front(guide, order)
    carrier_path = ROOT/'mouth_flared_cubic_n34'/'field.npz'
    archive = json.loads(ARCHIVE.read_text())
    assert digest(carrier_path) == archive['input_carriers_sha256'][str(carrier_path)]
    cfg = json.loads((carrier_path.parent/'stdout').read_text().splitlines()[0])
    domain = FlaredDomain(cfg['n'], cfg['profile'], cfg['roof'], cfg['feed_far'])
    with np.load(carrier_path) as data:
        np.testing.assert_array_equal(domain.points, data['points'])
        carrier = data['psi']
    for index, block in domain.blocks.items():
        for native, extended in zip(block, guide.order_blocks[index], strict=True):
            np.testing.assert_array_equal(native, extended)
    # Disjoint from all five offsets used for the empty and occupied slopes.
    positions = (.001, .003, .006)
    fields = {q:front.fields(q) for q in positions}
    totals = integrate(domain, fields, carrier, 6, .04)
    emit(test='preparation', intake=intake, order_sha256=digest(path),
         carrier_sha256=digest(carrier_path), force_units='per unit order tension',
         occupations_over_tension=[.01,.02], candidate_positions=positions)
    for q in positions:
        def frequency(displacement):
            values, _ = front.fields(displacement)
            reduced, chi, cut = truncate(guide, values, args.cut)
            return mode(reduced, chi, cut, args, False)['frequency']

        readings = []
        for h in (.0002, .0001):
            lower, upper = frequency(q-h), frequency(q+h)
            readings.append(-(upper-lower)/(2*h))
        wave = readings[-1]
        bare, advective = totals[q][2:4]
        emit(test='weak_occupation', q=q, energy_slopes_per_action=readings,
             bare_per_tension=bare, advection_per_tension=advective,
             net_drive_per_tension={str(eta):bare+advective+eta*wave for eta in (.01,.02)})
        if abs(readings[1]-readings[0]) > .0001:
            raise ArithmeticError('Wave energy slopes disagree at an unused displacement')


if __name__ == '__main__':
    main()
