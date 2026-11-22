"""Work requirements for a cold site crossing the cubic carrier's wall to half order.

The stationary guide supplies an ideal normal reaction. For constant specific
linear slip drag gamma, site work is gamma * (integral u_s ds - integral v_s² dt).
Only the first integral is available from the carrier. The rate bounds allow a
site to slip between the endpoints. A site comoving everywhere has zero slip work.
"""
import argparse

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_calibrated_flare_coating import read_inputs
from mouth_middle_tag_grip import cubic_trace
from neutral_patch_sites import BETWEEN, emit


def phase_integral(trace, edges, degree):
    nodes, weights = leggauss(degree)
    half_width = np.diff(edges)[:, None] / 2
    depths = edges[:-1, None] + half_width * (1 + nodes)
    slope = trace.geometry.radius(depths)[1]
    return float(np.sum(half_width * weights * trace.speed(depths) * (1 + slope**2)))


def run(mesh, supplied_rises):
    trace, intake = cubic_trace(mesh)
    depths = np.array([1.175, trace.crossing(.5)])
    order = trace.sample(trace.chi, depths)[0]
    if abs(intake - .25) > 1.e-8 or depths[1] <= depths[0] or abs(order[1] - .5) > 1.e-10:
        raise ValueError('Wrong carrier preparation or half-order cut')
    radius, slope = trace.geometry.radius(depths)
    speed_w = trace.speed(depths)
    kinetic = .5 * (1 + slope**2) * speed_w**2
    delta_k = float(kinetic[1] - kinetic[0])
    edges = np.r_[depths[0], trace.edges[(depths[0] < trace.edges) &
                                        (trace.edges < depths[1])], depths[1]]
    coarse, fine = (phase_integral(trace, edges, degree) for degree in (4, 8))
    phase = np.unwrap(np.angle(trace.sample(trace.psi, edges)[0]))
    phase_change = float(phase[-1] - phase[0])
    if fine <= 0 or not np.allclose([coarse, phase_change], fine, atol=2.e-12, rtol=1.e-10):
        raise ArithmeticError('Phase difference and refined wall quadratures disagree')
    emit(test='carrier_path', mesh=mesh, intake=intake, depths=depths.tolist(),
         shear_order=order.tolist(), radius=radius.tolist(), radius_slope=slope.tolist(),
         axial_speed=speed_w.tolist(), carrier_tangent_speed=(np.hypot(1, slope)*speed_w).tolist(),
         specific_kinetic_energy=kinetic.tolist(), delta_kinetic_energy=delta_k,
         integral_u_ds_degree4=coarse, integral_u_ds_degree8=fine, endpoint_phase_change=phase_change,
         specific_energy_unit='mouth radius squared per model time squared',
         integral_u_ds_unit='mouth radius squared per model time')

    if supplied_rises:
        cases = [('input', rise) for rise in supplied_rises]
    else:
        patch, _, _ = read_inputs(mesh)
        binding = 2 * patch['shared_scale'] * BETWEEN[0]
        cases = [('fixed_continuum_reference', binding * (order[0] - order[1])),
                 ('level_well_bottom', 0.)]
    for label, rise in cases:
        if not np.isfinite(rise):
            raise ValueError('Transfer work must be finite')
        matched = float(rise + delta_k)
        reach = float(rise - kinetic[0])
        # For positive drag a nontrivial transit has integral v_s² dt > 0.
        # The arrival kinetic energy can be zero in the reach-only comparison.
        emit(test='endpoints', mesh=mesh, transfer_law=label, delta_mu=rise,
             required_specific_work_at_carrier_endpoint_speeds=matched,
             reach_only_specific_work_floor=reach,
             rate_infimum_from_matched_speeds=max(0., matched/fine),
             rate_infimum_from_reach_only=max(0., reach/fine),
             condition='Positive constant specific drag must exceed the listed bound. '
                       'Zero bound adds no restriction. Endpoint conditions are necessary only.',
             rate_unit='inverse model time')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=int, required=True, choices=(26, 34))
    parser.add_argument('--delta-mu', type=float, nargs='+', action='extend',
                        help='Signed host-subtracted specific transfer work(s), in carrier units. '
                             'Default: fixed-continuum and level-bottom single-site controls.')
    args = parser.parse_args()
    run(args.mesh, args.delta_mu)
