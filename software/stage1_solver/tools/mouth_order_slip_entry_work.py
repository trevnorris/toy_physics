"""First-excursion energy test for the local equilibrium order-wake reduction.

Use the earlier one-coordinate adsorption model with reversible static dressing
and steady local wake drag. The result is a necessary static order excursion
at or before chi=.8, not a time evolution of the coupled front.
"""
import argparse

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_bound_site_carrier_work import phase_integral
from mouth_calibrated_flare_coating import read_inputs
from mouth_free_order_diffuse import NU
from mouth_middle_tag_grip import cubic_trace
from mouth_surface_site_transfer import BASIS, Surface
from neutral_patch_sites import emit


def run(mesh, side):
    trace, intake = cubic_trace(mesh)
    patch, _, _ = read_inputs(mesh)
    cuts = np.array([1.175, trace.crossing(.80)])
    order = trace.sample(trace.chi, cuts)[0]
    # Check that the whole interval lies on the positive-curvature side. The
    # derivative of the order interpolant is affine on each quadratic element.
    knots = np.r_[cuts[0], trace.edges[(trace.edges>cuts[0]) & (trace.edges<cuts[1])], cuts[1]]
    one_sided = np.r_[knots[:-1], np.nextafter(knots[1:],knots[:-1])]
    assert max(trace.sample(trace.chi, one_sided)[1]) < 0 and min(1-6*order*(1-order))>0
    theta = np.diff(np.unwrap(np.angle(trace.sample(trace.psi,cuts)[0]))).item()
    np.testing.assert_allclose(phase_integral(trace,knots,8), theta, rtol=1.e-9, atol=1.e-12)
    launch = trace.speed(cuts[:1])[0]*np.hypot(1.,trace.geometry.radius(cuts[:1])[1][0])
    kinetic = launch**2/2
    np.testing.assert_allclose(intake,.25, rtol=0.,atol=1.e-8)
    emit(test='segment', mesh=mesh, side=side, depths=cuts.tolist(), order=order.tolist(),
         intake=intake, max_one_sided_order_derivative=max(trace.sample(trace.chi,one_sided)[1]),
         min_scaled_local_curvature=min(1-6*order*(1-order)), launch_speed=launch,
         launch_kinetic=kinetic, carrier_phase_rise=theta, nu=NU)
    surface = Surface(side,5.,patch['pitch'],False)
    start = np.zeros(surface.dim)
    start[:2] = patch['pitch']*BASIS.sum(axis=0)/3
    start[2] = .75
    scale = 2*patch['shared_scale']
    nodes, weights = leggauss(6)
    order_path = np.r_[order[0], (order.sum()+(order[1]-order[0])*nodes)/2, order[1]]
    for averaged, law in [(False,'weaken_all_contacts'),(True,'retain_mean_attraction')]:
        x = start.copy()
        rows = []
        for chi in order_path:
            x, energy, slope, info = surface.equilibrate(chi,averaged,x)
            rows.append(dict(chi=chi, mu=scale*energy, slope=scale*slope,
                             min_hessian=surface.stability(x,chi,averaged), max_force=info['max_force']))
        for row in reversed(rows[:-1]):
            x, e_rev, _, _ = surface.equilibrate(row['chi'],averaged,x)
            np.testing.assert_allclose(scale*e_rev,row['mu'],rtol=0,atol=5.e-8)
        slopes = np.array([r['slope'] for r in rows])
        assert max(slopes)<0 and min(np.diff(slopes))>=-2.e-7 and min(r['min_hessian'] for r in rows)>0
        rise = rows[-1]['mu']-rows[0]['mu']
        integrated = (order[1]-order[0])/2*(weights @ slopes[1:-1])
        np.testing.assert_allclose(integrated, rise, atol=2.e-8,rtol=0.)
        work_per_shift = abs(slopes[0])*theta/NU
        dressing_per_shift = abs(slopes[-1])/2
        needed = (rise-kinetic)/(work_per_shift+dressing_per_shift)
        emit(test='local_equilibrium_first_excursion', law=law, branch_samples=rows,
             bare_specific_energy_rise=rise, integrated_specific_energy_rise=integrated,
             max_steady_stream_work_per_order_shift=work_per_shift,
             max_static_dressing_credit_per_order_shift=dressing_per_shift,
             deficit_at_order_shift_001=rise-kinetic-.01*(work_per_shift+dressing_per_shift),
             necessary_max_site_weighted_static_shift=needed)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh',type=int,choices=(26,34),default=34)
    parser.add_argument('--side',type=int,choices=(5,8),default=5)
    args=parser.parse_args()
    run(args.mesh,args.side)
