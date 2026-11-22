"""Stationary carrier momentum through nested cuts of the held cubic throat and receiver.

The drive is a weak lateral phase gradient in the remote brane; the remote
receiver has no imposed stream. Results apply to the reflecting material wall
of the carrier solver, with ordinary tangential slip and no extra tag potential.
"""
import argparse
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_carrier_lateral_response import Dipole, prepare, surface, surfaces, trace
from mouth_flared_receiver import EXIT, EXIT_RADIUS, emit
from mouth_segment_contact_hold import PITCH


def facet(dom, fields, block, along, fixed, normal, clip=None):
    """Physical angular integral of the advection and both stress contributions."""
    grid = dom.blocks[block][0 if along == 'p' else 1]
    if clip is not None:
        lo, hi = clip
        grid = np.r_[lo, grid[(grid > lo) & (grid < hi)], hi]
    roots, weights = leggauss(8)
    half = np.diff(grid)[:, None]/2
    s = (grid[:-1, None] + half*(1 + roots)).ravel()
    d = np.full_like(s, fixed)
    p, v = (s, d) if along == 'p' else (d, s)
    r, _, rp, rv, wp, wv, psi, dr, dw = trace(dom, block, p, v, fields)
    tr, tw = (rp, wp) if along == 'p' else (rv, wv)
    sign = np.sign(tw*normal[0] - tr*normal[1])
    if np.any(sign == 0):
        raise ValueError('Undefined outward normal')
    nr, nw = sign*tw, -sign*tr
    f, h = psi.T
    fr, hr = dr.T
    fw, hw = dw.T
    n = abs(f)**2
    amp = abs(f)
    radial, axial = fr/f, fw/f
    g = h/f
    gr, gw = (hr-g*fr)/f, (hw-g*fw)/f
    q = 4*np.pi/3 * r*r*(half*weights).ravel()
    un = nr*radial.imag + nw*axial.imag
    # d(v_x) weighted by the unperturbed normal flux. The angular
    # average of (1-mu²) is twice the average of mu² on an S² section.
    drift = q @ (n*un*(gr.imag + 2*g.imag/r))
    reweighted = q @ (n*radial.imag*(nr*gr.imag + nw*gw.imag + 2*un*g.real))
    ar, aw = amp*radial.real, amp*axial.real
    br, bw = ar*g.real + amp*gr.real, aw*g.real + amp*gw.real
    density_gradient = q @ ((br+2*amp*g.real/r)*(nr*ar+nw*aw) + ar*(nr*br+nw*bw))
    # Cross-check with the unsplit complex-field stress (same quadrature).
    total, isotropic = surface(dom, fields, block, along, fixed, normal, clip)
    np.testing.assert_allclose(total, drift+reweighted+density_gradient+isotropic, rtol=1e-9, atol=2e-12)
    return np.array([total, drift+reweighted, density_gradient, isotropic, drift, 3*q@(n*un)])


def read(mesh):
    folder = Path(f'_scratch/s11c/mouth_flared_cubic_n{mesh}')
    dom, f = prepare(SimpleNamespace(source=folder))
    model = Dipole(dom, f)
    speeds = np.array([1.e-4, -1.e-4, .5e-4, -.5e-4])
    cut_z = np.array([1.175, 1.4275, 2.4, EXIT-.015])
    emit(test='definitions', mesh=mesh, drives=speeds.tolist(), exit=EXIT,
         component_order=['outward_lateral_momentum','advection','density_gradient_stress',
                          'isotropic_stress','baseline_current_times_delta_vx','baseline_outward_current'],
         note='All first five components divided by signed far-brane stream speed; sixth is undriven number current.')
    for slot, U in enumerate(speeds):
        h, residual = model.steady(U, 0.)
        fields = np.column_stack((f, h/U))
        if slot == 0:
            reference = fields[:,1].copy()
        covariance = float(np.max(abs(fields[:,1]-reference)))
        wall, _, _ = surfaces(dom, fields)
        pressure = {k:v[1] for k,v in wall.items()}
        direct = sum(pressure.values())
        emit(test='signed_holder', speed=U, prescribed_rigid_external_force=-direct*U,
             plus_x_carrier_on_solid_per_speed=direct, pressure_per_speed=pressure,
             incremental_pde_residual=residual, normalized_field_difference=covariance)
        if covariance > 1.e-9:
            raise ArithmeticError('Signed solutions to the linear carrier equation disagree')
        if slot:
            continue
        for z in cut_z:
            whole = facet(dom, fields, 2, 'v', z, (0.,1.))
            ring = facet(dom, fields, 2, 'v', z, (0.,1.), (0., PITCH/float(dom.radius(z)[0])))
            flare = surface(dom, fields, 2, 'p', 0., (1.,0.), (z, EXIT))[1]
            emit(test='section', depth=z, whole=whole.tolist(), inner_one_pitch_at_wall=ring.tolist(),
                 downstream_surface_pressure_per_speed=flare,
                 drift_weighted_by_current_per_far_speed=whole[4]/whole[5],
                 near_wall_drift_weighted_by_current_per_far_speed=ring[4]/ring[5])
        for radius in (3.,6.,dom.far):
            rv = (radius-EXIT_RADIUS)/(dom.far-EXIT_RADIUS)
            feed = facet(dom, fields, 1, 'v', radius, (1.,0.))
            for z in (EXIT+.025, EXIT+.5, 5., 7., dom.roof):
                inner = facet(dom, fields, 3, 'v', z, (0.,1.))
                outer = facet(dom, fields, 4, 'v', z, (0.,1.), (0.,rv))
                side = facet(dom, fields, 4, 'p', rv, (1.,0.), (EXIT,z))
                received = inner+outer+side
                load = -feed[0]-received[0]
                emit(test='receiver_budget', radius=radius, depth=z, feed=feed.tolist(),
                     receiving_plane=(inner+outer).tolist(), receiving_side=side.tolist(),
                     receiving_total=received.tolist(), full_boundary_inferred_load_per_speed=load,
                     mismatch_from_direct_surface_pressure=load-direct,
                     reference_number_balance=feed[-1]+received[-1])


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh',type=int,required=True,choices=(26,34,48))
    read(parser.parse_args().mesh)
