"""Directional carrier supply and inventory at two cuts of a held cubic mouth.

The radial band is R(w)-b < r < R(w) at a fixed physical depth. Its axial
current counts passage per unit time and its density counts inventory per unit
depth. Either can be sampled under a declared capture law; neither identifies
bound sites. Positive w heads into the mouth. The angular axis is +x, the
imposed feeding velocity. All linear coefficients refer to d/dU at U=0.
"""
import argparse
import json

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_carrier_lateral_response import Dipole, trace
from mouth_curved_site_front import front_read, prepare, sha
from mouth_curved_site_pair import RECORD as MATERIAL
from mouth_flared_receiver import emit
from mouth_segment_contact_hold import PITCH


def gauss(edges, degree):
    nodes, weights = leggauss(degree)
    half = np.diff(edges)[:, None]/2
    return (edges[:-1, None] + half*(1+nodes)).ravel(), (half*weights).ravel()


def section(dom, fields, depth, width, degree):
    radius = float(dom.radius(depth)[0])
    outer_v = 1. if width is None else width/radius
    nodes = dom.blocks[2][1]
    v, dv = gauss(np.r_[0., nodes[(nodes>0) & (nodes<outer_v)], outer_v], degree)
    r, _, _, rv, _, _, f, _, fw = trace(dom, 2, np.full_like(v, depth), v, fields)
    measure = 4*np.pi*r*r*abs(rv)*dv
    return f, fw, measure, abs(rv)*dv


def angular_coefficients(f, fw, measure):
    psi, h = f.T
    pw, hw = fw.T
    rate0 = (psi.conj()*pw).imag
    rate1 = (psi.conj()*hw+h.conj()*pw).imag
    n0 = abs(psi)**2
    n1 = 2*(psi.conj()*h).real
    return measure @ np.column_stack([rate0, rate1, n0, n1])


def signed_moments(f, fw, measure, drive):
    # Reconstruct the tangent field at signed amplitude and form products
    # before angular integration. This controls projection, polarity and the
    # even products of the known tangent; no nonlinear carrier solve is implied.
    nodes, weights = leggauss(5)
    edges = np.array([-1., -1/3, 1/3, 1.])
    x = (edges[:-1, None]+(1+nodes)*np.diff(edges)[:, None]/2).ravel()
    weights = (np.diff(edges)[:, None]*weights/4).ravel()
    field = f[:, :1]+drive*f[:, 1:]*x
    gradient = fw[:, :1]+drive*fw[:, 1:]*x
    current = (field.conj()*gradient).imag
    weighted = np.array([measure@current, measure@(abs(field)**2)])*weights
    total = weighted.sum(axis=1)
    moments = (weighted@x)/total
    thirds = weighted.reshape(2, 3, -1).sum(axis=2)/total[:, None]
    # At fixed r the reconstructed current is quadratic in cos(theta).
    # Bound its minimum for every angular direction at each radial node.
    c0 = (f[:,0].conj()*fw[:,0]).imag
    c1 = drive*(f[:,1].conj()*fw[:,0]+f[:,0].conj()*fw[:,1]).imag
    c2 = drive**2*(f[:,1].conj()*fw[:,1]).imag
    vertex = np.clip(np.divide(-c1, 2*c2, out=np.zeros_like(c1), where=c2>0), -1., 1.)
    minimum = np.minimum(c0+c2-abs(c1), c0+c1*vertex+c2*vertex**2).min()
    return total, moments, thirds, minimum


def check_drive(f, fw, measure, keys, beta, drive):
    plus, positive, front, jp = signed_moments(f, fw, measure, drive)
    minus, negative, back, jm = signed_moments(f, fw, measure, -drive)
    by_first_moment = 3*(positive-negative)/(2*drive)
    by_thirds = (9/4)*((front[:, 2]-front[:, 0])-(back[:, 2]-back[:, 0]))/(2*drive)
    emit(test='signed_band', **keys, drive=drive, quantities=['crossing_rate','inventory_per_depth'],
         tangent_products_plus_totals=plus.tolist(), tangent_products_minus_totals=minus.tolist(),
         plus_mean_cosine=positive.tolist(), minus_mean_cosine=negative.tolist(),
         plus_back_middle_front_thirds=front.tolist(), minus_back_middle_front_thirds=back.tolist(),
         from_three_mean_cosine=by_first_moment.tolist(), from_polar_thirds=by_thirds.tolist(),
         minimum_reconstructed_inward_current_at_radial_nodes=min(jp,jm))
    np.testing.assert_allclose(plus, minus, rtol=1e-13)
    np.testing.assert_allclose(front, back[:, ::-1], atol=1e-13)
    np.testing.assert_allclose(by_first_moment, by_thirds, rtol=1e-8, atol=1e-11)
    np.testing.assert_allclose(by_first_moment, beta, rtol=1e-4, atol=1e-11)
    if min(jp,jm)<=0:
        raise ValueError('Sampled inward crossings reverse at this signed drive')


def radial_side(dom, fields, depths, width):
    radius, slope = dom.radius(depths)
    v = width/radius
    r, _, _, _, _, _, f, fr, fw = trace(dom, 2, depths, v, fields)
    fn = fr-slope[:, None]*fw
    j0 = (f[:,0].conj()*fn[:,0]).imag
    j1 = (f[:,1].conj()*fn[:,0]+f[:,0].conj()*fn[:,1]).imag
    return (4*np.pi*r[:, None]**2)*np.column_stack([j0, j1])


def balance(dom, fields, depths, width, currents, keys):
    # Positive side terms point toward increasing r. Angular divergence
    # supplies -2 Im(psi* h)/r² in the l=1 equation. The axial change must
    # equal angular import plus interior radial import minus the wall leak.
    lo, hi = depths
    pg = dom.blocks[2][0]
    z, dz = gauss(np.r_[lo, pg[(pg>lo) & (pg<hi)], hi], 10)
    outer = dz@radial_side(dom, fields, z, 0.)
    inner = np.zeros(2) if width is None else dz@radial_side(dom, fields, z, width)
    angular = []
    for w in z:
        f, _, _, dr = section(dom, fields, w, width, 8)
        angular.append(8*np.pi*dr @ (f[:,0].conj()*f[:,1]).imag)
    import_flux = np.array([0., dz@angular])
    change = currents[1]-currents[0]
    defect = change+outer-inner-import_flux
    emit(test='axial_angular_balance', **keys, interval=depths, modes=['base','dipole_per_U'],
         axial_current_gain=change.tolist(), wall_leak=outer.tolist(), inner_import=inner.tolist(),
         angular_import=import_flux.tolist(), residual=defect.tolist())


def main(mesh):
    dom, psi, chi = prepare(mesh)
    dpsi, residual = Dipole(dom, psi).steady(1.,0.)
    fields = np.column_stack([psi, dpsi])
    np.testing.assert_allclose(dpsi[dom.feed], 1j*psi[dom.feed]*dom.points[dom.feed,0], atol=1e-12)
    np.testing.assert_allclose(dpsi[dom.bulk], 0., atol=1e-12)
    depths = [1.175, float(front_read(dom, chi, [], 1.)[0])]
    material = json.loads(MATERIAL.read_text())['guarded_runs']['primary']
    assert material['limits_verified'] and material['guard_exit']==material['child_exit']==0
    sample, = [row for line in material['literal_stdout'].splitlines()
               if (row:=json.loads(line))['test']=='sample']
    bands = {'one_pitch_unit':PITCH, 'relaxed_lattice_pitch':sample['pitch']*PITCH,
             'half_pitch_unit':.5*PITCH, 'two_pitch_units':2*PITCH, 'full':None}
    emit(test='carrier_setup', mesh=mesh, depths=depths, native_to_mouth_length=PITCH,
         calibrated_native_pitch=sample['pitch'], calibration_sha256=sha(MATERIAL),
         dipole_free_absolute_residual=residual, feed_gradient=1., receiver_gradient=0.,
         nominal_intake=.25, widths=bands)
    for band, width in bands.items():
        currents = []
        for w, name in zip(depths, ('entry','half'), strict=True):
            f, fw, measure, _ = section(dom, fields, w, width, 8)
            vals = angular_coefficients(f, fw, measure)
            fg, wg, mg, _ = section(dom, fields, w, width, 5)
            alternate = angular_coefficients(fg, wg, mg)
            rate, rate1, n, n1 = vals
            beta = np.array([rate1/rate, n1/n])
            currents.append(vals[:2])
            radius = float(dom.radius(w)[0])
            independent = dom.section(psi, w)['integrated_jw']
            keys = dict(mesh=mesh, location=name, depth=w, band=band, width=width)
            emit(test='band', **keys, radius=radius, outer_radial_element=radius*dom.blocks[2][1][1],
                 rate=rate, rate_dipole_per_U=rate1, rate_beta=beta[0],
                 inventory_per_depth=n, inventory_dipole_per_depth_U=n1, inventory_beta=beta[1],
                 full_section_current=independent, band_fraction=rate/independent,
                 quadrature_change=(vals-alternate).tolist())
            np.testing.assert_allclose(vals, alternate, atol=2e-13, rtol=1e-11)
            if width is None:
                np.testing.assert_allclose(rate, independent, atol=2e-13, rtol=0)
            for drive in (.002, .001):
                check_drive(f, fw, measure, keys, beta, drive)
        balance(dom, fields, depths, width, np.array(currents), dict(mesh=mesh, band=band,width=width))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', required=True, type=int, choices=(26,34,48))
    main(parser.parse_args().mesh)
