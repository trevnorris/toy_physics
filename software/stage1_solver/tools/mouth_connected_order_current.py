"""Surface integrals for phase-compatible passive current across an order level.

The mouth cap comes from a saved curved conduit. An intact planar face is a
separate local preparation. Combining their weights is conditional on joining
them into one regular order level with the same reciprocal conversion law.
"""

import argparse
import json
from types import SimpleNamespace

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq

from mouth_curved_phase import quadratic
from mouth_flared_receiver import FlaredDomain
from mouth_free_order_diffuse import NU, SIGMA, read
from mouth_order_front_wave_force import Front
from mouth_order_wave_pulse import ARCHIVE, ROOT, digest, emit


def depth(domain, order, eta, level):
    # At fixed r/a(w), the order is a quadratic on each native depth cell.
    w, v, ids = domain.blocks[2]
    j = np.clip(np.searchsorted(v, 1-eta)-1, 0, len(v)-2)
    basis = quadratic(2*(1-eta-v[j])/(v[j+1]-v[j])-1)[0]
    grid = ids.reshape(len(w)-1, len(v)-1, 3, 3)[:, j]
    y = order[grid] @ basis - level
    a = (y[:, 0]+y[:, 2])/2-y[:, 1]
    b = (y[:, 2]-y[:, 0])/2
    roots = []
    for i, (left, center, right) in enumerate(y):
        # Split at any interior extremum to detect an extra crossing.
        edges = [-1., 1.]
        if abs(b[i]) < 2*abs(a[i]):
            edges.insert(1, -b[i]/(2*a[i]))
        poly = lambda t: (a[i]*t+b[i])*t + center
        for low, high in zip(edges[:-1], edges[1:]):
            if poly(low) == 0:
                root = low
            elif poly(low)*poly(high) < 0:
                root = brentq(poly, low, high)
            else:
                continue
            roots.append(w[i]+(root+1)*(w[i+1]-w[i])/2)
    if len(roots) != 1:
        raise ArithmeticError(f'Expected one regular crossing, eta={eta}, level={level}: {roots}')
    return roots[0]


def surface(domain, order, carrier, level, mobility, count):
    points, weights = leggauss(count)
    v = domain.blocks[2][1]
    # Integrate piecewise so derivative jumps at radial element boundaries are sampled correctly.
    eta_edges = 1-v[::-1]
    etas = (eta_edges[:-1,None] + np.diff(eta_edges)[:,None]*(points+1)/2).ravel()
    weights = (np.diff(eta_edges)[:,None]*weights/2).ravel()
    area = sheet_weight = number_flux = alternate_weight = 0.
    affinities = []
    for eta, weight in zip(etas, weights, strict=True):
        w = depth(domain, order, eta, level)
        radius, slope = domain.radius(w)
        r = eta*radius
        c, cr, cw = read(domain, order, r, w)
        psi, pr, pw = read(domain, carrier, r, w)
        n = abs(psi)**2
        derivative = cw + eta*slope*cr
        gradient = np.hypot(cr, cw)
        if abs(c-level)>1e-9 or derivative>=0 or n<=0:
            raise ArithmeticError('Non-regular falling order surface')
        dw = -radius*cr/derivative
        dr = radius+eta*slope*dw
        darea = weight*4*np.pi*r*r*np.hypot(dr,dw)
        area += darea
        sheet_weight += darea*mobility*n*n/gradient
        alternate_weight += weight*4*np.pi*r*r*radius*mobility*n*n/abs(derivative)
        adv = (pr/psi).imag*cr+(pw/psi).imag*cw
        number_flux += darea*(-n*adv/gradient)
        affinities.append(-adv/(mobility*n))
    np.testing.assert_allclose(sheet_weight, alternate_weight, rtol=2e-12)
    return dict(cap_three_area=area, cap_weight=sheet_weight,
                saved_flux=number_flux, cap_only_required_mu_over_n=number_flux/sheet_weight,
                saved_mu_over_n_range=[min(affinities), max(affinities)])


def main(mesh, count):
    manifest = json.loads(ARCHIVE.read_text())
    carrier_path = ROOT/f'mouth_flared_cubic_n{mesh}'/'field.npz'
    suffix = '_outer_residual' if mesh==48 else ''
    order_path = ROOT/f'mouth_free_order_cubic_n{mesh}_d04{suffix}'/'order_seed1.55.npz'
    assert digest(carrier_path)==manifest['input_carriers_sha256'][str(carrier_path)]
    assert digest(order_path)==manifest['diffuse_orders_sha256'][str(order_path)]
    cfg = json.loads((carrier_path.parent/'stdout').read_text().splitlines()[0])
    dom = FlaredDomain(mesh, cfg['profile'], cfg['roof'], cfg['feed_far'])
    with np.load(carrier_path) as data:
        np.testing.assert_array_equal(dom.points, data['points'])
        carrier = data['psi']
    with np.load(order_path) as data:
        front = Front(SimpleNamespace(order_blocks=dom.blocks), data['chi'])
    delta = .04
    mobility = NU/(6*SIGMA*delta)
    emit(test='preparation', mesh=mesh, quadrature=count, delta=delta, mobility=mobility,
         order_sha256=digest(order_path), carrier_sha256=digest(carrier_path),
         flat_kink_density=1., radius=1.)
    for shift in (0., .003, .006):
        field, _ = front.fields(shift)
        for c in (.25,.5,.75):
            result = surface(dom, field, carrier, c, mobility, count)
            flat_per_area = mobility*delta/(c*(1-c))
            equal_area = result['cap_weight']/flat_per_area
            emit(test='order_surface', q=shift, level=c, flat_weight_per_three_area=flat_per_area,
                 flat_three_area_for_equal_flow=equal_area,
                 flat_fraction_if_area_one_unit_ball=(4*np.pi/3)/(equal_area+4*np.pi/3), **result)
            if abs(result['saved_flux']-cfg['feed']) > 5e-4:
                raise ArithmeticError('Cap current and independent reservoir current disagree')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=int, choices=(26,34,48), required=True)
    parser.add_argument('--quadrature', type=int, default=6)
    args = parser.parse_args()
    main(args.mesh, args.quadrature)
