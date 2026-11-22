"""Power required for stationary passive order conversion on saved radius-one fields.

The carrier and translated order profile are read from the checked 4D preparation.
The lower bound uses a smooth order marker and the independently integrated number
flux between pure baths. The script measures costs on these fields only.
"""

import argparse
import json
from types import SimpleNamespace

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_curved_phase import quadratic, tensor
from mouth_flared_receiver import FlaredDomain
from mouth_free_order_diffuse import NU, SIGMA
from mouth_order_front_wave_force import Front
from mouth_order_wave_pulse import ARCHIVE, ROOT, digest, emit


def measure(mesh, quadrature):
    delta = .04
    mobility = NU / (6 * delta * SIGMA)
    archive = json.loads(ARCHIVE.read_text())
    source = ROOT / f'mouth_flared_cubic_n{mesh}' / 'field.npz'
    suffix = '_outer_residual' if mesh == 48 else ''
    layer = ROOT / f'mouth_free_order_cubic_n{mesh}_d04{suffix}' / 'order_seed1.55.npz'
    assert digest(source) == archive['input_carriers_sha256'][str(source)]
    assert digest(layer) == archive['diffuse_orders_sha256'][str(layer)]
    setup = json.loads((source.parent / 'stdout').read_text().splitlines()[0])
    domain = FlaredDomain(mesh, setup['profile'], setup['roof'], setup['feed_far'])
    with np.load(source) as saved:
        np.testing.assert_array_equal(saved['points'], domain.points)
        carrier = saved['psi']
    with np.load(layer) as saved:
        front = Front(SimpleNamespace(order_blocks=domain.blocks), saved['chi'])
    rate = setup['feed']
    np.testing.assert_allclose(carrier[domain.feed], 1., atol=1e-12, rtol=0)
    np.testing.assert_allclose(carrier[domain.bulk], np.exp(1j * setup['bias']), atol=1e-12, rtol=0)
    emit(test='preparation', mesh=mesh, quadrature=quadrature, order_sha256=digest(layer),
         carrier_sha256=digest(source), phase_bias=setup['bias'], feed=rate, bulk=setup['bulk'],
         mobility=mobility, tension=SIGMA, thickness=delta)

    x, weight = leggauss(quadrature)
    xi, eta = (arr.ravel() for arr in np.meshgrid(x, x, indexing='ij'))
    pval, pd = quadratic(xi)
    vval, vd = quadratic(eta)
    shape, dp, dv = tensor(pval, vval), tensor(pd, vval), tensor(pval, vd)
    fields = {q: front.fields(q)[0] for q in (0., .001, .003, .006)}
    for chi in fields.values():
        np.testing.assert_allclose(chi[domain.feed], 1., atol=1e-12, rtol=0)
        np.testing.assert_allclose(chi[domain.bulk], 0., atol=1e-12, rtol=0)
    totals = {q: np.zeros(3) for q in fields}
    for block, (p, v, cells) in domain.blocks.items():
        i, j = (a.ravel() for a in np.indices((len(p)-1, len(v)-1)))
        for start in range(0, len(i), 384):
            sel = slice(start, start + 384)
            a, b, ids = i[sel], j[sel], cells[sel]
            hp, hv = np.diff(p)[a, None]/2, np.diff(v)[b, None]/2
            r, _, rp, rv, wp, wv = domain.chart(block, p[a,None] + hp*(1+xi),
                                                  v[b,None] + hv*(1+eta))
            det = rp*wv - rv*wp
            wt = 4*np.pi*r**2*abs(det)*hp*hv*np.outer(weight, weight).ravel()

            def interpolate(field):
                local = field[ids]
                fp, fv = (local @ dp.T)/hp, (local @ dv.T)/hv
                return local @ shape.T, (wv*fp-wp*fv)/det, (rp*fv-rv*fp)/det

            psi, pr, pw = interpolate(carrier)
            n, ur, uw = abs(psi)**2, (pr/psi).imag, (pw/psi).imag
            for q, chi in fields.items():
                level, cr, cw = interpolate(chi)
                advective_rate = ur*cr + uw*cw
                marker_weight = 6*level*(1-level)
                totals[q] += [np.sum(wt*advective_rate**2/mobility),
                              np.sum(wt*(-n*advective_rate)*marker_weight),
                              np.sum(wt*mobility*(n*marker_weight)**2)]

    for q, (power, converted, denominator) in totals.items():
        bound = converted**2/denominator
        emit(test='conversion', q=q, loss=power, head_per_constituent=power/rate,
             loss_per_unit_tension=power/SIGMA, integrated_marker_conversion=converted,
             marker_norm=denominator, cauchy_lower_bound=bound,
             bound_from_bath_rate=rate**2/denominator, fraction_of_measured_loss=bound/power)
        if not (0 < bound <= power and abs(converted/rate - 1.) < .003):
            raise ArithmeticError('Conversion power bound or independent bath flux check failed')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=int, choices=(26, 34, 48), required=True)
    parser.add_argument('--quadrature', type=int, choices=(4, 6), required=True)
    args = parser.parse_args()
    measure(args.mesh, args.quadrature)
