"""Bare order work and passive advection on the wave's axial front displacement.

Configuration space is the native 4D throat with two angles integrated out.
Every force here is per unit order tension, positive for increasing q; the
longitudinal coherent velocity is held at its archived value. Integrals include
the complete physical order energy. Its derivative is supported in the conduit.
"""
import argparse
import json
from types import SimpleNamespace

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_curved_phase import quadratic, tensor
from mouth_flared_receiver import FlaredDomain
from mouth_free_order_diffuse import NU
from mouth_order_front_wave_force import Front
from mouth_order_wave_pulse import ARCHIVE, ROOT, digest, emit


def integrate(domain, fields, carrier, nquad, delta):
    nodes, weights = leggauss(nquad)
    xi, eta = (a.ravel() for a in np.meshgrid(nodes, nodes, indexing='ij'))
    ax, dx = quadratic(xi)
    ay, dy = quadratic(eta)
    shape, dp, dv = tensor(ax, ay), tensor(dx, ay), tensor(ax, dy)
    totals = {q: np.zeros(6) for q in fields}
    kappa, well = 6*delta, 3/delta  # per unit planar tension
    for block, (p, v, cells) in domain.blocks.items():
        i, j = (a.ravel() for a in np.indices((len(p)-1, len(v)-1)))
        for start in range(0, len(i), 512):
            sl = slice(start, start+512)
            a, b, ids = i[sl], j[sl], cells[sl]
            hp, hv = np.diff(p)[a, None]/2, np.diff(v)[b, None]/2
            r, _, rp, rv, wp, wv = domain.chart(block, p[a, None]+hp*(1+xi), v[b, None]+hv*(1+eta))
            det = rp*wv-rv*wp
            wt = 4*np.pi*r*r*abs(det)*hp*hv*np.outer(weights, weights).ravel()

            def value_grad(field):
                local = field[ids]
                fp, fv = (local@dp.T)/hp, (local@dv.T)/hv
                return local@shape.T, (wv*fp-wp*fv)/det, (rp*fv-rv*fp)/det

            psi, pr, pw = value_grad(carrier)
            ur, uw = (pr/psi).imag, (pw/psi).imag
            weight_phi = np.exp(-np.angle(psi)/NU)
            for q, (chi, tangent) in fields.items():
                c, cr, cw = value_grad(chi)
                d, dr, dw = value_grad(tangent)
                potential = well*(c*(1-c))**2
                gradient = 0.5*kappa*(cr*cr+cw*cw)
                de = 2*well*c*(1-c)*(1-2*c)*d + kappa*(cr*dr+cw*dw)
                # G(q) qdot = F_bare + F_adv under the passive projection.
                advection = -(kappa/NU)*(ur*cr+uw*cw)*d
                for index, density in enumerate((potential,gradient,-de,advection,
                                                 weight_phi*(potential+gradient),-weight_phi*de)):
                    totals[q][index] += np.sum(wt*density)
    return totals


def main(mesh, nquad):
    archive = json.loads(ARCHIVE.read_text())
    tag = '_outer_residual' if mesh == 48 else ''
    order_path = ROOT/f'mouth_free_order_cubic_n{mesh}_d04{tag}'/'order_seed1.55.npz'
    carrier_path = ROOT/f'mouth_flared_cubic_n{mesh}'/'field.npz'
    assert digest(order_path) == archive['diffuse_orders_sha256'][str(order_path)]
    assert digest(carrier_path) == archive['input_carriers_sha256'][str(carrier_path)]
    preparation = json.loads((carrier_path.parent/'stdout').read_text().splitlines()[0])
    dom = FlaredDomain(mesh, preparation['profile'], preparation['roof'], preparation['feed_far'])
    with np.load(carrier_path) as data:
        np.testing.assert_array_equal(dom.points, data['points'])
        carrier = data['psi']
    with np.load(order_path) as data:
        order = data['chi']
    front = Front(SimpleNamespace(order_blocks=dom.blocks), order)
    shifts = [0., -.02, -.01, .01, .02]
    fields = {q:front.fields(q) for q in shifts}
    for _, tangent in fields.values():
        assert not np.any(tangent[dom.feed | dom.bulk])
    emit(test='preparation', mesh=mesh, gauss_points=nquad, intake=preparation['feed'],
         carrier_sha256=digest(carrier_path), order_sha256=digest(order_path), delta=.04,
         nu=NU, tension_normalization=1., energy_by_quadrature='whole four-dimensional carrier domain')
    totals = integrate(dom, fields, carrier, nquad, .04)
    for q, (potential, gradient, fbare, fadv, weighted_energy, fweighted) in totals.items():
        emit(test='order_work', q=q, energy_per_tension=potential+gradient,
             force_bare_per_tension=fbare, force_advective_per_tension=fadv,
             projected_total_per_tension=fbare+fadv, weighted_energy_per_tension=weighted_energy,
             weighted_variational_force_per_tension=fweighted)
    for h in (.02,.01):
        low, high = totals[-h], totals[h]
        emit(test='secant', qstep=h, energy_force_per_tension=-sum((high-low)[:2])/(2*h),
             bare_tangent_stiffness_per_tension=-(high[2]-low[2])/(2*h),
             passive_projected_stiffness_per_tension=-sum((high-low)[2:4])/(2*h))
    # The original weighted nodal equilibrium was solved by 4-point quadrature.
    if nquad == 4 and abs(totals[0][5])>1e-7:
        raise ArithmeticError('Archived order field is not stationary on the shared displacement')
    if abs(totals[0][2]+totals[0][3])>.01*abs(totals[0][2]):
        raise ArithmeticError('Independently integrated advection does not balance empty order force')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', required=True, type=int, choices=(26,34,48))
    parser.add_argument('--quadrature', type=int, default=4, choices=(4,6))
    args=parser.parse_args()
    main(args.mesh,args.quadrature)
