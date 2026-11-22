"""Linear bare order load of a localized site on the fixed cubic mouth.

Use the unweighted energy e = kappa |grad chi|²/2 + W(chi). For the
medium's outward normal n, T n with T = e I - kappa grad chi ⊗ grad chi
is the external traction required by this sector. The force it exerts on
the wall has the opposite sign. The homogeneous Neumann condition leaves
only the normal term at the wall. A dipole pressure still has a nonzero
projection along the source's brane axis. A second observable translates
the Gaussian center with its confining wall, including their interaction
energy. Both exclude carrier stress, bath pressure response and mount.
"""
import argparse
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_carrier_lateral_response import trace
from mouth_curved_site_front import WIDTH, prepare, sha
from mouth_empty_bend_stiffness import quadrature
from mouth_flared_receiver import EXIT, FLARE_START, ROUNDING, emit
from mouth_free_order_diffuse import NU, SIGMA


def line_rule(edges, degree):
    nodes, weights = leggauss(degree)
    half = np.diff(edges)[:, None] / 2
    return (edges[:-1, None] + half * (1 + nodes)).ravel(), (half * weights).ravel()


def wall_integral(dom, fields, block, interval, kappa, delta):
    edges = dom.blocks[block][0]
    lo, hi = interval
    p, dp = line_rule(np.r_[lo, edges[(lo < edges) & (edges < hi)], hi], 12)
    r, _, rp, _, wp, _, values, grad_r, grad_w = trace(dom, block, p, np.zeros_like(p), fields)
    # p increases from feed to bulk. The reflecting solid lies outside the
    # medium in direction (w_p, -r_p) on both charts.
    ds = np.hypot(rp, wp)
    if np.any(wp <= 0):
        raise ValueError('Wall normal points away from chosen orientation')
    nr, nw = wp/ds, -rp/ds
    dt = (rp[:, None]*grad_r + wp[:, None]*grad_w) / ds[:, None]
    dn = nr[:, None]*grad_r + nw[:, None]*grad_w
    chi, eta = values.T
    well = kappa * chi*(1-chi)*(1-2*chi)*eta / delta**2
    gradient = kappa * dt[:, 0]*dt[:, 1]
    # ∫_{S²} cos²(theta) dΩ = 4π/3. The base has no angular gradient.
    measure = 4*np.pi/3 * r*r * dp
    natural = measure @ (wp*(gradient + well))
    # Full, unprojected Q2 stress trace. Normal derivatives are natural
    # only weakly on this mesh; this value diagnoses boundary accuracy.
    normal_product = kappa * dn[:, 0]*dn[:, 1]
    anisotropic = kappa * (grad_r[:, 0]*dn[:, 1] + (grad_r[:, 1]+2*eta/r)*dn[:, 0])
    raw = natural + measure @ (wp*normal_product - ds*anisotropic)
    return dict(holding=float(natural), gradient=float(measure @ (wp*gradient)),
                well=float(measure @ (wp*well)), full_q2_trace=float(raw))


def reservoir_integral(dom, fields, block, fixed_p, coordinate, normal, kappa):
    edges = dom.blocks[block][1 if fixed_p else 0]
    t, dt = line_rule(edges, 8)
    fixed = np.full_like(t, coordinate)
    p, v = (fixed, t) if fixed_p else (t, fixed)
    r, _, rp, rv, wp, wv, values, gr, gw = trace(dom, block, p, v, fields)
    ds = np.hypot(rv, wv) if fixed_p else np.hypot(rp, wp)
    nr, nw = normal
    # Fixed, spatially constant orders at the baths imply η=0 and that
    # the only nonzero gradient is normal. Then δ(T n)_x = -κ χ_n η_n n_x.
    np.testing.assert_allclose(values[:, 1], 0., atol=1.e-13, rtol=0.)
    gradients = nr*gr + nw*gw
    return float(-kappa * nr * (4*np.pi/3 * r*r * dt * ds)
                 @ (gradients[:, 0] * gradients[:, 1]))


def bulk_advection(dom, points_r, gr, gw, carrier, fields, mobility):
    phase_field = dom.N @ carrier
    velocity_r = ((gr @ (carrier-1))/phase_field).imag
    velocity_w = ((gw @ (carrier-1))/phase_field).imag
    grad_r, grad_w = gr @ fields, gw @ fields
    # div(T) = μ grad chi. At fixed carrier and E_int=-h∫f chi,
    # μ0 = -u.grad chi/M, (δμ)_1 = f1 - u.grad eta/M.
    mu_adv = -(velocity_r[:, None]*grad_r + velocity_w[:, None]*grad_w)/mobility
    eta = dom.N @ fields[:, 1]
    body = mu_adv[:, 1]*grad_r[:, 0] + mu_adv[:, 0]*(grad_r[:, 1] + 2*eta/points_r)
    return float((dom.weight/3) @ body)


def run(mesh, fields_path):
    dom, carrier, chi = prepare(mesh)
    kappa = 6*SIGMA*WIDTH
    mobility = NU/kappa
    fields_hash = sha(fields_path)
    with np.load(fields_path) as data:
        np.testing.assert_array_equal(dom.points, data['points'])
        np.testing.assert_array_equal(chi, data['chi'])
        r, _, gr, gw = quadrature(dom)
        for location in ('entry', 'half'):
            fields = np.column_stack((data['chi'], data[f'eta_1_{location}']))
            if not np.isfinite(fields).all():
                raise ValueError('Nonfinite dipole field')
            pieces = dict(rounded_lip=wall_integral(dom, fields, 0, (dom.blocks[0][0][0], 0.), kappa, WIDTH),
                          collar=wall_integral(dom, fields, 2, (ROUNDING, FLARE_START), kappa, WIDTH),
                          flare=wall_integral(dom, fields, 2, (FLARE_START, EXIT), kappa, WIDTH))
            # Other reflecting surfaces have n_x=0 pointwise. The axis is
            # regular, and chart seams in the reservoir are interior.
            baths = dict(feed=reservoir_integral(dom, fields, 1, True, dom.far, (1., 0.), kappa),
                         bulk_side=reservoir_integral(dom, fields, 4, False, 1., (1., 0.), kappa),
                         roof_inner=reservoir_integral(dom, fields, 3, True, dom.roof, (0., 1.), kappa),
                         roof_outer=reservoir_integral(dom, fields, 4, True, dom.roof, (0., 1.), kappa))
            source = float(data[f'source_dipole_grad_{location}'])
            adv = bulk_advection(dom, r, gr, gw, carrier, fields, mobility)
            reconstructed = adv + source - sum(baths.values())
            holding = sum(x['holding'] for x in pieces.values())
            # Holding the spatial order fixed, translate the wall and the
            # center together. In translated coordinates the truncated,
            # normalized Gaussian is stationary while chi(y + X) varies.
            # Its explicit energy derivative per h is -∫f ∂x chi, including
            # the moving boundary and normalization. A rigid translation
            # over the source support suffices up to Gaussian tail errors.
            composite = holding - source
            emit(test='linear_bare_order_wall_load', mesh=mesh, location=location,
                 kappa=kappa, mobility=mobility, delta=WIDTH, fields_sha256=fields_hash,
                 positive_axis='from mouth axis toward site', per_unit='positive order-favoring h',
                 direct_external_holding=holding, direct_medium_on_wall=-holding,
                 wall_pieces_holding=pieces, bulk_source=source, bulk_advection=adv,
                 reservoir_outward_flux=baths, reconstructed_holding=reconstructed,
                 direct_minus_reconstructed_holding=holding-reconstructed,
                 attached_source_interaction_holding=-source,
                 held_wall_and_attached_source=composite,
                 medium_on_wall_and_attached_source=-composite,
                 reconstructed_wall_and_attached_source=adv-sum(baths.values()))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=int, choices=(26, 34, 48), required=True)
    parser.add_argument('--fields', type=Path, required=True)
    args = parser.parse_args()
    run(args.mesh, args.fields)
