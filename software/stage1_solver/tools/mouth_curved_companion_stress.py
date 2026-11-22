"""Mixed site-strength/incident-flow derivative of scalar-order force.

Write chi_total = chi + h a + U b cos(theta) + h U c cos(theta). The
isotropic population has total strength h and translates with the wall.
Positive x follows the incident phase gradient. All force integrals use
physical volume and the unweighted T = e I - kappa grad(chi) ⊗ grad(chi).
"""
import numpy as np

from mouth_carrier_lateral_response import trace
from mouth_curved_site_front import WIDTH
from mouth_curved_site_stress import line_rule, reservoir_integral, wall_integral
from mouth_empty_bend_stiffness import quadrature
from mouth_flared_receiver import EXIT, FLARE_START, ROUNDING
from mouth_free_order_diffuse import NU, SIGMA


def wall(dom, fields, block, along_p, fixed, outward, interval, kappa):
    edges = dom.blocks[block][0 if along_p else 1]
    if interval is not None:
        lo, hi = interval
        edges = np.r_[lo, edges[(lo < edges) & (edges < hi)], hi]
    q, dq = line_rule(edges, 12)
    p, v = (q, np.full_like(q, fixed)) if along_p else (np.full_like(q, fixed), q)
    r, _, rp, rv, wp, wv, values, gr, gw = trace(dom, block, p, v, fields)
    tr, tw = (rp, wp) if along_p else (rv, wv)
    ds = np.hypot(tr, tw)
    projection = tw*outward[0] - tr*outward[1]
    if np.any(abs(projection) < 1e-12):
        raise ValueError('Ambiguous outward wall normal')
    nr, nw = np.sign(projection)*tw/ds, -np.sign(projection)*tr/ds
    tangent = (tr[:, None]*gr + tw[:, None]*gw)/ds[:, None]
    gn = nr[:, None]*gr + nw[:, None]*gw
    chi, a, b, c = values.T
    tc0, ta, tb, tc = tangent.T
    # n·grad = 0 for each physical coefficient. n_x = n_r cos(theta)
    # and the mixed pressure carries another cos(theta).
    measure = 4*np.pi/3 * r*r * ds * dq
    weight = measure*nr
    curvature = 1 - 6*chi*(1-chi)
    linear = float(weight @ (kappa*(tc0*tc + chi*(1-chi)*(1-2*chi)*c/WIDTH**2)))
    cross = float(weight @ (kappa*(ta*tb + curvature*a*b/WIDTH**2)))
    # The full finite-element stress trace retains the normal derivatives
    # that the natural condition enforces only weakly on reflecting walls.
    normal_energy = gn[:, 0]*gn[:, 3] + gn[:, 1]*gn[:, 2]
    anisotropy = (gr[:, 0]*gn[:, 3] + (gr[:, 3]+2*c/r)*gn[:, 0] +
                  gr[:, 1]*gn[:, 2] + (gr[:, 2]+2*b/r)*gn[:, 1])
    raw = linear + cross + kappa*float(measure @ (nr*normal_energy - anisotropy))
    return dict(holding=linear+cross, base_with_mixed=linear, source_with_flow=cross,
                full_q2_trace=raw)


def shell_source_gradient(dom, b, depth, width, degree):
    """Return ∫ f0 ∂x(b cos theta) and the kernel's pre-normalization mass."""
    pg, vg, _ = dom.blocks[2]
    lo, hi = depth-12*width, depth+12*width
    if not pg[0] < lo < hi < pg[-1]:
        raise ValueError('Shell support extends outside the conduit chart')
    # Use the same entire axial elements as source() in mouth_curved_site_front.
    edges = pg[np.searchsorted(pg, lo, side='right')-1:np.searchsorted(pg, hi)+1]
    p, dp = line_rule(edges, degree)
    v, dv = line_rule(vg, degree)
    p_grid, v_grid = np.meshgrid(p, v, indexing='ij')
    r, _, _, rv, _, _, val, gr, _ = trace(dom, 2, p_grid.ravel(), v_grid.ravel(), b[:, None])
    center = dom.radius(depth)[0]
    alpha = r*center/width**2
    # Angular mean of the 4D Gaussian about any point on the S² shell.
    # Integrating exp[-alpha(1-cos theta)] over cos theta yields this f0.
    gaussian = np.exp(-((r-center)**2+(p_grid.ravel()-depth)**2)/(2*width**2))
    f0 = gaussian * (-np.expm1(-2*alpha)/(2*alpha)) / (2*np.pi*width**2)**2
    weights = 4*np.pi*r*r*abs(rv)*np.outer(dp, dv).ravel()*f0
    norm = float(weights.sum())
    # ∂x(b cos theta) = b_r cos²(theta) + b/r sin²(theta).
    grad = float(weights @ ((gr[:, 0] + 2*val[:, 0]/r)/3) / norm)
    return grad, norm


def anchor_gradient(dom, b, depth, width=.01):
    """Physical ∫ f0 ∂x(b P1) for a translating, normalized isotropic shell."""
    return shell_source_gradient(dom, b, depth, width, 12)[0]


def volume_advection(dom, psi, dpsi, fields, mobility):
    """Expand -(u·grad chi) ∂x chi/M directly at physical quadrature points."""
    r, _, gr, gw = quadrature(dom)
    base = dom.N @ psi
    ratio = (dom.N @ dpsi)/base
    # The carrier's dipole phase coefficient is Im(delta psi / psi). Its
    # angular gradient does not contract with either monopole gradient.
    def speeds(derivative):
        log_gradient = (derivative @ (psi-1))/base
        delta_velocity = ((derivative @ dpsi)/base - ratio*log_gradient).imag
        return log_gradient.imag, delta_velocity
    ur, dur = speeds(gr)
    uw, duw = speeds(gw)
    values = dom.N @ fields
    radial, axial = gr @ fields, gw @ fields
    adv = ur[:, None]*radial + uw[:, None]*axial
    # The carrier responds to U alone. Its velocity derivative acts on
    # the empty order in the U coefficient and on a in the hU coefficient.
    adv[:, 2:] += dur[:, None]*radial[:, :2] + duw[:, None]*axial[:, :2]
    projected = radial.copy()
    projected[:, 2:] += 2*values[:, 2:]/r[:, None]
    # In this product each monopole (chi,a) pairs with a dipole (b,c).
    # The cos² average is 1/3, including the angular derivative above.
    integrals = (-dom.weight/(3*mobility)) @ (adv*projected[:, ::-1])
    return dict(zip(('adv_chi_grad_c', 'adv_a_grad_b', 'adv_b_grad_a',
                     'adv_c_grad_chi'), map(float, integrals), strict=True))


def budget(dom, psi, chi, a, b, c, dpsi, depth, width=.01):
    """Compute outward scalar traction and the force on the held wall and shell."""
    kappa = 6*SIGMA*WIDTH
    fields = np.column_stack((chi, a, b, c))
    if not np.isrealobj(fields) or not np.all(np.isfinite(fields)):
        raise ValueError('The scalar-order coefficients must be real and finite')
    segments = dict(rounded_lip=(0, (dom.blocks[0][0][0], 0.)),
                    collar=(2, (ROUNDING, FLARE_START)), flare=(2, (FLARE_START, EXIT)))
    surfaces = {name: (block, True, 0., (1., 0.), interval)
                for name, (block, interval) in segments.items()}
    # On the remaining reflecting surfaces n_x is zero. Include them to
    # diagnose weak Neumann errors in the unprojected Q2 stress trace.
    surfaces.update(feed_plane=(1, True, 0., (0., 1.), None),
                    feed_floor=(1, True, 1., (0., -1.), None),
                    inner_floor=(0, True, 1., (0., -1.), (-np.pi/2, dom.bottom_turn)),
                    bulk_plane=(4, False, EXIT, (0., -1.), None))
    pieces = {name: wall(dom, fields, *surface, kappa) for name, surface in surfaces.items()}
    holding = sum(piece['holding'] for piece in pieces.values())
    # Check both wall terms against the first-order wall readout at four
    # shifted fields. Its unweighted traction is quartic at most in h.
    step = 1.e-3/max(1., np.max(np.abs(a)))
    sampled = []
    for h in step*np.array([-2., -1., 1., 2.]):
        shifted = np.column_stack((chi+h*a, b+h*c))
        sampled.append(sum(wall_integral(dom, shifted, *seg, kappa, WIDTH)['holding']
                           for seg in segments.values()))
    wall_check = float(np.dot([1., -8., 8., -1.], sampled)/(12*step))
    reservoirs = dict(feed=(1, True, dom.far, (1., 0.)),
                      bulk_side=(4, False, 1., (1., 0.)),
                      roof_inner=(3, True, dom.roof, (0., 1.)),
                      roof_outer=(4, True, dom.roof, (0., 1.)))
    # On fixed uniform-order baths only the normal gradients survive in
    # T n. These are the two mixed products in -kappa n_x chi_n²/2.
    bath_pairs = {name: dict(base_with_mixed=reservoir_integral(dom, fields[:, [0, 3]], *face, kappa),
                             source_with_flow=reservoir_integral(dom, fields[:, [1, 2]], *face, kappa))
                  for name, face in reservoirs.items()}
    baths = {name: sum(pair.values()) for name, pair in bath_pairs.items()}
    source, norm = shell_source_gradient(dom, b, depth, width, 12)
    source_check, norm_check = shell_source_gradient(dom, b, depth, width, 9)
    adv_terms = volume_advection(dom, psi, dpsi, fields, NU/kappa)
    adv = sum(adv_terms.values())
    # div T = (h f0 - u·grad chi/M) grad chi. Translating the wall with
    # its normalized attached shell contributes -∫ f0 ∂x(b cos theta)
    # to the holding load from -h∫f0 chi, canceling the source body term.
    composite = holding - source
    reconstructed = adv - sum(baths.values())
    return dict(direct_external_holding=holding, wall_pieces_holding=pieces,
                full_q2_wall_holding=sum(piece['full_q2_trace'] for piece in pieces.values()),
                wall_stencil_holding=wall_check, wall_stencil_step=step,
                wall_minus_stencil=holding-wall_check, reservoir_outward_flux=baths,
                reservoir_gradient_pairs=bath_pairs, bulk_source=source,
                source_check_degree9=source_check, gaussian_fraction_inside=norm,
                gaussian_fraction_inside_degree9=norm_check, bulk_advection=adv,
                bulk_advection_terms=adv_terms, attached_source_interaction_holding=-source,
                held_wall_and_attached_source=composite, medium_on_wall_and_attached_source=-composite,
                reconstructed_holding=reconstructed+source,
                reconstructed_wall_and_attached_source=reconstructed,
                direct_minus_reconstructed_holding=composite-reconstructed)
