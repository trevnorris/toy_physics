"""Match the leading static sink/source tails at fixed intake on archived mouth meshes.

The outer phases have open monopole asymptotics. Secondary multipoles, transverse
slab structure and quantum pressure beyond this cut are not prescribed exactly.
Interior hemisphere measurements are independent of the boundary phase data.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.sparse.linalg import splu

from mouth_carrier_lateral_response import prepare, trace, weighted_gram
from mouth_curved_phase import RATE
from mouth_flared_receiver import EXIT, EXIT_RADIUS, emit
from mouth_occupied_push_baths import load_bath


def equilibrium(dom, original):
    field = original.copy()
    interp, free = dom.N, dom.free
    half_grad = dom.stiff * .5
    f = interp @ original
    n = abs(f)**2
    interp_free = interp[:, free]
    # A fixed Jacobian at the archived solution preconditions the small
    # change to the boundary data. Iterate the full nonlinear residual.
    grad = half_grad[free][:, free] + weighted_gram(interp_free, dom.weight*(n**4-1)/4)
    slope = dom.weight*2*n**3
    cross = weighted_gram(interp_free, slope*f.real*f.imag)
    jac = sparse.bmat([[grad+weighted_gram(interp_free, slope*f.real**2), cross],
                       [cross, grad+weighted_gram(interp_free, slope*f.imag**2)]], format="csc")
    linear = splu(jac)
    feed_radius = dom.far
    r, w = dom.points[dom.bulk].T
    distance = np.hypot(r, w-EXIT)
    # 4πr² times unit slab thickness, and π²s³ for a 4D hemisphere.
    feed_phase = RATE/(4*np.pi*feed_radius)
    source_phase = -RATE/(2*np.pi**2*distance**2)
    feed_speed = RATE/(4*np.pi*feed_radius**2)
    source_speed = RATE/(np.pi**2*distance**3)
    # Leading stationary Bernoulli amplitude. Dispersion and multipoles
    # are tested by moving the truncating faces and sampling fixed inner cuts.
    field[dom.feed] = (1-2*feed_speed**2)**.125*np.exp(1j*feed_phase)
    original_bias = float(np.angle(original[dom.bulk][0]))
    bias = original_bias + feed_phase - np.mean(source_phase)
    derivative = RATE/original_bias
    previous = None
    for _ in range(10):
        field[dom.bulk] = (1-2*source_speed**2)**.125*np.exp(1j*(bias+source_phase))
        for _ in range(24):
            f = interp @ field
            residual = half_grad @ (field-1) + interp.T @ (dom.weight*(abs(f)**8-1)*f/4)
            error = float(np.max(abs(residual[free])))
            # Use the stationary carrier's established absolute residual tolerance.
            if error < 3e-11:
                break
            increment = linear.solve(-np.r_[residual[free].real, residual[free].imag])
            da, db = np.split(increment, 2)
            field[free] += da+1j*db
        else:
            raise ArithmeticError(f"Stationary interior residual {error}")
        flux = 2*(field.conj()*residual).imag
        feed = -flux[dom.feed].sum()
        bulk = flux[dom.bulk].sum()
        if abs(bulk-RATE) < 8e-11:
            emit(test="matched_static", radius=dom.far, roof=dom.roof, bias=bias,
                 boundary_phase_span=float(np.ptp(source_phase)), feed=feed, bulk=bulk,
                 residual=error, min_density=float(np.min(abs(f)**2)),
                 section=dom.section(field, 2.4))
            if abs(feed-bulk) > 3e-9 or np.min(abs(f)**2) < .5:
                raise ArithmeticError("Intake or finite-density check failed")
            return field, bias
        if previous is not None:
            derivative = (bulk-previous[1])/(bias-previous[0])
        previous = bias, bulk
        bias += (RATE-bulk)/derivative
    raise ArithmeticError("Intake calibration did not converge")


def hemisphere(dom, field, radius):
    # Split quadrature at each Q2 cell crossing, including the seam.
    inner_depth, inner_v, _ = dom.blocks[3]
    _, outer_v, _ = dom.blocks[4]
    radii = np.r_[EXIT_RADIUS*(1-inner_v), EXIT_RADIUS+(dom.far-EXIT_RADIUS)*outer_v]
    depth = inner_depth-EXIT
    edges = np.unique(np.r_[0., np.pi/2, np.arcsin(radii[(radii>0)&(radii<radius)]/radius),
                           np.arccos(depth[(depth>0)&(depth<radius)]/radius)])
    root, weights = leggauss(8)
    half = np.diff(edges)[:, None]/2
    angle = (edges[:-1,None] + half*(root+1)).ravel()
    r, z = radius*np.sin(angle), EXIT+radius*np.cos(angle)
    frad = np.empty_like(angle)
    phase = np.empty_like(angle)
    for block, mask in ((3, r<EXIT_RADIUS), (4, r>=EXIT_RADIUS)):
        v = 1-r[mask]/EXIT_RADIUS if block==3 else (r[mask]-EXIT_RADIUS)/(dom.far-EXIT_RADIUS)
        _, _, _, _, _, _, f, fr, fw = trace(dom, block, z[mask], v, field[:,None])
        derivative = fr[:,0]*np.sin(angle[mask]) + fw[:,0]*np.cos(angle[mask])
        frad[mask] = (f[:,0].conj()*derivative).imag
        phase[mask] = np.angle(f[:,0])
    measure = 4/np.pi*np.sin(angle)**2*(half*weights).ravel()
    np.testing.assert_allclose(measure.sum(), 1., atol=2e-13, rtol=0)
    normalized = frad*(np.pi**2*radius**3/RATE)
    return dict(radius=radius, flux=RATE*(measure@normalized),
                fractional_rms_from_monopole=np.sqrt(measure@((normalized-1)**2)),
                quadrupole_flux_over_monopole=measure@(normalized*(4*np.cos(angle)**2-1)),
                inferred_far_phase=measure@phase + RATE/(2*np.pi**2*radius**2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sources = parser.add_mutually_exclusive_group(required=True)
    sources.add_argument("--source", type=Path)
    sources.add_argument("--bath", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    folder = args.bath or args.source
    dom, field = load_bath(args.bath) if args.bath else prepare(args)
    emit(test="source", file=str(folder/'field.npz'),
         sha256=sha256((folder/'field.npz').read_bytes()).hexdigest())
    matched, bias = equilibrium(dom, field)
    cuts = []
    for label, f in (("archived",field), ("matched",matched)):
        for radius in (4.,6.):
            row = hemisphere(dom, f, radius)
            emit(test="inner_hemisphere", state=label, **row)
            cuts.append(row)
    if any(abs(row["flux"]-RATE)>.001 for row in cuts):
        raise ArithmeticError("Interior hemisphere flux differs from target")
    np.savez_compressed(args.output, points=dom.points, psi=matched, bias=bias)
    emit(test="stored", file=str(args.output), sha256=sha256(args.output.read_bytes()).hexdigest())


if __name__ == "__main__":
    main()
