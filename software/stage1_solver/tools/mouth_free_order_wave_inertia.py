"""Matched axial shear pulses on the two archived free-order preparations.

Use the trial radially uniform scalar displacement u(w,t) in the rest-frame
shear action. K=integral chi dA, M=integral chi**m dA with dA=4*pi*r² dr.
The specimen runs from the common feed cut to a prescribed shear-free outlet.
This Ritz projection omits transverse strains, through-flow and mode conversion.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import cumulative_trapezoid

from mouth_curved_phase import quadratic
from mouth_flared_receiver import FlaredDomain, ROUNDING

ARCHIVE = Path('research/pde_ledger_v3/_measurements/mouth_free_order_depth.json')
FEED = .10
PACKET_LO, PACKET_HI = .16, .44
TIMES = (.5, 1., 2., 4., 8., 16.)
CUTS = (1., 1.65)
CFL = .72


def emit(**data):
    print(json.dumps(data, allow_nan=False), flush=True)


def verified(path, digests):
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != digests[str(path)]:
        raise ValueError(f'Input digest mismatch for {path}')
    return actual


def sections(dom, field, w):
    """Area average of the archived Q2 chi on each physical cross section."""
    a = dom.radius(w)[0]
    lip = w < ROUNDING
    a[lip] = dom.rc - np.sqrt(ROUNDING**2 - (w[lip] - ROUNDING)**2)
    avg = np.empty_like(w)
    # Rounded-end charts are curvilinear. Composite quadrature in r avoids
    # treating an entire section as one interpolation cell.
    x, g = leggauss(6)
    eta = ((np.arange(24)[:, None] + (x+1)/2)/24).ravel()
    weights = np.tile(g/48, 24) * 3 * eta**2
    np.testing.assert_allclose(weights.sum(), 1., rtol=0., atol=1.e-14)
    indices = np.flatnonzero(lip)
    for i in range(0, len(indices), 100):
        rows = indices[i:i+100]
        r = a[rows, None] * eta
        values = dom.read(field, r, np.broadcast_to(w[rows, None], r.shape))[0]
        avg[rows] = values @ weights
    # In the flared chart, integral r² chi dr is a fourth-degree polynomial
    # in v in each radial element. Gauss quadrature is exact there.
    pg, vg, dofs = dom.blocks[2]
    hv = np.diff(vg)/2
    v = vg[:-1, None] + hv[:, None]*(1+x)
    radial_weights = np.einsum('eq,qj->ej', 3*(1-v)**2*hv[:, None]*g, quadratic(x)[0])
    planes = np.einsum('abij,bj->ai', field[dofs].reshape(len(pg)-1, len(vg)-1, 3, 3), radial_weights)
    i = np.clip(np.searchsorted(pg, w[~lip])-1, 0, len(pg)-2)
    p = 2*(w[~lip]-pg[i])/np.diff(pg)[i]-1
    avg[~lip] = np.sum(quadratic(p)[0]*planes[i], axis=1)
    if np.any(avg <= 0) or np.any(avg > 1 + 1.e-10):
        raise ValueError('Area projection leaves the positive order interval')
    return (4*np.pi/3)*a**3, avg


def packet(w):
    where = np.clip((w-PACKET_LO)/(PACKET_HI-PACKET_LO), 0, 1)
    return np.where((w > PACKET_LO) & (w < PACKET_HI), np.sin(np.pi*where)**4, 0.)


def propagate(w, area, order, m, spacing, outlet_load, common_feed_impedance):
    # ds/dw=sqrt(M/K) gives unit characteristics and impedance sqrt(K M).
    # Truncate only outside the finite *discrete* dependence cone of the
    # source. A free edge there cannot affect any physical readout.
    speed_inverse = order**((m-1)/2)
    s_original = cumulative_trapezoid(speed_inverse, w, initial=0.)
    impedance = area*order**((1+m)/2)
    packet_end = np.interp(PACKET_HI, w, s_original)
    max_s = min(s_original[-1], packet_end + TIMES[-1]/(CFL*.95) + 1.)
    n = int(np.ceil(max_s/spacing))
    s = np.linspace(0., max_s, n+1)
    ds = s[1]
    grid = np.interp(s, s_original, w)
    impedance_edge = np.exp(np.interp((s[:-1]+s[1:])/2, s_original, np.log(impedance)))
    mass_edge = impedance_edge*ds
    mass = np.r_[mass_edge[0], mass_edge[:-1]+mass_edge[1:], mass_edge[-1]]/2
    # Initial displacement is zero, with the same compact axial velocity
    # support and shape, normalized to one unit of total energy in each law.
    v = packet(grid)
    launch_scale = 1/np.sqrt(.5*np.dot(mass, v*v))
    v *= launch_scale
    d = np.zeros(n)   # displacement gradient in travel time s
    # Identify the first recorded window exactly. Staggered energy lives at
    # half steps. With zero initial strain and no pulse at either end, both
    # velocities bracketing t=0 are the identical specified velocity.
    substeps = int(np.ceil(.5/(CFL*ds)))
    dt = .5/substeps
    if max_s < s_original[-1] and max_s - packet_end <= TIMES[-1]*ds/dt:
        raise ArithmeticError('Truncated grid entered the discrete influence cone')
    w_weights = [np.clip((grid[1:]-cut)/np.diff(grid), 0, 1) for cut in CUTS]
    targets = {round(t/dt):t for t in TIMES}
    bath = impedance[-1] if outlet_load == 'matched' else 0.
    loss = np.array([common_feed_impedance, bath])
    indices = [0, -1]
    port_half = .5*dt*loss/mass[indices]
    factors = dt/mass
    traction = np.empty_like(d)
    dv = np.empty_like(v)
    peaks = np.zeros(len(CUTS))
    samples = []
    crossings = {}
    escaped = np.zeros(2)
    max_error = 0.
    sample_every = max(1, round(.025/dt))
    for step in range(1, max(targets)+1):
        # Advance strain from integer time to the next, using v at the
        # intervening half step. Save d_old*d_new only when reading energy.
        next_d = d + (dt/ds)*np.diff(v)
        sample = step in targets or step % sample_every == 0
        if sample:
            energy_edge = mass_edge*(.25*(v[:-1]**2+v[1:]**2) + .5*d*next_d)
            remaining = energy_edge.sum()
            beyond = [np.dot(q, energy_edge) for q in w_weights]
            peaks = np.maximum(peaks, beyond)
            defect = float(remaining + sum(escaped) - 1.)
            max_error = max(abs(defect), max_error)
            if step in targets:
                samples.append(dict(t_half=(step-.5)*dt, remaining=remaining,
                                    escaped_feed=escaped[0], escaped_bath=escaped[1],
                                    beyond_w1=beyond[0], beyond_w165=beyond[1], balance=defect))
        d = next_d
        np.multiply(impedance_edge, d, out=traction)
        np.subtract(traction[1:], traction[:-1], out=dv[1:-1])
        dv[0], dv[-1] = traction[0], -traction[-1]
        old_port = v[indices]
        v += factors*dv
        v[indices] = (v[indices] - port_half*old_port)/(1+port_half)
        escaped += dt*loss*((v[indices]+old_port)/2)**2
        for level in (.5, .75, .9):
            if level not in crossings and escaped[0] >= level:
                crossings[level] = step*dt
    return dict(order_inertia_exponent=m, ds=ds, dt=dt, initial_amplitude=launch_scale,
                feed_impedance=loss[0], outlet_impedance=loss[1], s_bath=s_original[-1],
                last_simulated_w=grid[-1], peak_sampled_beyond=peaks.tolist(),
                feed_fraction_arrival=crossings, maximum_sample_balance_error=max_error,
                samples=samples)


def main(args):
    archive = json.loads(ARCHIVE.read_text())
    # Same discretization and launch; in a homogeneous line half the pulse
    # travels to each matched end. A free end sends its half back to the feed.
    w = np.linspace(FEED, 2., 1901)
    for load in ('matched', 'free'):
        res = propagate(w, np.ones_like(w), np.ones_like(w), 0, min(args.spacing), load, 1.)
        emit(test='uniform', load=load, **res)
        last = res['samples'][-1]
        expected = .5 if load == 'matched' else 1.
        if abs(last['escaped_feed']-expected) > 5.e-4 or res['maximum_sample_balance_error'] > 1.e-9:
            raise ArithmeticError('Uniform guide flux calibration failed')
    for profile in args.profiles:
        source = Path(f'_scratch/s11c/mouth_free_order_{profile}_n34_d04')
        header = json.loads((source/'stdout').read_text().splitlines()[0])
        carrier = Path(header['source'])
        verified(carrier, archive['input_carriers_sha256'])
        cb = json.loads((carrier.parent/'stdout').read_text().splitlines()[0])
        dom = FlaredDomain(header['n'], profile, header['roof'], cb['feed_far'])
        with np.load(carrier) as data:
            np.testing.assert_array_equal(data['points'], dom.points)
        for seed in (.44, 1.55):
            path = source/f'order_seed{seed}.npz'
            digest = verified(path, archive['diffuse_orders_sha256'])
            with np.load(path) as data:
                chi = data['chi']
            # Dense axial readout supports the travel-time map and coefficient
            # interpolation, especially in the slowly propagating deep tail.
            w = np.linspace(FEED, max(args.outlets), round((max(args.outlets)-FEED)*4000)+1)
            a, avg = sections(dom, chi, w)
            for bath_w in args.outlets:
                end = np.searchsorted(w, bath_w-1.e-12) + 1
                np.testing.assert_allclose(w[end-1], bath_w, rtol=0, atol=1.e-10)
                for m in (0, 1):
                    for spacing in args.spacing:
                        result = propagate(w[:end], a[:end], avg[:end], m, spacing, 'free', a[0])
                        emit(test='pulse', profile=profile, seed=seed, order_file=str(path), sha256=digest,
                             requested_spacing=spacing, bath_w=bath_w, mean_chi_at_bath=avg[end-1], **result)
                        if result['maximum_sample_balance_error'] > 1.e-9:
                            raise ArithmeticError('Discrete energy ledger failed')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--profiles', nargs='+', choices=('cubic','delayed'), default=['cubic','delayed'])
    p.add_argument('--outlets', type=float, nargs='+', default=[2., 2.2])
    p.add_argument('--spacing', type=float, nargs='+', default=[.004, .002])
    main(p.parse_args())
