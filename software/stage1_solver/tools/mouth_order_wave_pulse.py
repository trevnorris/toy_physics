"""Torsional pulse on archive-verified passive order fields in one flared 4D mouth.

X = f(r,w) (e_z cross x), with x in R^3.  The angularly integrated mass is
(8 pi/3) rho r^4.  The three angularly integrated curl squares are proportional
to r^2 [2 f^2 + (2 f+r f_r)^2 + r^2 f_w^2].  Weak variation sets the mechanical
traction to zero except for the distant feed's normal-impedance termination.
The order and carrier are held fixed; advection of the wave is omitted.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse

from mouth_curved_phase import divided, quadratic, tensor
from mouth_flared_receiver import EXIT, FlaredDomain


ROOT = Path('_scratch/s11c')
ARCHIVE = Path('research/pde_ledger_v3/_measurements/mouth_free_order_depth.json')
REGIONS = ('feed', 'ordered', 'upper_layer', 'lower_layer', 'tail', 'receiver')
AMPLITUDE = 1e-4
SAMPLE = .05


def emit(**fields):
    print(json.dumps(fields, allow_nan=False), flush=True)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def split(edges, subdivisions):
    t = np.arange(subdivisions) / subdivisions
    return np.r_[(edges[:-1, None] + np.diff(edges)[:, None]*t).ravel(), edges[-1]]


def mesh_geometry(dom, grids, subdivisions, quadratic_nodes, axial_counts=None, feed_spacing=None):
    blocks, chart_points = [], []
    count = 0
    for ib, (pg, vg) in enumerate(grids):
        counts = None
        if ib == 1 and feed_spacing:
            counts = np.maximum(subdivisions, np.ceil(np.diff(pg)/feed_spacing).astype(int))
        if ib == 2 and axial_counts is not None:
            counts = axial_counts
        if counts is not None:
            pg = np.r_[np.concatenate([np.linspace(a, b, n, endpoint=False)
                                       for a, b, n in zip(pg[:-1], pg[1:], counts, strict=True)]),
                       pg[-1]]
        else:
            pg = split(pg, subdivisions)
        vg = split(vg, subdivisions)
        pn, vn = (divided(pg), divided(vg)) if quadratic_nodes else (pg, vg)
        pp, vv = np.meshgrid(pn, vn, indexing='ij')
        r, z, *_ = dom.chart(ib, pp, vv)
        chart_points.append(np.c_[r.ravel(), z.ravel()])
        ids = np.arange(count, count+pp.size).reshape(pp.shape)
        count += pp.size
        i, j = (a.ravel() for a in np.indices((len(pg)-1, len(vg)-1)))
        step = 2 if quadratic_nodes else 1
        a = np.arange(step+1)
        dofs = ids[step*i[:, None, None]+a[None, :, None],
                   step*j[:, None, None]+a[None, None, :]].reshape(len(i), -1)
        blocks.append((pg, vg, dofs))
    points, mapping = np.unique(np.round(np.concatenate(chart_points), 11), axis=0,
                                return_inverse=True)
    return points, [(p, v, mapping[d]) for p, v, d in blocks]


class Geometry(FlaredDomain):
    """Parent builds the original charts. Only element topology is needed here."""
    def assemble(self, grids):
        self.grids = grids
        self.points, self.order_blocks = mesh_geometry(self, grids, 1, True)


def product(a, b):
    return (a[..., :, None]*b[..., None, :]).reshape(*a.shape[:-1], 4)


def linear(x):
    return np.stack((1-x, 1+x), axis=-1)/2


def assemble(dom, order, subdivisions, tail_refinement, feed_spacing):
    # Resolve wavelength shortening in the constant-inertia tail. Use the local
    # Q2 nodal envelope to select the number of axial subdivisions within each
    # original element row. The ceiling controls refinement, never stiffness.
    po, vo, archived = dom.order_blocks[2]
    envelope = order[archived].reshape(len(po)-1, len(vo)-1, 9).max(axis=(1, 2))
    speed = np.sqrt(np.clip(envelope, 0., 1.))
    boost = np.minimum(tail_refinement, 1/np.maximum(speed, 1e-30))
    # The ultra-soft far end lies beyond finite-time reach; keep several
    # well-resolved rows beyond the late packet before returning to the base grid.
    boost = np.where(speed < 4e-4, 1., boost)
    counts = np.ceil(subdivisions*boost).astype(int)
    points, blocks = mesh_geometry(dom, dom.grids, subdivisions, False, counts, feed_spacing)
    count = len(points)
    roots, weights = leggauss(3)
    xi, eta = (a.ravel() for a in np.meshgrid(roots, roots, indexing='ij'))
    shape = product(linear(xi), linear(eta))
    deriv = np.tile([-.5, .5], (len(xi), 1))
    shape_p, shape_v = product(deriv, linear(eta)), product(linear(xi), deriv)
    mass = np.zeros((len(REGIONS), count))
    stiffness_parts = [[] for _ in REGIONS]
    chi_min, chi_max = 1., 0.
    for ib, ((pg, vg, cells), (po, vo, archived)) in enumerate(
            zip(blocks, dom.order_blocks, strict=True)):
        ii, jj = (a.ravel() for a in np.indices((len(pg)-1, len(vg)-1)))
        for start in range(0, len(ii), 1536):
            chunk = slice(start, start+1536)
            i, j, dofs = ii[chunk], jj[chunk], cells[chunk]
            hp, hv = np.diff(pg)[i, None]/2, np.diff(vg)[j, None]/2
            p, v = pg[i, None]+hp*(1+xi), vg[j, None]+hv*(1+eta)
            oi = np.searchsorted(po, .5*(pg[i]+pg[i+1]))-1
            oj = j//subdivisions
            ax, _ = quadratic(2*(p-po[oi, None])/np.diff(po)[oi, None]-1)
            ay, _ = quadratic(2*(v-vo[oj, None])/np.diff(vo)[oj, None]-1)
            chi = np.einsum('eqi,ei->eq', tensor(ax, ay),
                            order[archived[oi*(len(vo)-1)+oj]])
            chi_min, chi_max = min(chi_min, chi.min()), max(chi_max, chi.max())
            chi = np.clip(chi, 0., 1.)
            r, w, rp, rv, wp, wv = dom.chart(ib, p, v)
            det = rp*wv-rv*wp
            gp, gv = shape_p/hp[:, :, None], shape_v/hv[:, :, None]
            gr = (wv[:, :, None]*gp-wp[:, :, None]*gv)/det[:, :, None]
            gz = (rp[:, :, None]*gv-rv[:, :, None]*gp)/det[:, :, None]
            wt = (8*np.pi/3)*r*r*abs(det)*hp*hv*np.outer(weights, weights).ravel()
            curl_r, curl_z = 2*shape+r[:, :, None]*gr, r[:, :, None]*gz
            region = np.select((w < 0., w >= EXIT, chi >= .9, chi >= .5, chi >= .1),
                               (0, 5, 1, 2, 3), default=4)
            for group in np.unique(region):
                keep = np.any(region == group, axis=1)
                ids = dofs[keep]
                part = wt[keep]*(region[keep] == group)
                # Positive row-sum lumping of the Q1 kinetic matrix. Retain
                # local mass partitions so the discrete region energies add.
                local_m = (part*r[keep]**2) @ shape
                np.add.at(mass[group], ids.ravel(), local_m.ravel())
                stiffness_w = part*chi[keep]
                local_k = 2*np.einsum('eq,qi,qj->eij', stiffness_w, shape, shape)
                for term in (curl_r[keep], curl_z[keep]):
                    local_k += np.einsum('eq,eqi,eqj->eij', stiffness_w, term, term)
                stiffness_parts[group].append((ids, local_k))
    stiff = []
    for parts in stiffness_parts:
        ids = np.concatenate([p[0] for p in parts])
        local = np.concatenate([p[1] for p in parts])
        sm = sparse.coo_matrix((local.ravel(),
                                (np.repeat(ids, 4, axis=1).ravel(),
                                 np.tile(ids, (1, 4)).ravel())), shape=(count, count)).tocsr()
        sm.eliminate_zeros()
        stiff.append(sm)
    # rho=mu=1 at the remote ordered feed. The saved order is fixed to zero
    # on all distant receiving surfaces, so their normal shear impedance is zero.
    b = np.zeros(count)
    remote = np.flatnonzero((points[:, 0] == dom.far) & (points[:, 1] <= 0))
    remote = remote[np.argsort(points[remote, 1])]
    edge_m = (8*np.pi/3)*dom.far**4*np.diff(points[remote, 1])/2
    b[remote[:-1]] += edge_m
    b[remote[1:]] += edge_m
    if chi_min < -1e-9 or chi_max > 1+1e-9 or mass.sum(axis=0).min() <= 0:
        raise ArithmeticError('Archive order range or wave mass invalid')
    return points, blocks, mass, stiff, b, [float(chi_min), float(chi_max)]


def pulse(points, center, halfwidth):
    r, w = points.T
    return AMPLITUDE*(np.maximum(0., 1-((w-center)/halfwidth)**2)*
                      np.maximum(0., 1-(r/.48)**2))**4


def probe(dom, blocks, r, z):
    if z < .5:
        ib = 0
        p = np.arctan2(z-dom.rounding, dom.rc-r)
        v = (np.hypot(z-dom.rounding, dom.rc-r)-dom.rounding)/(dom.reach(p)[0]-dom.rounding)
    else:
        ib, p, v = 2, z, 1-r/dom.radius(z)[0]
    pg, vg, dofs = blocks[ib]
    i, j = np.searchsorted(pg, p)-1, np.searchsorted(vg, v)-1
    ap = linear(2*(p-pg[i])/np.diff(pg)[i]-1)
    av = linear(2*(v-vg[j])/np.diff(vg)[j]-1)
    return product(ap, av), dofs[i*(len(vg)-1)+j]


def energies(mass, stiff, v, previous, current):
    kinetic = .5*(mass @ (v*v))
    strain = np.array([.5*(previous @ (k @ current)) for k in stiff])
    return kinetic, strain


def evolve(args, scheme, config, state, state_sha, excitation):
    points, blocks, mass, stiff, damp, chi_range = scheme
    total_m = mass.sum(axis=0)
    K = sum(stiff).tocsr()
    # Gershgorin bound on the generalized eigenvalues of M^{-1} K.
    bound = np.max(np.asarray(abs(K).sum(axis=1)).ravel()/total_m)
    steps_per_bin = int(np.ceil(SAMPLE*np.sqrt(bound)/(2*args.cfl)))
    dt = SAMPLE/steps_per_bin
    feed_nodes = np.flatnonzero(damp)
    if excitation == 'entrance_velocity':
        u = np.zeros(len(points))
        initial_v = pulse(points, .25, .20)
    else:
        u = pulse(points, .875, .225)
        initial_v = np.zeros(len(points))
    pv, sv = (probe(config, blocks, .27, z) for z in (.25, .875))
    phys_kinetic, phys_strain = energies(mass, stiff, initial_v, u, u)
    energy0 = float(phys_kinetic.sum()+phys_strain.sum())
    v = initial_v + dt/2*((K @ u)+damp*initial_v)/total_m
    # Starting state at t=-dt/2. The centered initial velocity is initial_v.
    a, s = energies(mass, stiff, v, u-dt*v, u)
    invariant0 = (a+s).sum()
    inv = dt/(total_m+.5*dt*damp)
    back = (total_m-.5*dt*damp)/(total_m+.5*dt*damp)
    escaped, maxerr = 0., 0.
    rows = []
    physical_initial = (phys_kinetic+phys_strain)/energy0
    stop = int(round(args.end/SAMPLE))
    if not np.isclose(stop*SAMPLE, args.end):
        raise ValueError('end must be a multiple of the sampling interval')

    def record(t, kin, strain, displacement, velocity):
        nonlocal maxerr
        regional = (kin+strain)/invariant0
        residual = regional.sum()+escaped/invariant0-1
        maxerr = max(maxerr, abs(residual))
        row = dict(t=t, **dict(zip(REGIONS, regional, strict=True)),
                   feed_absorbed=escaped/invariant0,
                   throat_kinetic=kin[1:5].sum()/invariant0,
                   throat_strain=strain[1:5].sum()/invariant0,
                   f_entrance=pv[0] @ displacement[pv[1]],
                   f_cavity=sv[0] @ displacement[sv[1]],
                   fdot_entrance=pv[0] @ velocity[pv[1]],
                   fdot_cavity=sv[0] @ velocity[sv[1]], balance_error=residual)
        rows.append(row)

    # All evolved samples are at half times. The initial row uses the
    # physical integer-time observables, with its own physical normalization.
    rows.append(dict(t=0., **dict(zip(REGIONS, physical_initial, strict=True)),
                     feed_absorbed=0., throat_kinetic=phys_kinetic[1:5].sum()/energy0,
                     throat_strain=phys_strain[1:5].sum()/energy0,
                     f_entrance=pv[0] @ u[pv[1]], f_cavity=sv[0] @ u[sv[1]],
                     fdot_entrance=pv[0] @ initial_v[pv[1]],
                     fdot_cavity=sv[0] @ initial_v[sv[1]], balance_error=0.))
    for step in range(1, stop*steps_per_bin+1):
        vf = v[feed_nodes]
        v = back*v - inv*(K @ u)
        vcenter = .5*(vf+v[feed_nodes])
        escaped += dt*(damp[feed_nodes] @ (vcenter*vcenter))
        next_u = u + dt*v
        if step % steps_per_bin == 0:
            kinetic, strain = energies(mass, stiff, v, u, next_u)
            record((step-.5)*dt, kinetic, strain, .5*(u+next_u), v)
        u = next_u
    out = args.outdir/f'{state}_{excitation}.csv'
    with out.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    # Derivatives of cumulative deliveries give signed window-averaged powers.
    # Positive on each outward leg. The initialized packet is in the throat.
    milestones = [0, 1, 2, 4, 6, 8, 10, 12, 16, args.end]
    snapshots = [rows[round(t/SAMPLE)] for t in sorted(set(milestones)) if t <= args.end]
    peaks = {}
    for name, region_names in [('feed', ('feed', 'feed_absorbed')), ('receiver', ('receiver',))]:
        fractions = np.array([sum(row[key] for key in region_names) for row in rows])
        powers = np.diff(fractions)/np.diff([row['t'] for row in rows])
        peaks[name] = dict(peak_net_delivery_per_time=powers.max(),
                           max_energy_delivered=fractions.max(),
                           min_net_delivery_per_time=powers.min())
    signals = {}
    for key in ('f_entrance', 'f_cavity'):
        times = np.array([r['t'] for r in rows])
        readings = np.array([r[key] for r in rows])
        negative = np.flatnonzero(readings < 0.)
        late = times >= args.end-2
        signals[key] = dict(first_negative_sample=times[negative[0]] if negative.size else None,
                            initial=readings[0], late_min=readings[late].min(),
                            late_max=readings[late].max())
    emit(test='torsional_pulse', profile=args.profile, n=args.n, subdivisions=args.subdivisions,
         state=state, order_sha256=state_sha, pulse=excitation, amplitude=AMPLITUDE,
         tail_refinement=args.tail_refinement, feed_spacing=args.feed_spacing,
         wave_dofs=len(points), dt=dt, cfl=args.cfl, rho=1., chi_quadrature_range=chi_range,
         initial_physical_energy=energy0, initial_kinetic_energy=float(phys_kinetic.sum()),
         initial_strain_energy=float(phys_strain.sum()), initial_invariant=invariant0,
         maximum_abs_relative_balance_residual=maxerr, history=str(out),
         history_sha256=digest(out), outward_flux=peaks, signals=signals, readings=snapshots)
    if not np.isfinite(maxerr) or maxerr > 1e-8:
        raise ArithmeticError('Discrete energy balance failed')


def exit_transit_bound(dom, order):
    # Convert each tensor-Q2 Lagrange representation on full cross-sections to
    # tensor Bernstein control values. Nonnegative basis and partition of unity
    # bound chi on the entire cell, including between Gauss nodes. Every path
    # to the receiver must cross this w interval with |dw/dt| <= sqrt(chi).
    pg, vg, dofs = dom.order_blocks[2]
    rows = np.flatnonzero((pg[:-1] >= 2.) & (pg[1:] <= 2.5))
    local = order[dofs].reshape(len(pg)-1, len(vg)-1, 3, 3)[rows]
    transform = np.array([[1., 0., 0.], [-.5, 2., -.5], [0., 0., 1.]])
    upper = float(np.clip(np.max(transform @ local @ transform.T), 0., 1.))
    length = pg[rows[-1]+1]-pg[rows[0]]
    return dict(w_min=pg[rows[0]], w_max=pg[rows[-1]+1], order_upper=upper,
                speed_upper=np.sqrt(upper), minimum_transit_time=length/np.sqrt(upper))


def run(args):
    archive = json.loads(ARCHIVE.read_text())
    source = ROOT/f'mouth_flared_{args.profile}_n{args.n}'/'field.npz'
    source_hash = digest(source)
    if archive['input_carriers_sha256'][str(source)] != source_hash:
        raise ValueError('Carrier digest mismatch')
    cfg = json.loads((source.parent/'stdout').read_text().splitlines()[0])
    dom = Geometry(cfg['n'], cfg['profile'], cfg['roof'], cfg['feed_far'])
    with np.load(source) as data:
        np.testing.assert_array_equal(dom.points, data['points'])
    if abs(cfg['bulk']-.25) > 1e-8:
        raise ValueError('Carrier intake mismatch')
    for state, seed in (('near', '0.44'), ('deep', '1.55')):
        folder = ROOT/f'mouth_free_order_{args.profile}_n{args.n}_d04'
        if args.profile == 'cubic' and args.n == 48 and state == 'deep':
            folder = folder.with_name(folder.name+'_outer_residual')
        path = folder/f'order_seed{seed}.npz'
        state_sha = digest(path)
        if archive['diffuse_orders_sha256'][str(path)] != state_sha:
            raise ValueError('Stationary order digest mismatch')
        with np.load(path) as data:
            order = data['chi']
        np.testing.assert_array_equal(order[dom.feed], 1.)
        np.testing.assert_array_equal(order[dom.bulk], 0.)
        emit(test='exit_transit_bound', n=args.n, profile=args.profile, state=state,
             **exit_transit_bound(dom, order))
        scheme = assemble(dom, order, args.subdivisions, args.tail_refinement, args.feed_spacing)
        for excitation in ('entrance_velocity', 'cavity_displacement'):
            evolve(args, scheme, dom, state, state_sha, excitation)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--n', type=int, choices=(26, 34, 48), required=True)
    ap.add_argument('--profile', choices=('cubic', 'delayed'), required=True)
    ap.add_argument('--subdivisions', type=int, default=1)
    ap.add_argument('--tail-refinement', type=int, default=1)
    ap.add_argument('--feed-spacing', type=float, help='Maximum radial cell width in remote feed chart')
    ap.add_argument('--cfl', type=float, default=.72)
    ap.add_argument('--end', type=float, default=12.)
    ap.add_argument('--outdir', type=Path, required=True)
    args = ap.parse_args()
    if min(args.subdivisions, args.tail_refinement) < 1 or not (0 < args.cfl < 1) or args.end < SAMPLE:
        ap.error('Positive subdivision count, 0<cfl<1 and end >= one sample required')
    if args.feed_spacing is not None and args.feed_spacing <= 0:
        ap.error('Feed cell width must be positive')
    run(args)
