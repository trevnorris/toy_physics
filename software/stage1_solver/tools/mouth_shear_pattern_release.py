"""Release one waist-clamped toroidal cavity mode into the open feeding slab.

The harmonic is normalized to the angular norm of the earlier ell=1 torsion.
Order is frozen, both shear stiffness and transverse inertia are proportional
to it, and no flow advection, driving or moving wall is included.
"""
import argparse
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.sparse.linalg import eigsh

from mouth_shear_cut import impedance, truncate
from mouth_shear_ringdown import assemble, digest, emit, energies, linear, prepared, product, quadratic, tensor, EXIT
from mouth_shear_ringdown_burst import momentum_and_mass


def angular_stiffness(dom, blocks, chi, size):
    """Integral of (8 pi/3) chi r^2 N_i N_j; ell=1 is already in assemble()."""
    roots, weights = leggauss(3)
    xi, eta = (x.ravel() for x in np.meshgrid(roots, roots, indexing='ij'))
    shape = product(linear(xi), linear(eta))
    parts = [[] for _ in range(6)]
    for ib, (pg, vg, cells) in enumerate(blocks):
        po, vo, old = dom.order_blocks[ib]
        rows, cols = (x.ravel() for x in np.indices((len(pg)-1, len(vg)-1)))
        for start in range(0, len(rows), 1024):
            batch = slice(start, start+1024)
            i, j, ids = rows[batch], cols[batch], cells[batch]
            hp, hv = np.diff(pg)[i, None]/2, np.diff(vg)[j, None]/2
            p, v = pg[i, None]+hp*(1+xi), vg[j, None]+hv*(1+eta)
            oi = np.searchsorted(po, (pg[i]+pg[i+1])/2)-1
            oj = np.searchsorted(vo, (vg[j]+vg[j+1])/2)-1
            ax = quadratic(2*(p-po[oi, None])/np.diff(po)[oi, None]-1)[0]
            ay = quadratic(2*(v-vo[oj, None])/np.diff(vo)[oj, None]-1)[0]
            order = np.clip(np.einsum('eqi,ei->eq', tensor(ax, ay), chi[old[oi*(len(vo)-1)+oj]]), 0, 1)
            r, w, rp, rv, wp, wv = dom.chart(ib, p, v)
            wt = (8*np.pi/3)*r**2*abs(rp*wv-rv*wp)*hp*hv*order*np.outer(weights, weights).ravel()
            groups = np.select([w < 0, w >= EXIT, order >= .9, order >= .5, order >= .1], [0, 5, 1, 2, 3], default=4)
            for region in np.unique(groups):
                take = (groups == region).any(axis=1)
                local = np.einsum('eq,qi,qj->eij', wt[take]*(groups[take] == region), shape, shape)
                parts[region].append((ids[take], local))
    out = []
    for pieces in parts:
        ids = np.concatenate([p[0] for p in pieces])
        val = np.concatenate([p[1] for p in pieces])
        mat = sparse.coo_matrix((val.ravel(), (np.repeat(ids, 4, axis=1).ravel(),
                                np.tile(ids, (1, 4)).ravel())), shape=(size, size)).tocsr()
        mat.eliminate_zeros()
        out.append(mat)
    return out


def run(a):
    dom, chi, intake, source = prepared(a)
    dom, chi, cut = truncate(dom, chi, a.cut)
    pts, blocks, _, parts, feed_z, _ = assemble(dom, chi, a.sub, 1, a.feed_spacing)
    mass_parts, _ = momentum_and_mass(dom, blocks, chi, len(pts), 0.)
    mass = mass_parts.sum(axis=0)
    if a.degree > 1:
        extra = angular_stiffness(dom, blocks, chi, len(pts))
        factor = a.degree*(a.degree+1)-2
        parts = [p+factor*b for p, b in zip(parts, extra, strict=True)]
    cut_z, _, _ = impedance(dom, blocks, pts, chi, cut, 0.)
    if a.termination == 'free':
        cut_z.fill(0.)
    active = (mass > 0) & ((pts[:, 0] > 0) | (a.degree == 1))
    mass, mass_parts, pts = mass[active], mass_parts[:, active], pts[active]
    feed_z, cut_z = feed_z[active], cut_z[active]
    parts = [p[active][:, active] for p in parts]
    k = sum(parts).tocsr()
    # Free soft end and fixed waist for preparation. Remove the waist clamp
    # and use the requested soft-end condition for the entire subsequent run.
    pocket = pts[:, 1] > .5 + 1e-10
    kk = k[pocket][:, pocket]
    invsqrt = sparse.diags(1/np.sqrt(mass[pocket]))
    norm_k = (invsqrt@kk@invsqrt).tocsr()
    eig, basis = eigsh(norm_k, k=max(3, a.branch), sigma=0, which='LM',
                       v0=np.random.default_rng(1729).standard_normal(norm_k.shape[0]), tol=2e-12)
    pick = np.argsort(eig)[a.branch-1]
    b = basis[:, pick]
    residual = np.linalg.norm(norm_k@b-eig[pick]*b)/eig[pick]
    vel = np.zeros_like(mass)
    vel[pocket] = (invsqrt@b)*np.sqrt(2)
    u = np.zeros_like(mass)
    e0 = .5*(mass@(vel**2))
    if residual > 1e-8 or abs(e0-1) > 1e-10:
        raise ArithmeticError('Unverified preparation eigenmode or normalization')
    initial_ordered = .5*(mass_parts[1:4].sum(axis=0)@(vel**2))/e0
    limit = np.max(np.asarray(abs(k).sum(axis=1)).ravel()/mass)
    sample = .2
    stride = int(np.ceil(sample*np.sqrt(limit)/1.4))
    dt = sample/stride
    emit(test='preparation', degree=a.degree, branch=a.branch, eigen_frequency=float(np.sqrt(eig[pick])),
         relative_eigen_residual=residual, initial_ordered=initial_ordered, energy=e0, intake=intake,
         order_file=str(source), order_sha256=digest(source), mesh=a.mesh, sub=a.sub, cut=cut,
         termination=a.termination, outer=a.outer, spacing=a.feed_spacing, end=a.end,
         nodes=len(mass), dt=dt, radius=a.radius)
    damping = feed_z+cut_z
    inv = dt/(mass+.5*dt*damping)
    ratio = (mass-.5*dt*damping)/(mass+.5*dt*damping)
    edge, far = np.flatnonzero(cut_z), np.flatnonzero(feed_z)
    lost_feed = lost_cut = 0.
    trace = []
    for step in range(round(a.end/sample)*stride):
        vn = ratio*vel-inv*(k@u)
        un = u+dt*vn
        lost_feed += .25*dt*(feed_z[far]@(vn[far]+vel[far])**2)
        lost_cut += .25*dt*(cut_z[edge]@(vn[edge]+vel[edge])**2)
        if (step+1)%stride == 0:
            kin, strain = energies(mass_parts, parts, vn, u, un)
            en = (kin+strain)/e0
            trace.append([(step+.5)*dt, en[0]+lost_feed/e0, en[1:4].sum(), en[4], en[5], lost_cut/e0,
                          en.sum()+(lost_feed+lost_cut)/e0-1])
        u, vel = un, vn
    trace = np.array(trace)
    err = float(np.max(abs(trace[:, -1])))
    if err > 1e-8 or not np.isfinite(trace).all():
        raise ArithmeticError('Open-cavity energy ledger failed')
    out = a.outdir/'trace.npz'
    cols = ['time', 'feed_net', 'active_throat', 'weak_throat', 'receiver', 'cut_export', 'closure']
    np.savez_compressed(out, columns=cols, trace=trace)
    windows = {}
    for lo, hi in ((4, 8), (14, 18), (24, 28), (44, 48), (94, 98), (144, 148), (194, 198)):
        if hi <= a.end:
            within = (trace[:, 0] >= lo) & (trace[:, 0] < hi)
            windows[f'{lo}-{hi}'] = dict(zip(cols[1:6], trace[within, 1:6].mean(axis=0), strict=True))
    emit(test='release', windows=windows, max_energy_error=err, trace_sha256=digest(out))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.set_defaults(radius=1., state=None)
    p.add_argument('--degree', type=int, required=True)
    p.add_argument('--branch', type=int, default=1)
    p.add_argument('--mesh', type=int, default=34)
    p.add_argument('--sub', type=int, default=1)
    p.add_argument('--cut', type=float, default=1.85)
    p.add_argument('--termination', choices=('free', 'matched'), default='matched')
    p.add_argument('--outer', type=float, default=28.)
    p.add_argument('--feed-spacing', type=float, default=.08)
    p.add_argument('--end', type=float, default=48.)
    p.add_argument('--outdir', type=Path, required=True)
    a = p.parse_args()
    if min(a.degree, a.branch, a.sub) < 1 or a.outer <= 9 or a.end <= 0:
        p.error('Positive harmonic, branch, refinement and duration and an extended feed are required.')
    run(a)
