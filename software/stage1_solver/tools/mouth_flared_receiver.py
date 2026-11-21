"""Stationary single-phase flow through a fixed rounded mouth, collar and flared conduit.

The third brane coordinate is represented by the exact S² integration measure
4*pi*r² dr dw. Both finite reservoirs belong to the nonlinear field solve.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse

from mouth_curved_phase import (Domain, K_FRONT, RATE, RIM_SLOPE, divided, front,
                                integrated_radial_current, quadratic, readings, tensor)
from mouth_curved_phase_stability import half_plane_bound


ROUNDING = .5
COLLAR_LENGTH = .7
BORE = 1.
EXIT_RADIUS = 1.5
FLARE_LENGTH = 1.85
FLARE_START = ROUNDING + COLLAR_LENGTH
EXIT = FLARE_START + FLARE_LENGTH


class FlaredDomain(Domain):
    """Three conforming parts: rounded feed, shaped conduit, receiving chamber."""

    def __init__(self, n, profile, roof, far=9.):
        self.far, self.rounding, self.slope = far, ROUNDING, RIM_SLOPE
        self.profile, self.roof = profile, roof
        angle = np.arctan(self.slope)
        self.rc = BORE + ROUNDING*np.sin(angle)
        self.lip = ROUNDING*(1 - np.cos(angle))
        self.neck = self.rc - ROUNDING
        self.bottom_turn = -np.arctan((1 + ROUNDING)/self.rc)
        self.top_turn = np.arctan((far + ROUNDING)/self.rc)
        breaks = [-np.pi/2, angle - np.pi/2, self.bottom_turn, 0.]
        angular = np.r_[np.concatenate([np.linspace(a, b, max(1, int(np.ceil((b-a)*n))),
                                                  endpoint=False)
                                       for a, b in zip(breaks[:-1], breaks[1:], strict=True)]), 0.]
        stretch = np.sinh(2.4*np.linspace(0, 1, n+1)) / np.sinh(2.4)
        radial = self.rc + (far-self.rc)*stretch
        axial = np.r_[np.linspace(ROUNDING, FLARE_START, max(2, int(np.ceil(n*COLLAR_LENGTH)))+1)[:-1],
                      np.linspace(FLARE_START, EXIT, max(2, int(np.ceil(n*FLARE_LENGTH)))+1)]
        # Refine near the common exit corner. The roof/depth control keeps the
        # complete conduit mesh and its entrance independent of bath placement.
        chamber = EXIT + (roof-EXIT)*stretch
        self.assemble([(angular, stretch), (radial, stretch), (axial, stretch),
                       (chamber, stretch), (chamber, stretch)])
        self.feed = (self.points[:, 0] == far) & (self.points[:, 1] <= 0)
        self.bulk = ((self.points[:, 0] == far) & (self.points[:, 1] >= round(EXIT, 11))) | (self.points[:, 1] == roof)
        self.free = ~(self.feed | self.bulk)
        assert np.any(self.feed) and np.any(self.bulk)

    def radius(self, z):
        """Common C² collar; cubic historic or delayed quartic flare at fixed ends."""
        t = np.clip((z-ROUNDING)/COLLAR_LENGTH, 0, 1)
        rise = BORE-self.neck
        # Quintic endpoint interpolant. At the waist a''=1/rounding
        # matches the original circular lip. At the bore a'=a''=0.
        bend = COLLAR_LENGTH**2/(2*ROUNDING)
        collar = self.neck + rise*t**3*(10-15*t+6*t*t) + bend*t*t*(1-t)**3
        dc = (rise*30*t*t*(1-t)**2 + bend*t*(1-t)**2*(2-5*t))/COLLAR_LENGTH
        u = np.clip((z-FLARE_START)/FLARE_LENGTH, 0, 1)
        if self.profile == 'cubic':
            flare, df = u*u*(3-2*u), 6*u*(1-u)
        else:
            flare, df = u**3*(4-3*u), 12*u*u*(1-u)
        return (np.where(z <= FLARE_START, collar, BORE+(EXIT_RADIUS-BORE)*flare),
                np.where(z <= FLARE_START, dc, (EXIT_RADIUS-BORE)*df/FLARE_LENGTH))

    def chart(self, block, p, v):
        if block < 2:
            return super().chart(block, p, v)
        one, zero = np.ones_like(p), np.zeros_like(p)
        if block == 4:
            width = self.far-EXIT_RADIUS
            return EXIT_RADIUS+v*width, p, zero, width*one, one, zero
        radius, slope = self.radius(p) if block == 2 else (EXIT_RADIUS*one, zero)
        return radius*(1-v), p, slope*(1-v), -radius, one, zero

    def assemble(self, grids):
        """Use the original tensor-Q2 basis and four-point Gauss rule on each chart."""
        nodes, cell_dofs, stiffness, weights = [], [], [], []
        root, gw = leggauss(4)
        xi, eta = [x.ravel() for x in np.meshgrid(root, root, indexing='ij')]
        xp, dx = quadratic(xi)
        yp, dy = quadratic(eta)
        basis = tensor(xp, yp)
        offset = 0
        self.blocks = {}
        for block, (pg, vg) in enumerate(grids):
            pp, vv = np.meshgrid(divided(pg), divided(vg), indexing='ij')
            r, z, *_ = self.chart(block, pp, vv)
            nodes.append(np.c_[r.ravel(), z.ravel()])
            ids = np.arange(pp.size).reshape(pp.shape)+offset
            offset += pp.size
            i, j = [a.ravel() for a in np.indices((len(pg)-1, len(vg)-1))]
            dofs = ids[2*i[:, None, None]+np.arange(3)[None, :, None],
                       2*j[:, None, None]+np.arange(3)[None, None, :]].reshape(-1, 9)
            self.blocks[block] = pg, vg, dofs
            cell_dofs.append(dofs)
            hp, hv = np.diff(pg)[i, None]/2, np.diff(vg)[j, None]/2
            r, _, rp, rv, zp, zv = self.chart(block, pg[i, None]+hp*(1+xi), vg[j, None]+hv*(1+eta))
            det = rp*zv-rv*zp
            if np.any(abs(det) < 1e-13) or np.any(np.sign(det) != np.sign(det.flat[0])):
                raise ValueError("Folded or singular integration chart")
            gp, gv = tensor(dx, yp)/hp[:, :, None], tensor(xp, dy)/hv[:, :, None]
            gr = (zv[:, :, None]*gp-zp[:, :, None]*gv)/det[:, :, None]
            gz = (rp[:, :, None]*gv-rv[:, :, None]*gp)/det[:, :, None]
            wt = 4*np.pi*r*r*abs(det)*hp*hv*np.outer(gw, gw).ravel()
            stiffness.append(np.einsum('eq,eqi,eqj->eij', wt, gr, gr) +
                             np.einsum('eq,eqi,eqj->eij', wt, gz, gz))
            weights.append(wt)
        self.points, remap = np.unique(np.round(np.concatenate(nodes), 11), axis=0, return_inverse=True)
        for block, (pg, vg, dofs) in self.blocks.items():
            self.blocks[block] = pg, vg, remap[dofs]
        dofs = remap[np.concatenate(cell_dofs)]
        count = len(self.points)
        self.stiff = sparse.coo_matrix((np.concatenate(stiffness).ravel(),
            (np.repeat(dofs, 9, axis=1).ravel(), np.tile(dofs, (1, 9)).ravel())),
            shape=(count, count)).tocsr()
        wt = np.concatenate(weights)
        self.weight = wt.ravel()
        shape = (*wt.shape, 9)
        qnodes = np.arange(wt.size).reshape(*wt.shape, 1)
        self.N = sparse.coo_matrix((np.broadcast_to(basis, shape).ravel(),
            (np.broadcast_to(qnodes, shape).ravel(), np.broadcast_to(dofs[:, None, :], shape).ravel())),
            shape=(wt.size, count)).tocsr()
        if np.max(abs(self.stiff @ np.ones(count))) > 1e-8:
            raise ValueError("Constant field failed the gradient check")
        self.check_seams()

    def check_seams(self):
        def indices(block):
            pg, vg, dofs = self.blocks[block]
            return dofs.reshape(len(pg)-1, len(vg)-1, 3, 3)
        a, b, c, d, e = [indices(i) for i in range(5)]
        for left, right in ((a[0, :, 0, :], b[0, :, 0, :]),
                            (a[-1, :, -1, :], c[0, :, 0, :]),
                            (c[-1, :, -1, :], d[0, :, 0, :]),
                            (d[:, 0, :, 0], e[:, 0, :, 0])):
            np.testing.assert_array_equal(left, right)

    def section(self, field, z):
        """Direct axial flux on a cross-section in the conduit using elementwise Gauss integration."""
        pg, vg, dofs = self.blocks[2]
        i = np.searchsorted(pg, z)-1
        hp = (pg[i+1]-pg[i])/2
        ap, dap = quadratic(np.array((z-pg[i])/hp-1))
        root, gw = leggauss(6)
        av, dav = quadratic(root)
        hv = np.diff(vg)/2
        v = vg[:-1, None]+hv[:, None]*(1+root)
        values = field[dofs[i*(len(vg)-1):(i+1)*(len(vg)-1)]].reshape(-1, 3, 3)
        f = np.einsum('i,qj,eij->eq', ap, av, values)
        fp = np.einsum('i,qj,eij->eq', dap, av, values)/hp
        fv = np.einsum('i,qj,eij->eq', ap, dav, values)/hv[:, None]
        radius, slope = self.radius(z)
        fz = fp+(1-v)*slope*fv/radius
        jw = (f.conj()*fz).imag
        flux = np.sum(4*np.pi*(radius*(1-v))**2*jw*radius*hv[:, None]*gw)
        return dict(w=z, radius=float(radius), integrated_jw=flux, min_sampled_jw=jw.min())


def emit(**fields):
    print(json.dumps(fields, allow_nan=False), flush=True)


def run(args):
    dom = FlaredDomain(args.n, args.profile, args.roof)
    field, state = dom.solve(RATE)
    certificate = half_plane_bound(dom, field, state['bias'])
    np.testing.assert_allclose(field[dom.feed], 1., rtol=0, atol=1e-14)
    np.testing.assert_allclose(field[dom.bulk], np.exp(1j*state['bias']), rtol=0, atol=1e-14)
    emit(test='carrier', n=args.n, profile=args.profile, roof=args.roof, feed_far=dom.far,
         neck=dom.neck, neck_w=ROUNDING, bore=BORE, flare_start=FLARE_START,
         exit_radius=EXIT_RADIUS, exit_w=EXIT, midpoint_radius=float(dom.radius(FLARE_START+FLARE_LENGTH/2)[0]),
         **state, **certificate)
    sections = [dom.section(field, z) for z in (.638, .993, 1.417, 2.083, 2.793, 2.97)]
    emit(test='sections', values=sections)
    hcurve = front(dom, field, K_FRONT, args.steps, args.panels)
    r = np.linspace(0, 1, 2001)
    h, hp = hcurve(r), hcurve(r, 1)
    jr, jw = dom.current(field, r, h)
    q = jw-hp*jr
    result = readings(dom, field, r, h, hp)
    sr = np.r_[np.linspace(0, .98, 50), .99, 1.]
    cyl = -4*np.pi*sr**2*integrated_radial_current(dom, field, sr, hcurve(sr), 2*args.panels)
    shell = np.diff(cyl)
    wall_curvature = (1+dom.slope**2)**1.5/ROUNDING
    emit(test='front', K=K_FRONT, steps=args.steps, panels=args.panels,
         positive_interior_count=int(np.sum(q[:-1] > 0)), interior_count=len(q)-1,
         min_interior_q=q[:-1].min(), positive_shell_count=int(np.sum(shell > 0)),
         shell_count=len(shell), min_shell=shell.min(), last_shell=shell[-1],
         cylinder_total=cyl[-1], grazing_qprime_from_slip=jr[-1]*(wall_curvature-2*dom.slope),
         **result)
    if args.reference:
        with np.load(args.reference) as reference:
            np.testing.assert_array_equal(r, reference['r'])
            q0, h0 = reference['q_0.4'], reference['h_0.4']
        emit(test='reference', reference=str(args.reference),
             sha256=hashlib.sha256(args.reference.read_bytes()).hexdigest(),
             delta_center_q=q[0]-q0[0], delta_q_at_099=q[1980]-q0[1980],
             max_abs_delta_q=np.max(abs(q-q0)), delta_center_depth=h[0]-h0[0],
             max_abs_delta_depth=np.max(abs(h-h0)))
    if (certificate['half_plane_lower'] <= 0 or h.max() >= ROUNDING or
            result['min_wall_clearance'] < -1e-10 or abs(result['total']-RATE) > .00125 or
            max(abs(s['integrated_jw']-RATE) for s in sections) > .0025):
        raise ArithmeticError('Carrier branch, front placement or integrated current failed verification')
    if args.output:
        np.savez_compressed(args.output, points=dom.points, psi=field, r=r, h=h,
                            q=q, hprime=hp, phase=np.angle(dom.read(field, r, h)[0]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, required=True)
    parser.add_argument('--profile', choices=('cubic', 'delayed'), required=True)
    parser.add_argument('--roof', type=float, default=10.)
    parser.add_argument('--steps', type=int, default=80)
    parser.add_argument('--panels', type=int, default=60)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--reference', type=Path)
    args = parser.parse_args()
    if args.n < 4 or args.roof <= EXIT or min(args.steps, args.panels) < 1:
        parser.error('Require n >= 4, roof above exit, and positive integration discretizations')
    run(args)
