"""Free stationary advected equal-well order on an archived single-phase flared mouth.

For u = grad(theta), exp(-theta/nu) is an integrating factor for the
Allen–Cahn equation. Minimize its weighted energy with reflecting natural
boundary and fixed distant phases on the original S²-weighted Q2 mesh.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.optimize import brentq
from scipy.sparse.linalg import eigsh, splu
from scipy.special import expit

from mouth_curved_phase import quadratic, tensor
from mouth_flared_receiver import EXIT, ROUNDING, FlaredDomain

NU = .4
SIGMA = 1 / 3
ARCHIVE = Path('research/pde_ledger_v3/_measurements/mouth_order_and_flare.json')


def emit(**values):
    print(json.dumps(values, allow_nan=False), flush=True)


def weighted_stiffness(dom, field):
    roots, _ = leggauss(4)
    xi, eta = (x.ravel() for x in np.meshgrid(roots, roots, indexing='ij'))
    ap, dp = quadratic(xi)
    av, dv = quadratic(eta)
    value, grad_p, grad_v = tensor(ap, av), tensor(dp, av), tensor(ap, dv)
    nodes, local_stiffness = [], []
    offset = 0
    phase_q = np.angle(dom.N @ field)
    if np.max(abs(phase_q)) >= np.pi / 2:
        raise ValueError('carrier phase leaves chosen principal-branch interval')
    wt = dom.weight * np.exp(-phase_q / NU)
    for block, (pg, vg, dofs) in dom.blocks.items():
        i, j = (x.ravel() for x in np.indices((len(pg)-1, len(vg)-1)))
        hp, hv = np.diff(pg)[i, None] / 2, np.diff(vg)[j, None] / 2
        _, _, rp, rv, zp, zv = dom.chart(block, pg[i, None]+hp*(1+xi), vg[j, None]+hv*(1+eta))
        det = rp*zv-rv*zp
        gp, gv = grad_p / hp[:, :, None], grad_v / hv[:, :, None]
        gr = (zv[:, :, None]*gp-zp[:, :, None]*gv)/det[:, :, None]
        gz = (rp[:, :, None]*gv-rv[:, :, None]*gp)/det[:, :, None]
        w = wt.reshape(-1, len(xi))[offset:offset+len(dofs)]
        local_stiffness.append(np.einsum('eq,eqi,eqj->eij', w, gr, gr) +
                               np.einsum('eq,eqi,eqj->eij', w, gz, gz))
        nodes.append(dofs)
        offset += len(dofs)
    dofs = np.concatenate(nodes)
    indices = np.repeat(dofs, 9, axis=1).ravel(), np.tile(dofs, (1, 9)).ravel()
    stiff = sparse.coo_matrix((np.concatenate(local_stiffness).ravel(), indices),
                              shape=dom.stiff.shape).tocsr()
    assert np.max(abs(stiff @ np.ones(len(field)))) < 1e-8
    return stiff, wt


def positive_factor(matrix):
    # With identical permutations and no numerical row exchanges, these are
    # the diagonal pivots of a symmetric LDLᵀ elimination (up to roundoff).
    fac = splu(matrix.tocsc(), permc_spec='MMD_AT_PLUS_A', diag_pivot_thresh=0,
               options={'SymmetricMode': True, 'Equil': False})
    if np.array_equal(fac.perm_r, fac.perm_c) and np.min(fac.U.diagonal()) > 0:
        return fac
    return None


class Order:
    def __init__(self, dom, field, delta):
        self.dom, self.delta = dom, delta
        self.stiff, self.wt = weighted_stiffness(dom, field)
        self.Nf = dom.N[:, dom.free]
        self.K = NU * self.stiff[dom.free][:, dom.free]
        self.scale = NU / delta**2
        self.mass = self.mass_weight(self.wt)

    def mass_weight(self, weights):
        return self.Nf.T @ self.Nf.multiply(weights[:, None])

    def energy(self, x):
        c = self.dom.N @ x
        return .5 * (NU * (x @ (self.stiff @ x)) + self.scale * (self.wt @ (c*(1-c))**2))

    def measure(self, x):
        c = self.dom.N @ x
        resid = NU*(self.stiff @ x)[self.dom.free] + self.Nf.T @ (
            self.scale * self.wt * c*(1-c)*(1-2*c))
        jac = self.K + self.mass_weight(self.scale * self.wt * (1-6*c*(1-c)))
        return resid, jac

    def minimize(self, seed):
        dom = self.dom
        x = expit((seed-dom.points[:, 1])/self.delta)
        x[dom.feed], x[dom.bulk] = 1., 0.
        energy = self.energy(x)
        shift = 0.
        for iteration in range(1, 501):
            resid, jac = self.measure(x)
            factor = positive_factor(jac)
            if np.max(abs(resid)) < 2e-11 and factor is not None:
                # Positive pivots first establish all eigenvalues > 0. Shift
                # invert then returns the smallest relaxation rate of this mesh.
                eig = eigsh(jac, k=1, M=self.mass, sigma=0., which='LM',
                            v0=np.ones(len(resid)), tol=1e-7, return_eigenvectors=False)[0]
                return x, dict(iterations=iteration, energy=energy,
                    max_weak_residual=np.max(abs(resid)), min_LDL_pivot=factor.U.diagonal().min(),
                    minimum_meridional_decay_rate=eig, nodal_min=x.min(), nodal_max=x.max())
            if factor is None:
                shift = max(1., shift/2)
                while factor is None:
                    factor = positive_factor(jac + shift*self.mass)
                    if factor is None:
                        shift *= 2
            else:
                shift = 0.
            direction = factor.solve(-resid)
            slope = resid @ direction
            if slope >= 0:
                raise ArithmeticError('positive factor yielded non-descent')
            step = min(1., .22 / np.max(abs(direction)))
            while True:
                candidate = x.copy()
                candidate[dom.free] += step*direction
                e_next = self.energy(candidate)
                if e_next <= energy + 1e-4*step*slope + 3e-13*(1+abs(energy)):
                    break
                step *= .5
                if step < 1e-10:
                    raise ArithmeticError('line search exhausted floating-point descent')
            x, energy = candidate, e_next
        raise ArithmeticError('stationary solver did not settle after 500 descent iterations')


def section_geometry(dom, z):
    if z >= ROUNDING:
        return dom.radius(z)
    width = np.sqrt(ROUNDING**2 - (z-ROUNDING)**2)
    return dom.rc-width, (z-ROUNDING)/width


def read(dom, field, r, z):
    if z < ROUNDING:
        return dom.read(field, np.asarray(r), np.asarray(z))
    radius, slope = dom.radius(z)
    p, v = z, 1-r/radius
    pg, vg, dofs = dom.blocks[2]
    i = np.clip(np.searchsorted(pg, p)-1, 0, len(pg)-2)
    j = np.clip(np.searchsorted(vg, v)-1, 0, len(vg)-2)
    hp, hv = np.diff(pg)[i]/2, np.diff(vg)[j]/2
    ap, dp = quadratic((p-pg[i])/hp-1)
    av, dv = quadratic((v-vg[j])/hv-1)
    values = field[dofs[i*(len(vg)-1)+j]]
    val = tensor(ap,av) @ values
    fp, fv = tensor(dp,av) @ values / hp, tensor(ap,dv) @ values / hv
    return val, -fv/radius, fp + fv*(1-v)*slope/radius


def diagnostics(dom, carrier, x, delta):
    mobility = NU/(6*SIGMA*delta)
    def sample(eta, level, brackets):
        def err(z):
            radius, _ = section_geometry(dom, z)
            return read(dom, x, eta*radius, z)[0] - level
        z = brentq(err, *brackets, xtol=2e-13)
        radius, slope = section_geometry(dom, z)
        r = eta*radius
        _, cr, cz = read(dom, x, r, z)
        psi, pr, pz = read(dom, carrier, r, z)
        n = abs(psi)**2
        ur, uz = (psi.conj()*pr).imag/n, (psi.conj()*pz).imag/n
        g = np.hypot(cr, cz)
        un = -(ur*cr+uz*cz)/g
        res = dict(eta=eta, r=r, w=z, n=n, grad_order=g, u_normal=un,
                   j_normal=n*un, mu_per_n=un*g/(mobility*n))
        if eta == 1.:
            normal_dot = -(cr-slope*cz)/(g*np.hypot(1, slope))
            res['wall_normal_dot'] = normal_dot
            res['angle_between_normals_deg'] = np.degrees(np.arccos(np.clip(normal_dot, -1, 1)))
        return res
    # Detect sign changes across the whole rounded-to-flared conduit section.
    # One threshold per ray is needed to interpret the closed-loop legs.
    zscan = np.linspace(.002, EXIT-.002, 340)
    def bracket(eta):
        vals = np.array([read(dom, x, eta*section_geometry(dom,z)[0], z)[0] for z in zscan])
        for level in (.45,.5,.55):
            count = np.count_nonzero(np.diff(vals > level))
            if count != 1:
                raise ArithmeticError(f'found {count} crossings at eta {eta}, level {level}')
        return (zscan[0], zscan[-1])
    etas = (0., .25,.5,.75,1.)
    bands = {eta:bracket(eta) for eta in etas}
    points = [sample(eta, .5, bands[eta]) for eta in etas]
    root, weights = leggauss(16)
    work = [.05 * weights @ [sample(eta, .5+.05*t, bands[eta])['mu_per_n'] for t in root]
            for eta in (.25,.75)]
    return dict(midpoint=points, leg_chi_interval=[.45,.55], leg_eta=[.25,.75],
                leg_work=work, signed_closed_work=work[0]-work[1])


def run(args):
    path = args.source/'field.npz'
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != json.loads(ARCHIVE.read_text())['output_field_sha256'][str(path)]:
        raise ValueError('carrier archive digest disagrees')
    data = json.loads((args.source/'stdout').read_text().splitlines()[0])
    dom = FlaredDomain(data['n'], data['profile'], data['roof'], data['feed_far'])
    with np.load(path) as saved:
        np.testing.assert_array_equal(dom.points, saved['points'])
        field = saved['psi']
    model = Order(dom, field, args.delta)
    emit(test='domain', n=data['n'], profile=data['profile'], roof=data['roof'],
         source=str(path), source_sha256=digest, delta=args.delta, sigma=SIGMA, nu=NU,
         mobility=NU/(6*SIGMA*args.delta), dofs=len(field), source_flux=data['bulk'])
    for seed in args.seeds:
        order, solver = model.minimize(seed)
        np.savez_compressed(args.outdir/f'order_seed{seed}.npz', chi=order)
        emit(test='stationary', seed=seed, **solver, **diagnostics(dom, field, order, args.delta))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--delta', type=float, required=True)
    parser.add_argument('--seeds', nargs='+', type=float, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    run(parser.parse_args())
