"""Evaluate the original free-order FE energy with differences on constant plateaus.

Near chi=1, tiny roundoff in K@ones and 1-(N@chi) can dominate an absolute
weak residual. Evaluate the same gradient on pairs of nodal differences and
interpolate ordered and disordered fractions independently. Then check the
solution against the unmodified assembly, including its full Newton correction.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from mouth_flared_receiver import FlaredDomain
from mouth_free_order_diffuse import ARCHIVE, NU, SIGMA, Order, diagnostics, emit, positive_factor


class DifferenceOrder(Order):
    def __init__(self, dom, carrier, delta):
        super().__init__(dom, carrier, delta)
        self.row = np.repeat(np.arange(len(carrier)), np.diff(self.stiff.indptr))
        self.col = self.stiff.indices
        self.quad_sum = dom.N @ np.ones(len(carrier))
        self.Nf = sparse.diags(1/self.quad_sum) @ self.Nf
        corrected = self.stiff - sparse.diags(self.stiff @ np.ones(len(carrier)))
        self.K = NU*corrected[dom.free][:,dom.free]
        self.mass = self.mass_weight(self.wt)

    def energy(self, x):
        diff = x[self.col]-x[self.row]
        c = (self.dom.N @ x)/self.quad_sum
        d = (self.dom.N @ (1-x))/self.quad_sum
        return -.25*NU*(self.stiff.data @ diff**2) + .5*self.scale*(self.wt @ (c*d)**2)

    def measure(self, x):
        c = (self.dom.N @ x)/self.quad_sum
        d = (self.dom.N @ (1-x))/self.quad_sum
        diff = x[self.col]-x[self.row]
        grad = np.add.reduceat(self.stiff.data*diff, self.stiff.indptr[:-1])
        residual = NU*grad[self.dom.free] + self.Nf.T @ (self.scale*self.wt*c*d*(d-c))
        jac = self.K + self.mass_weight(self.scale*self.wt*(1-6*c*d))
        return residual, jac


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--seed', required=True, type=float)
    parser.add_argument('--delta', required=True, type=float)
    parser.add_argument('--outdir', required=True, type=Path)
    args = parser.parse_args()
    path = args.source/'field.npz'
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != json.loads(ARCHIVE.read_text())['output_field_sha256'][str(path)]:
        raise ValueError('carrier archive digest disagrees')
    state = json.loads((args.source/'stdout').read_text().splitlines()[0])
    dom = FlaredDomain(state['n'], state['profile'], state['roof'], state['feed_far'])
    with np.load(path) as saved:
        np.testing.assert_array_equal(dom.points,saved['points'])
        carrier = saved['psi']
    model = DifferenceOrder(dom,carrier,args.delta)
    x, solver = model.minimize(args.seed)
    # Evaluate the original code and its Newton correction on the recovered field.
    original = Order(dom,carrier,args.delta)
    residual, jac = original.measure(x)
    fac = positive_factor(jac)
    correction = np.max(abs(fac.solve(residual))) if fac is not None else None
    report = dict(original_max_weak_residual=np.max(abs(residual)), original_max_Newton_step=correction,
                  max_basis_sum_error=np.max(abs(model.quad_sum-1)),
                  original_minus_difference_energy=original.energy(x)-model.energy(x))
    emit(test='recovery_validation', source=str(path), source_sha256=digest, seed=args.seed,
         n=state['n'], profile=state['profile'], delta=args.delta, nu=NU, sigma=SIGMA, **report)
    if fac is None or correction > 1e-8 or report['original_max_weak_residual'] > 1e-9:
        raise ArithmeticError('original assembly failed recovery cross-check')
    np.savez_compressed(args.outdir/f'order_seed{args.seed}.npz', chi=x)
    emit(test='stationary', seed=args.seed, **solver, **diagnostics(dom,carrier,x,args.delta))
