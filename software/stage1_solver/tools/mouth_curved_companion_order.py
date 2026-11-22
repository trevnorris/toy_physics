"""Companion-induced holding load of a dilute spherical population on the curved wall.

U is the imposed feeding-bath phase gradient in +x. The centers stay attached
at fixed depth. h is their total coefficient in -h times the shell-average
order. The existing carrier has no reciprocal response to this trial site law.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_carrier_lateral_response import Dipole
from mouth_curved_site_front import WIDTH, prepare, sha, source
from mouth_curved_site_pair import BASIS, material_calibration
from mouth_curved_site_stress import wall_integral
from mouth_empty_bend_stiffness import quadrature
from mouth_flared_receiver import EXIT, FLARE_START, ROUNDING
from mouth_free_order_diffuse import NU, SIGMA, Order, emit, positive_factor
from mouth_curved_companion_stress import anchor_gradient, budget


PRIOR = Path('research/pde_ledger_v3/_measurements/mouth_curved_site_front.json')
SAVED = Path('_scratch/s11c/curved_site01')
KAPPA = 6 * SIGMA * WIDTH


class Response:
    def __init__(self, dom, psi, chi):
        self.dom = dom
        self.order = Order(dom, psi, WIDTH)
        self.dpsi, residual = Dipole(dom, psi).steady(1., 0.)
        phase = dom.N @ psi
        numerator = dom.N @ self.dpsi
        radius, _, self.gr, self.gw = quadrature(dom)
        self.drive = [((op @ self.dpsi) / phase
                       - (numerator / phase) * (op @ (psi-1)) / phase).imag
                      for op in (self.gr, self.gw)]
        self.angular = self.order.mass_weight(2*NU*self.order.wt/radius**2)
        self.regular = dom.points[dom.free, 0] > 0
        self.active = np.flatnonzero(dom.free)[self.regular]
        empty_residual, jac = self.order.measure(chi)
        if np.max(abs(empty_residual)) > 5e-9:
            raise ArithmeticError('Empty order residual')
        self.factor = self.dipole_factor(jac)
        self.b = self.solve(self.factor, -self.advection(chi))
        np.testing.assert_allclose(self.dpsi[dom.feed], 1j*psi[dom.feed]*dom.points[dom.feed,0],
                                   atol=1e-12, rtol=0)
        np.testing.assert_allclose(self.dpsi[dom.bulk], 0., atol=1e-12, rtol=0)
        emit(test='incident', carrier_max_free_residual=residual,
             empty_order_max_free_residual=float(max(abs(empty_residual))),
             feed_phase_gradient=1., receiving_phase_gradient=0.,
             angular_basis='cosine toward the feed velocity; zero first-order monopole')

    def dipole_factor(self, jac):
        mat = (jac+self.angular)[self.regular][:,self.regular]
        factor = positive_factor(mat)
        if factor is None:
            raise ArithmeticError('Unstable dipole')
        return mat, factor

    def advection(self, field):
        return self.drive[0]*(self.gr @ field) + self.drive[1]*(self.gw @ field)

    def solve(self, system, rhs_q):
        mat, factor = system
        rhs = (self.order.Nf.T @ (self.order.wt*rhs_q))[self.regular]
        out = np.zeros(len(self.dom.points))
        out[self.active] = factor.solve(rhs)
        residual = np.linalg.norm(mat@out[self.active]-rhs)/np.linalg.norm(rhs)
        if not np.isfinite(residual) or residual > 2e-8:
            raise ArithmeticError(f'Dipole relative residual {residual}')
        return out

    def mixed(self, chi, a):
        N = self.dom.N
        quadratic = self.order.scale*(12*(N@chi)-6)*(N@a)*(N@self.b)
        return self.solve(self.factor, -self.advection(a)-quadratic)

    def finite_shell(self, chi, weighted_source, strength):
        fields = chi.copy()
        rhs = strength * (NU/KAPPA) * weighted_source[self.dom.free]
        points = self.dom.points[self.dom.free]
        near = (points[:,0]<2.) & (points[:,1]>0.) & (points[:,1]<3.)
        for _ in range(15):
            residual, jac = self.order.measure(fields)
            error = residual-rhs
            if max(abs(error)) < 3e-10 and max(abs(error[near])) < 3e-12:
                break
            factor = positive_factor(jac)
            if factor is None:
                raise ArithmeticError('Finite shell left stable deep branch')
            fields[self.dom.free] -= factor.solve(error)
        else:
            raise ArithmeticError(f'Finite shell residual {max(abs(error))}, near {max(abs(error[near]))}')
        emit(test='finite_residual', h=strength, max_absolute=float(max(abs(error))),
             max_near_conversion=float(max(abs(error[near]))))
        dipole = self.solve(self.dipole_factor(jac), -self.advection(fields))
        return fields, dipole


def finite_check(response, chi, a, c, depth, width, psi, target):
    dom = response.dom
    load = source(dom, psi, chi, depth, width, 12)[0][0]
    intervals = [(0, (dom.blocks[0][0][0], 0.)), (2, (ROUNDING, FLARE_START)),
                 (2, (FLARE_START, EXIT))]
    def holding(z, b, strength):
        fields = np.column_stack((z, b))
        wall = sum(wall_integral(dom, fields, block, interval, KAPPA, WIDTH)['holding']
                   for block, interval in intervals)
        return wall-strength*anchor_gradient(dom,b,depth,width)
    empty = holding(chi,response.b,0.)
    for h in (.0002, .0001):
        states = [response.finite_shell(chi, load, s) for s in (h, -h)]
        positive, negative = [holding(z,b,s) for (z,b),s in zip(states,(h,-h))]
        da, dc = [(states[0][i]-states[1][i])/(2*h) for i in (0,1)]
        errors = [np.linalg.norm(measured-reference)/np.linalg.norm(reference)
                  for measured,reference in [(da,a),(dc,c)]]
        central = (positive-negative)/(2*h)
        emit(test='finite_shell', depth=depth, h=h, empty_holding_per_U=empty,
             occupied_holding_per_U=positive, reversed_coupling_holding_per_U=negative,
             positive_subtraction_per_h_U=(positive-empty)/h, centered_per_h_U=central,
             analytic_per_h_U=target, relative_l2_mean_and_mixed=errors)
        if max(errors)>0.008 or abs(central-target)>0.008*abs(target):
            raise ArithmeticError('Finite coupling check disagrees with mixed derivative')


def energy_derivatives(surface, scale, order, law):
    start = np.zeros(surface.dim)
    start[:2] = surface.pitch*BASIS.sum(axis=0)/3
    start[2] = .75
    x, _, slope, _ = surface.equilibrate(order, law, start)
    # The unrelaxed potential is affine in order. Internal displacements
    # relax, giving the envelope curvature -g_chi^T Hessian^{-1} g_chi.
    steps = 2e-5*np.eye(surface.dim)
    grad = lambda y,z: surface.evaluate(y,z,law)[1]
    H = np.column_stack([(grad(x+d,order)-grad(x-d,order))/(4e-5) for d in steps])
    H = (H+H.T)/2
    if np.linalg.eigvalsh(H)[0] <= 0 or slope>=0:
        raise ArithmeticError('Unstable adsorption branch or sign')
    g = grad(x,1.)-grad(x,0.)
    curvature = -g@np.linalg.solve(H,g)
    return scale*slope, scale*curvature


def run(args):
    dom, psi, chi = prepare(args.mesh)
    label = str(args.mesh) if args.width == .01 else 'width34'
    if label == 'width34' and args.mesh!=34:
        raise ValueError('Doubled footprint archived on mesh 34 only')
    fields_path = SAVED/('fields_wide.npz' if label=='width34' else f'fields{label}.npz')
    prior = json.loads(PRIOR.read_text())
    if sha(fields_path)!=prior['result_and_input_sha256'][str(fields_path)]:
        raise ValueError('Archived order response changed')
    stress = prior['guarded_runs']['stress'+label]
    assert stress['child_exit']==stress['guard_exit']==0 and stress['limits_verified']
    rows = [json.loads(line) for line in stress['literal_stdout'].splitlines()]
    isolated = {r['location']:r['held_wall_and_attached_source'] for r in rows
                if r['test']=='linear_bare_order_wall_load' and r['fields_sha256']==sha(fields_path)}
    response = Response(dom,psi,chi)
    surface,scale = material_calibration(5)
    with np.load(fields_path) as fields:
        np.testing.assert_array_equal(dom.points,fields['points'])
        np.testing.assert_array_equal(chi,fields['chi'])
        for loc in ('entry','half'):
            depth = float(fields['depth_'+loc])
            a = fields['eta_0_'+loc]
            c = response.mixed(chi,a)
            terms = budget(dom,psi,chi,a,response.b,c,response.dpsi,depth,args.width)
            fixed = terms['held_wall_and_attached_source']
            emit(test='force', mesh=args.mesh, location=loc, depth=depth, width=args.width,
                 fields_sha256=sha(fields_path), **terms)
            if args.finite_check:
                finite_check(response,chi,a,c,depth,args.width,psi,fixed)
            loads = source(dom,np.ones_like(psi),chi,depth,args.width,12)[0]
            np.testing.assert_allclose(loads[0].sum(),1.,atol=2e-12)
            mean = loads[0]@chi
            polar = (loads[1]@response.b)/3
            for law,name in [(False,'weaken_all_contacts'),(True,'retain_mean_attraction')]:
                mu1,mu2 = energy_derivatives(surface,scale,mean,law)
                strength_change = (mu2/mu1)*polar
                # A dipole in site strength or count has mean cos² = 1/3.
                count_coefficient = -abs(mu1)*isolated[loc]/3
                fixed_force = -abs(mu1)*fixed
                reweighted_force = count_coefficient*strength_change
                net = fixed_force+reweighted_force
                emit(test='binding', mesh=args.mesh, location=loc, law=name, width=args.width,
                     order_average=mean, order_dipole_overlap_per_U=polar,
                     energy_slope=mu1, energy_curvature=mu2, isolated_holding=isolated[loc],
                     fixed_strength_medium_force_per_mass_U=fixed_force,
                     changing_strength_medium_force_per_mass_U=reweighted_force,
                     fixed_count_medium_force_per_mass_U=net,
                     count_dipole_force_coefficient=count_coefficient,
                     fractional_count_dipole_per_U_for_zero_load=-net/count_coefficient)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh', type=int, choices=(26,34,48), required=True)
    p.add_argument('--width', type=float, choices=(.01,.02), default=.01)
    p.add_argument('--finite-check', action='store_true')
    run(p.parse_args())
