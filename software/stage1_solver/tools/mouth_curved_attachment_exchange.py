"""Dilute tracer occupancy at common bath activity above a clamped triangular host.

K is the integral of the Boltzmann factor over a finite height window and
one primitive cell, divided by cell area. The reservoir number density times
K gives surface density. All lengths and kBT in the integral are microscopic.
The polar order shift is taken from the archived weak-flow mouth solution.
"""
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import logsumexp

from mouth_curved_companion_order import energy_derivatives
from mouth_curved_site_front import sha
from mouth_surface_site_transfer import BASIS, Surface, sheet_mean
from neutral_patch_sites import BETWEEN, TEMPERATURE, emit, morse


RECORDS = Path('research/pde_ledger_v3/_measurements')
Z_EDGES = (0., .25, .54, .64, .69, .72, .75, .78, .81, .84, .87, .90,
           .94, 1., 1.15, 1.3, 1.5, 2.)
WINDOWS = ((.25, 1.5), (.25, 2.), (0., 2.))
LEVELS = ((10, 14, 14), (15, 20, 18))


def archived_rows(name, run):
    path = RECORDS / (name+'.json')
    r = json.loads(path.read_text())['guarded_runs'][run]
    if r['guard_exit'] or r['child_exit'] or not r['limits_verified']:
        raise ValueError(f'Unverified input {path}/{run}')
    emit(test='input', record=str(path), sha256=sha(path), run=run,
         source_stdout_sha256=r['stdout_sha256'])
    return list(map(json.loads, r['literal_stdout'].splitlines()))


def composite_gauss(edges, degree):
    node, weight = leggauss(degree)
    a, b = np.array(edges[:-1]), np.array(edges[1:])
    return (((a+b)[:, None]+(b-a)[:, None]*node)/2).ravel(), (
        (b-a)[:, None]*weight/2).ravel()


class Cell:
    def __init__(self, surface, levels):
        radial, angular, normal = levels
        self.surface = surface
        # The lattice permutes the twelve triangles within one primitive
        # area. In each, start at the hollow, continue toward the atom/edge.
        s, ds = composite_gauss((0., .002, .01, .04, .14, .4, 1.), radial)
        t, dt = composite_gauss((0., 1.), angular)
        hollow = surface.pitch*BASIS.sum(axis=0)/3
        edge = surface.pitch*BASIS[0]/2
        xy = hollow + s[:, None, None]*(t[None, :, None]*edge-hollow)
        self.xy = xy.reshape(-1, 2)
        self.weight = np.outer(2*s*ds, dt).ravel()
        np.testing.assert_allclose(self.weight.sum(), 1., atol=1e-14, rtol=0.)
        delta = self.xy[:, None, :]-surface.origins
        self.r2 = np.sum(delta**2, axis=-1)
        self.z, self.dz = composite_gauss(Z_EDGES, normal)
        self.coherent = np.array([self.at_height(z) for z in self.z])
        self.mean = sheet_mean(self.z, surface.pitch)[0][:, None]
        self.log_measure = np.log(self.weight[None, :]*self.dz[:, None])
        self.boundaries = {z:self.at_height(z) for z in {v for pair in WINDOWS for v in pair}}

    def at_height(self, z):
        # These are exactly the finite image list and pair energy of Surface.
        r = np.sqrt(self.r2+z*z)
        return morse(r, *BETWEEN)[0].sum(axis=1)

    def boundary(self, z, order, mixed):
        v = order*self.boundaries[z] + (1-order)*mixed*sheet_mean(z, self.surface.pitch)[0]
        return logsumexp(-v/TEMPERATURE, b=self.weight)

    def check(self):
        # Full periodic potential check, including the triangle decomposition.
        # Its lateral mean is the exact uniform sheet, checked at probe height.
        z = 1.1
        integrated = self.weight @ self.at_height(z)
        exact = sheet_mean(z, self.surface.pitch)[0]
        indices = [0, len(self.xy)//2, len(self.xy)-1]
        direct = [self.surface.evaluate(np.r_[self.xy[i], z], 1., False)[0] for i in indices]
        vector = self.at_height(z)[indices]
        emit(test='cell', nodes_lateral=len(self.xy), nodes_normal=len(self.z),
             area_weight=float(self.weight.sum()), sheet_mean=exact, lattice_mean=integrated,
             full_surface_difference=(np.array(direct)-vector).tolist())
        np.testing.assert_allclose(direct, vector, atol=1e-11, rtol=0.)
        np.testing.assert_allclose(integrated, exact, atol=1e-9, rtol=0.)


class Layer:
    def __init__(self, cell, lo, hi, mixed):
        ids = (cell.z>lo) & (cell.z<hi)
        self.intercept = mixed*cell.mean[ids]
        self.slope = cell.coherent[ids]-self.intercept
        self.log_measure = cell.log_measure[ids]

    def partition(self, order, moments=False):
        log_weight = self.log_measure-(self.intercept+order*self.slope)/TEMPERATURE
        log_K = logsumexp(log_weight)
        if not moments:
            return log_K
        p = np.exp(log_weight-log_weight.max())
        p /= p.sum()
        mean = np.sum(p*self.slope)
        variance = np.sum(p*(self.slope-mean)**2)
        return log_K, mean, variance

    def signed_shell(self, order, dipole, speed, base):
        x, weights = leggauss(10)
        # Independent reservoir exchange at each patch. Preserve the mean
        # population change before normalizing the observable 3 <cos(theta)>.
        n = np.exp([self.partition(order+speed*dipole*xx)-base for xx in x])
        total = weights @ n/2
        return 3*(weights @ (x*n))/(2*total), total


def measure(cell, binding, force, scale, fine, cold):
    loc, name = binding['location'], binding['law']
    mixed = name=='retain_mean_attraction'
    order, dipole = binding['order_average'], binding['order_dipole_overlap_per_U']
    frozen_h = force['held_wall_and_attached_source']
    isolated = binding['isolated_holding']
    all_results = []
    for lo, hi in WINDOWS:
        layer = Layer(cell, lo, hi, mixed)
        log_K, mean, variance = layer.partition(order, True)
        first = scale*mean
        second = -scale*variance/TEMPERATURE
        log_derivative = -mean/TEMPERATURE
        beta = log_derivative*dipole
        # The wall source is -dF/dchi per occupant. Its change includes
        # both the response per occupied site and its independent abundance.
        count_coefficient = first*isolated/3
        fixed = first*frozen_h + count_coefficient*(second/first)*dipole
        net = fixed + count_coefficient*beta
        lower = cell.boundary(lo, order, mixed)
        upper = cell.boundary(hi, order, mixed)
        row = dict(law=name, location=loc, z_lo=lo, z_hi=hi, chi=order,
                   order_dipole_per_U=dipole, log_K=log_K, K=np.exp(log_K),
                   dlogK_dchi=log_derivative, free_energy_slope_specific=first,
                   free_energy_curvature_specific=second,
                   log_boundary_lower=lower, log_boundary_upper=upper,
                   count_dipole_per_U=beta, fixed_count_order_force=fixed,
                   count_force_coefficient=count_coefficient,
                   equilibrium_order_force=net, zero_force_count_dipole=-fixed/count_coefficient,
                   cold_rigid_slope_specific=cold[0], cold_rigid_curvature_specific=cold[1],
                   cold_relaxed_slope_specific=binding['energy_slope'],
                   cold_relaxed_curvature_specific=binding['energy_curvature'])
        emit(test='partition', resolution='fine' if fine else 'coarse', **row)
        all_results.append(row)
        if fine and (lo,hi)==WINDOWS[0]:
            # The small same-law curvature subtracts almost equal mean
            # slopes. A fourth-order stencil permits a resolved step.
            h = .0001
            moments = np.array([layer.partition(order+j*h, True)[:2] for j in (-2,-1,1,2)])
            fd1, fd2 = np.array([1,-8,8,-1]) @ moments/(12*h)
            fd2 *= scale
            emit(test='order_difference', law=name, location=loc, stencil_step=h,
                 log_K_difference_derivative=fd1, moment_derivative=log_derivative,
                 slope_difference_curvature=fd2, variance_curvature=second)
            np.testing.assert_allclose([fd1,fd2], [log_derivative,second], rtol=3e-8, atol=1e-10)
            for speed in (.002, .001):
                plus, n1 = layer.signed_shell(order, dipole, speed, log_K)
                minus, n2 = layer.signed_shell(order, dipole, -speed, log_K)
                measured = (plus-minus)/(2*speed)
                emit(test='signed_exchange', law=name, location=loc, speed=speed,
                     positive_count_dipole=plus, negative_count_dipole=minus,
                     mean_density_over_undriven_plus=n1, mean_density_over_undriven_minus=n2,
                     measured_slope=measured, weak_drive_slope=beta)
                np.testing.assert_allclose(n1, n2, rtol=3e-13)
                np.testing.assert_allclose(measured, beta, rtol=2e-5)
    return all_results


def run():
    sample, = [r for r in archived_rows('mouth_bound_site_transfer', 'primary') if r['test']=='sample']
    mouth = archived_rows('mouth_curved_companion_order', 'mesh48')
    bindings = [r for r in mouth if r['test']=='binding']
    forces = {r['location']:r for r in mouth if r['test']=='force'}
    surface = Surface(sample['side'], sample['shell'], sample['pitch'], True)
    scale = sample['scale']
    if len(bindings)!=4 or len(forces)!=2:
        raise ValueError('Incomplete shell inputs')
    emit(test='geometry', pitch=surface.pitch, cell_area=abs(np.linalg.det(BASIS))*surface.pitch**2,
         kBT=TEMPERATURE, specific_scale=scale, specific_kBT=TEMPERATURE*scale,
         image_count=len(surface.origins), windows=WINDOWS, quadrature_levels=LEVELS)
    cold = [energy_derivatives(surface, scale, row['order_average'],
                               row['law']=='retain_mean_attraction') for row in bindings]
    results = []
    for i, level in enumerate(LEVELS):
        cell = Cell(surface, level)
        cell.check()
        results.append([r for b, c in zip(bindings, cold, strict=True) for r in measure(
            cell, b, forces[b['location']], scale, i==1, c)])
    for a,b in zip(*results, strict=True):
        keys = ('log_K', 'dlogK_dchi', 'free_energy_curvature_specific',
                'log_boundary_lower', 'log_boundary_upper')
        delta = {k:a[k]-b[k] for k in keys}
        emit(test='refinement', law=b['law'], location=b['location'], z_lo=b['z_lo'],
             z_hi=b['z_hi'], coarse_minus_fine=delta)
        if (abs(delta['log_K'])>1e-5 or abs(delta['dlogK_dchi'])>3e-4
                or abs(delta['free_energy_curvature_specific'])>1e-5
                or max(abs(delta[k]) for k in keys[-2:])>.02):
            raise ArithmeticError('Partition or boundary quadrature unresolved')


if __name__=='__main__':
    run()
