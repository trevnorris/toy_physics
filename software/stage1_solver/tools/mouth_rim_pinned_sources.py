"""Physical constraint data for a local, axial carrier mesh displacement at R4.

A duct-interior extension may move the mesh around fixed material sites, but
cannot move the bonded wall through those sites. Use a zero-trace difference
of two extensions for a fixed physical solid. On the other bank use the
prescribed receiving-order field evaluated at the moved physical points.
"""
import numpy as np

from mouth_last_ring_solid import interpolation
from mouth_rim_shared_crystal import TerminalQuadrature, assemble, deposit


def amplitude(bank, edge, knots, x, z):
    """Product of the axial tent and the radial taper at original chart points."""
    lo, mid, hi = knots
    hat = np.maximum(0., np.minimum((z - lo) / (mid - lo), (hi - z) / (hi - mid)))

    def taper(bound):
        if bank == 1:
            return np.clip(1. - x / bound, 0., 1.)
        return np.clip((x - bound) / (1. - bound), 0., 1.)

    if np.ndim(edge) == 0:
        return taper(edge) * hat
    return (taper(edge[0]) - taper(edge[1])) * hat


class Pinned:
    def __init__(self, model, state_path, wall_path, reference_order, baseline_load):
        self.model = model
        self.state_path = state_path
        self.wall_path = wall_path
        self.order = np.asarray(reference_order)
        self.load = np.asarray(baseline_load)
        self.sites = None

    def ref_change(self, bank, edge, knots, step):
        m = self.model
        x, z = np.meshgrid(m.rr[bank], m.z[1], indexing='ij')
        b = m.shape.radius(z)[0]
        r = b + x * (28. - b) if bank else x * b
        new_z = z + step * amplitude(bank, edge, knots, x, z)
        if bank:
            profile = m.terminal_profile.outside
            delta = profile(r, new_z, False) - profile(r, z, False)
        else:
            new_x = r / m.shape.radius(new_z)[0]
            if np.min(new_x) < -1e-13 or np.max(new_x) > 1. + 1e-13:
                raise ValueError('Physical duct pin would need material outside its domain')
            nt, nx = self.order.shape
            axes = (np.linspace(0., 1., nt), np.linspace(0., 1., nx), self.order)
            ref = interpolation(*axes, z / m.shape.length, x)[0]
            new = interpolation(*axes, new_z / m.shape.length, np.clip(new_x, 0., 1.))[0]
            delta = new - ref
        delta[new_z == z] = 0.
        return delta

    def material_sites(self):
        if self.sites is None:
            m = self.model
            quadrature = TerminalQuadrature(self.state_path, self.wall_path, m.rr[0], m.z[1], 64)
            records = []
            for block in quadrature.blocks():
                records.append(np.stack((block['r'], block['z'], block['measure'] * block['W'])))
            self.sites = np.concatenate(records, axis=1)
        return self.sites

    def crystal_change(self, edge, knots, step):
        """Integrate at the same physical sites and strain; invert the moved FE map."""
        if np.ndim(edge) == 0:
            raise ValueError('Axial wall displacement would move the bond relative to the fixed crystal')
        m = self.model
        r, z, weight = self.material_sites()
        x = r / m.shape.radius(z)[0]
        active = amplitude(0, edge, knots, x, z) != 0.
        r, z, weight, x = (v[active] for v in (r, z, weight, x))
        z0 = z.copy()
        for _ in range(10):
            x0 = r / m.shape.radius(z0)[0]
            z0 = z - step * amplitude(0, edge, knots, x0, z0)
        x0 = r / m.shape.radius(z0)[0]
        inverse_error = np.max(abs(z0 + step * amplitude(0, edge, knots, x0, z0) - z))
        if inverse_error > 2e-13 or x0.max() > 1. or x0.min() < 0.:
            raise ArithmeticError('Fixed physical crystal sites no longer covered by the carrier map')
        plane = np.zeros((len(m.rr[0]), len(m.z[1])))
        deposit(plane, m.rr[0], m.z[1], x0, z0, weight)
        deposit(plane, m.rr[0], m.z[1], x, z, -weight)
        if abs(plane.sum()) > 2e-11:
            raise ArithmeticError('Lost partition of unity while transferring fixed material')
        return assemble(m, plane), dict(site_count=int(active.sum()), mass_moment_error=float(plane.sum()),
                                        inverse_map_error=float(inverse_error))

    def corrections(self, bank, edge, knots, step):
        if bank == 0 and np.ndim(edge) == 0 and step != 0.:
            raise ValueError('A physically fixed, bonded duct face cannot follow the imposed wall shift')
        if bank not in (0, 1):
            raise ValueError('Expected either duct or receiving bank')
        delta = self.ref_change(bank, edge, knots, step)
        load = self.load
        info = dict(max_reference_order_change=float(np.max(abs(delta))))
        if bank == 0 and step != 0.:
            change, metrics = self.crystal_change(edge, knots, step)
            load = load + change
            info.update(metrics)
        return {(1, bank): delta}, load, info
