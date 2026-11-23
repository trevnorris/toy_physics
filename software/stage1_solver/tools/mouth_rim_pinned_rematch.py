"""Stationary signed wall and pure-mesh probes of the held one-quantum R4 rim.

The energy is the original two-point carrier FE functional, including its
oriented vortex surface term and registered material load. All three fields
are varied wherever the original R4 constitutive law makes them mobile.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from mouth_exchange_carrier import BIAS
from mouth_joined_front_stationarity import solve as stationary
from mouth_rim_consistent_fields import from_archive, material_load
from mouth_rim_pinned_sources import Pinned
from mouth_supercore_joined_front import file_sha


def radial_weights(x, side, edge):
    if side:
        value = np.clip(1 - x / edge, 0, 1)
        slope = np.where(x < edge, -1 / edge, 0)
    else:
        value = np.clip((x - edge) / (1 - edge), 0, 1)
        slope = np.where(x > edge, 1 / (1 - edge), 0)
    return value, slope


def depth_weights(z, knots):
    lo, mid, hi = knots
    value = np.maximum(0, np.minimum((z - lo) / (mid - lo), (hi - z) / (hi - mid)))
    slope = np.where((z > lo) & (z < mid), 1 / (mid-lo),
                     np.where((z > mid) & (z < hi), -1 / (hi-mid), 0))
    return value, slope


def slit_matrix_delta(model, side, knots, relative):
    # On the side itself the axial map is z' = z + step*tent(z). Its
    # contribution to the radial oriented area is exactly 1 + step*tent'.
    matrix = sparse.csr_matrix(model.B.shape)
    if relative:
        return matrix
    z_nodes = model.z[1]
    dz = np.diff(z_nodes)
    cells = np.flatnonzero((z_nodes[:-1] >= knots[0]-1e-13) &
                           (z_nodes[1:] <= knots[2]+1e-13))
    da = np.diff(model.a)[:, None]
    face = model.ids[1, side][0 if side else -1]
    ids = np.stack([face[i:i+len(da), cells+j] for i in (0, 1)
                    for j in (0, 1)], axis=-1)
    quadrature = (1 + np.array([-1., 1.])/np.sqrt(3.))/2
    for angular in quadrature:
        sn = np.sin(model.a[:-1, None]+angular*da)
        for depth in quadrature:
            z = z_nodes[cells] + depth*dz[cells]
            r = model.shape.radius(z)[0]
            F, Fprime = model._vortex(r*sn)
            weights = ((1-2*side)*np.pi*r*r*sn**2*da*dz[cells] *
                       F*Fprime*depth_weights(z, knots)[1])
            basis = np.kron([1-angular, angular], [1-depth, depth])
            local = weights[..., None, None]*np.outer(basis, basis)
            rows = np.broadcast_to(ids[..., :, None], local.shape)
            cols = np.broadcast_to(ids[..., None, :], local.shape)
            matrix += sparse.coo_matrix((local.ravel(), (rows.ravel(), cols.ravel())),
                                        shape=matrix.shape).tocsr()
    return matrix


class MovedMesh:
    def __init__(self, model, pinned, side, edge, knots):
        self.model, self.pinned = model, pinned
        self.side, self.edge, self.knots = side, edge, knots
        self.baseline = model.D, model.G, model.ref, model.B, model.material_load
        D = model.D
        r, z = D[:, 1], D[:, 3]
        radius, slope = model.shape.radius(z)
        width = 28-radius if side else radius
        x = (r-radius)/width if side else r/radius
        self.relative = isinstance(edge, tuple)
        first = edge[0] if self.relative else edge
        v, dv = radial_weights(x, side, first)
        if self.relative:
            other, dother = radial_weights(x, side, edge[1])
            v -= other
            dv -= dother
        tent, dt = depth_weights(z, knots)
        inside = ((r > radius) if side else (r < radius)) & (tent > 0)
        # Force bank selection geometrically only in the terminal band, where
        # the independently measured polar deformation has exactly zero support.
        self.eta = np.where(inside, tent*v, 0)
        self.er = np.where(inside, tent*dv/width, 0)
        self.ez = np.where(inside, dt*v-slope*(1-x if side else x)*self.er, 0)
        self.wall = slit_matrix_delta(model, side, knots, self.relative)

    def baseline_model(self):
        m = self.model
        m.D, m.G, m.ref, m.B, m.material_load = self.baseline
        m.hq = m.D[:, 0]

    def apply(self, step):
        m = self.model
        self.baseline_model()
        D, G, reference, B, load = self.baseline
        if step == 0:
            self.unpinned_reference = reference
            return {'minimum_axial_jacobian': 1.}
        delta_ref, physical_load, info = self.pinned.corrections(
            self.side, self.edge, self.knots, step)
        jacobian = 1+step*self.ez
        if jacobian.min() < .95:
            raise ArithmeticError('Excessive axial mesh strain')
        scale = 1/jacobian
        shear = -step*self.er*scale
        m.G = [(G[0]+G[2].multiply(shear[:, None])).tocsr(), G[1],
               G[2].multiply(scale[:, None]).tocsr()]
        m.D = D.copy()
        m.D[:, 0] *= jacobian
        m.D[:, 3] += step*self.eta
        m.D[:, 9] *= jacobian
        m.hq = m.D[:, 0]
        m.B = B + step*self.wall
        self.unpinned_reference = reference.copy()
        self.unpinned_reference[:, 1] += shear*reference[:, 3]
        self.unpinned_reference[:, 3] *= scale
        m.ref = self.unpinned_reference.copy()
        for key, increments in delta_ref.items():
            # An independently prescribed order may have different limits on
            # opposite banks of the slit, including at shared field vertices.
            # Keep its contributions restricted to its own element block.
            ids, start = m.blocks[key]
            end = start+8
            offset = sum(8*len(n) for (n,k) in m.blocks.values() if k < start)
            rows = slice(offset, offset+len(ids)*(end-start))
            values = np.zeros(m.size)
            values[m.ids[key]] = increments[:, None, :]
            m.ref[rows, 0] += m.V[rows]@values
            for i, derivative in enumerate(m.G):
                m.ref[rows, 1+i] += derivative[rows]@values
        m.material_load = physical_load
        return {**info, 'minimum_axial_jacobian': float(jacobian.min())}

    def energies(self, u):
        m = self.model
        energy = m.evaluate(u)[0]
        # Isolate the physical order pin at the *same re-equilibrated field*
        # before attributing any remaining energy difference to the wall.
        pinned_reference = m.ref
        try:
            m.ref = self.unpinned_reference
            carried_order = m.evaluate(u)[0]
        finally:
            m.ref = pinned_reference
        crystal_correction = (m.material_load-self.baseline[4])@u[2]
        return {'energy': float(energy), 'order_pin_relative_to_carried': float(energy-carried_order),
                'crystal_pin_relative_to_fixed_nodal_load': float(crystal_correction)}


def record(folder, name, data):
    path = folder/(name+'.json')
    tmp = path.with_suffix('.pending')
    tmp.write_text(json.dumps(data, allow_nan=False, indent=2)+'\n')
    tmp.replace(path)


def rematch(m, initial, rest, tolerance):
    m.seed = initial.copy()
    m.seed[0, ~m.masks[0]] = 1.
    m.seed[1, m.feed] = 0.
    m.seed[1, m.bulk] = 0. if rest else BIAS
    u, info = stationary(m, rest, tolerance=tolerance)
    energy, gradient = m.evaluate(u)
    residual = m.vars(gradient)
    inverse = m.norm(residual)
    stationarity = float(np.max(np.abs(inverse)))
    incoming = -float(gradient[1, m.feed].sum())
    outgoing = float(gradient[1, m.bulk].sum())
    metric_energy = float(.5*residual@inverse)
    if (not np.isfinite(u).all() or stationarity >= tolerance or
            (not rest and min(incoming, outgoing) <= 0) or
            metric_energy < -1e-16 or abs(incoming-outgoing) >= 3e-5):
        raise ArithmeticError(f'Failed stationary field: {stationarity}, currents {incoming}, {outgoing}')
    return u, {**info, 'rechecked_energy': float(energy), 'stationarity': stationarity,
               'residual_metric_energy': metric_energy, 'feed': incoming, 'bulk': outgoing,
               'port_error': abs(incoming-outgoing)}


def run(a):
    base, tr, m, archived, level, digest, _, _ = from_archive(
        a.state, a.wall, a.fields, True)
    m.material_load = material_load(m, a.load, a.wall)
    with np.load(a.fields) as saved:
        if str(saved['load_sha']) != file_sha(a.load):
            raise ValueError('Field states and registered load do not match')
    with np.load(a.axes) as saved:
        if saved['level'] != 'double':
            raise ValueError('Expected original doubled axes for identical physical motion')
        knots = saved['duct_z'][[-9, -5, -1]]
        receiver = tuple(saved['outside_x'][[4, 8]])
        duct = tuple(saved['inside_x'][[-5, -9]])
    pinned = Pinned(m, a.state, a.wall, base.c, m.material_load)
    cases = [('four', 1, receiver[0]), ('eight', 1, receiver[1]), ('four_minus_eight', 0, duct)]
    meta = dict(field_sha256=digest, wall_sha256=file_sha(a.wall), load_sha256=file_sha(a.load),
                state_sha256=file_sha(a.state), reference_axes_sha256=file_sha(a.axes),
                grid=level, terminal_knots=knots.tolist(), receiver_bounds=receiver, duct_bounds=duct,
                full_step=a.step, stationarity_tolerance=a.tolerance)
    record(a.folder, 'inputs', meta)
    zero, initial_rows = [], []
    for rest, u in enumerate(archived):
        settled, info = rematch(m, u, bool(rest), a.tolerance)
        zero.append(settled)
        initial_rows.append(info)
    record(a.folder, 'zero', initial_rows)
    np.savez_compressed(a.folder/'zero.npz', wound=zero[0], rest=zero[1])
    rows = []
    for name, side, edge in cases:
        chart = MovedMesh(m, pinned, side, edge, knots)
        for mult, label in [(1,'plus'),(-1,'minus'),(.5,'plus_half'),(-.5,'minus_half')]:
            step = mult*a.step
            pin_info = chart.apply(step)
            fields, pair = [], []
            for rest, u0 in enumerate(zero):
                u, info = rematch(m, u0, bool(rest), a.tolerance)
                fields.append(u)
                pair.append({**info, **chart.energies(u), 'frozen_energy': float(m.evaluate(u0)[0]),
                             'response_max': np.max(abs(u-u0), axis=1).tolist()})
            file = a.folder/(name+'_'+label+'.npz')
            np.savez_compressed(file, wound=fields[0], rest=fields[1])
            value = dict(name=name, side=side, step=step, pin_info=pin_info, pair=pair, file=file.name,
                         sha256=file_sha(file))
            record(a.folder, name+'_'+label, value)
            rows.append(value)
            print(f'{level} {name} {label}: stationary wound/rest', flush=True)
        chart.baseline_model()
    record(a.folder, 'result', dict(inputs=meta, zero=initial_rows, signed=rows))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('state', 'wall', 'fields', 'load', 'axes', 'folder'):
        p.add_argument(key, type=Path)
    p.add_argument('--step', type=float, default=2e-5)
    p.add_argument('--tolerance', type=float, default=5e-9)
    run(p.parse_args())
