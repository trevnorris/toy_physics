"""Common bonded-rim displacement and two distinct carrier-mesh extensions.

The physical lift is the original near enclosure (axial band 2, radial band
4). All four blocks agree on the disk, all slit traces and both order pins.
The carrier coordinate choices differ within the terminal elements.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from mouth_last_ring_solid import interpolation
from mouth_rim_consistent_fields import from_archive, material_load
from mouth_rim_pinned_rematch import radial_weights, record, rematch
from mouth_rim_shared_crystal import TerminalQuadrature, assemble, deposit, moving_energy
from mouth_supercore_joined_front import file_sha


def linear_profile(z, knots, values):
    """Endpoint-clamped interpolation and its derivative away from vertices."""
    z = np.asarray(z)
    index = np.clip(np.searchsorted(knots, z, side="right") - 1, 0, len(knots) - 2)
    slopes = np.diff(values) / np.diff(knots)
    derivative = np.where((z > knots[0]) & (z < knots[-1]), slopes[index], 0.)
    return np.interp(z, knots, values), derivative


class SharedRimChart:
    def __init__(self, model, reference_native, axes_archive):
        self.shape = model.shape
        length = self.shape.length
        below, above = reference_native.z[1:]
        lower = np.flatnonzero(np.isclose(below, length-model.h, atol=1e-12, rtol=0)).item()
        upper = np.flatnonzero(np.isclose(above, length+model.h, atol=1e-12, rtol=0)).item()
        self.window = np.r_[below[lower-2:lower], above[upper+1:upper+3]]
        self.outer_taper = np.array(reference_native.rr[1][4:6])
        with np.load(axes_archive, allow_pickle=False) as grid:
            if grid["level"].item() != "double":
                raise ValueError("The extension must use the same doubled reference axes")
            self.tents = np.r_[grid["duct_z"][[-9, -5, -1]], grid["receiver_z"][[4, 8]]]
            self.edges = {n: (grid["inside_x"][-1-n], grid["outside_x"][n]) for n in (4, 8)}
        if abs(self.tents[2]-length) > 1e-13 or not np.all(np.diff(self.tents) > 0):
            raise ValueError("Malformed exit disk supports")
        if self.window[1] > self.tents[0] or self.window[2] < self.tents[-1]:
            raise ValueError("The carrier tent overlaps the physical taper")

    def weights(self, bank, r, z):
        radius, slope = self.shape.radius(z)
        width = 28.-radius if bank else radius
        x = (r-radius)/width if bank else r/radius
        xz = -(1-x if bank else x)*slope/width
        if bank:
            a, b = self.outer_taper
            p = np.clip((b-x)/(b-a), 0., 1.)
            px = -((x > a) & (x < b)).astype(float)/(b-a)
        else:
            p, px = np.ones_like(x), np.zeros_like(x)
        return x, width, xz, p, px

    def physical(self, bank, r, z):
        a, az = linear_profile(z, self.window, [0., 1., 1., 0.])
        _, width, xz, p, px = self.weights(bank, r, z)
        return a*p, a*px/width, az*p+a*px*xz

    def carrier(self, bank, r, z, n):
        a, az = linear_profile(z, self.window, [0., 1., 1., 0.])
        b, bz = linear_profile(z, self.tents, [0., 1., 0., 1., 0.])
        x, width, xz, p, px = self.weights(bank, r, z)
        g, gx = radial_weights(x, bank, self.edges[n][bank])
        derivative_x = (a-b)*px+b*gx
        return (a-b)*p+b*g, derivative_x/width, (az-bz)*p+bz*g+derivative_x*xz

    def stationary(self, bank, r, z):
        x, width, xz, _, _ = self.weights(bank, r, z)
        b, bz = linear_profile(z, self.tents[:3], [0., 1., 0.])
        g4, d4 = radial_weights(x, bank, self.edges[4][bank])
        g8, d8 = radial_weights(x, bank, self.edges[8][bank])
        if bank:
            return np.zeros_like(x), np.zeros_like(x), np.zeros_like(x)
        return b*(g4-g8), b*(d4-d8)/width, bz*(g4-g8)+b*(d4-d8)*xz

    @staticmethod
    def still(bank, r, z):
        return np.zeros_like(r), np.zeros_like(r), np.zeros_like(r)


def slit_derivative(model, physical):
    """Factorized vortex area contribution, for both independently held banks."""
    result = sparse.csr_matrix(model.B.shape)
    z = model.z[1]
    a = model.a
    for bank in (0, 1):
        face = model.ids[1, bank][0 if bank else -1]
        ids = np.stack((face[:-1, :-1], face[:-1, 1:], face[1:, :-1], face[1:, 1:]), axis=-1)
        rr = np.broadcast_to(ids[..., :, None], ids.shape + (4,)).ravel()
        cc = np.broadcast_to(ids[..., None, :], ids.shape + (4,)).ravel()
        gauss = (1+np.array([-1., 1.])/np.sqrt(3))/2
        for u in gauss:
            angle = a[:-1, None]+u*np.diff(a)[:, None]
            sine = np.sin(angle)
            for v in gauss:
                w = z[:-1]+v*np.diff(z)
                r, slope = model.shape.radius(w)
                _, dr, dw = physical(bank, r, w)
                derivative = dw+dr*slope
                factor, prime = model._vortex(r*sine)
                weight = (1-2*bank)*np.pi*r*r*sine*sine*factor*prime*derivative*np.diff(z)*np.diff(a)[:, None]
                basis = np.kron([1-u, u], [1-v, v])
                local = weight[..., None, None]*np.outer(basis, basis)
                result += sparse.coo_matrix((local.ravel(), (rr, cc)), shape=result.shape).tocsr()
    return result


class CrystalEnergy:
    """Material point motion is prescribed in space, independently of carrier width."""
    def __init__(self, model, chart, state, wall, reference_order):
        self.model = model
        self.material = TerminalQuadrature(state, wall, model.rr[0], model.z[1], 64,
                                           halfwidth=model.shape.length-chart.window[0])
        self.chart = chart
        self.order = reference_order
        self.cache = {}

    def increment(self, step):
        if step in self.cache:
            return self.cache[step]
        q = self.material
        plane = np.zeros((len(q.ax), len(q.z)))
        abs_moment = 0.
        reference_energy_change = 0.
        t0, x0 = (np.linspace(0, 1, n) for n in self.order.shape)
        for v in q.blocks():
            _, dr, dz = self.chart.physical(0, v["r"], v["z"])
            active = (dr != 0) | (dz != 0)
            if not active.any():
                continue
            delta = moving_energy(v["K"][active], v["J"][active], v["hoop"][active],
                                  dr[active], dz[active], step) - v["W"][active]
            moment = delta*v["measure"][active]
            deposit(plane, q.ax, q.z, v["x"][active], v["z"][active], moment)
            abs_moment += abs(moment).sum()
            reference_order = interpolation(t0, x0, self.order, v["z"][active]/q.L, v["x"][active])[0]
            reference_energy_change += moment @ reference_order
        self.cache[step] = (assemble(self.model, plane), float(abs_moment), float(reference_energy_change))
        return self.cache[step]


class MovedCarrier:
    def __init__(self, m, chart, carrier, physical):
        self.m = m
        self.base = m.D, m.G, m.B, m.hq, m.material_load, m.ref
        data = m.D
        self.u = np.zeros((3, len(data)))
        roi = (data[:, 3] > chart.window[0]) & (data[:, 3] < chart.window[-1])
        r, z = data[:, 1], data[:, 3]
        outer = r[roi] > m.shape.radius(z[roi])[0]
        for bank in (0, 1):
            selected = np.flatnonzero(roi)[outer == bank]
            self.u[:, selected] = carrier(bank, r[selected], z[selected])
        self.slit = slit_derivative(m, physical)

    def apply(self, step):
        m = self.m
        D, G, B, _, _, _ = self.base
        eta, ur, uz = self.u
        jac = 1+step*uz
        if jac.min() < .95:
            raise ValueError("Carrier mesh determinant below .95")
        m.G = [(G[0]+G[2].multiply((-step*ur/jac)[:, None])).tocsr(), G[1],
               G[2].multiply((1/jac)[:, None]).tocsr()]
        m.D = D.copy()
        m.D[:, 0] *= jac
        m.D[:, 3] += step*eta
        m.D[:, 9] *= jac
        m.hq = m.D[:, 0]
        m.B = B+step*self.slit
        return float(jac.min())

    def reset(self):
        self.m.D, self.m.G, self.m.B, self.m.hq, self.m.material_load, self.m.ref = self.base


def run(a):
    from mouth_material_order_transfer import MixedOrder
    from mouth_rim_swept_cells import SweptCells

    base, trial, m, _, level, sha, _, _ = from_archive(a.state, a.wall, a.fields, True)
    m.material_load = material_load(m, a.load, a.wall)
    initial_load = m.material_load
    report = json.loads(a.evidence.read_text())
    accepted = next(c for c in report["cases"] if c["inputs"]["grid"] == level)
    if accepted["zero_state_sha256"] != file_sha(a.zero) or accepted["inputs"]["field_sha256"] != sha:
        raise ValueError("Different zero-displacement material or carrier fields")
    if accepted["inputs"]["reference_axes_sha256"] != file_sha(a.axes):
        raise ValueError("Reference axes mismatch")
    with np.load(a.zero, allow_pickle=False) as archive:
        zero = [archive["wound"], archive["rest"]]
    chart = SharedRimChart(m, trial.model, a.axes)
    mixed = MixedOrder(m, base.c)
    crystal = CrystalEnergy(m, chart, a.state, a.wall, base.c)
    zero_stats = []
    for u, ref in zip(zero, accepted["zero"], strict=True):
        e, grad = m.evaluate(u)
        norm = np.max(abs(m.norm(m.vars(grad))))
        if abs(e-ref["energy"]) > 2e-12 or norm >= 5e-9:
            raise ValueError("Archived zero no longer stationary in native material")
        zero_stats.append(dict(energy=e, stationarity=norm))
    info = dict(grid=level, fields=sha, zero_sha=file_sha(a.zero), state_sha=file_sha(a.state),
                load_sha=file_sha(a.load), wall_sha=file_sha(a.wall), axes_sha=file_sha(a.axes),
                window=chart.window.tolist(), tents=chart.tents.tolist(), radial_taper=chart.outer_taper.tolist(),
                edges={n: list(v) for n,v in chart.edges.items()}, step=a.step, tolerance=a.tolerance,
                zero=zero_stats)
    record(a.folder, "input", info)
    rows = []
    for name, n in (("stationary", 0), ("four", 4), ("eight", 8)):
        physical = chart.physical if n else chart.still
        carrier = (lambda bank, r, z: chart.carrier(bank, r, z, n)) if n else chart.stationary
        moving = MovedCarrier(m, chart, carrier, physical)
        swept = SweptCells(m, base.c, chart.tents, physical, carrier)
        for label, h in (("plus", a.step), ("minus", -a.step), ("plus_half", a.step/2),
                         ("minus_half", -a.step/2)):
            jac = moving.apply(h)
            delta, abs_moment, common_crystal = crystal.increment(h) if n else (None, 0., 0.)
            load = initial_load + delta if n else initial_load
            m.material_load = load
            mixed.install(h, physical, carrier, load)
            mixed_blocks = mixed.active
            m.evaluate = swept.attach(m.evaluate, h)
            pair = []
            fields = []
            for rest, u0 in enumerate(zero):
                e0 = m.evaluate(u0)[0]
                u, stats = rematch(m, u0, bool(rest), a.tolerance)
                stats.update(frozen_energy=float(e0), change_max=np.max(abs(u-u0),axis=1).tolist(),
                             crystal_strain_increment=float(delta@u[2]) if n else 0.,
                             swept_layer_energy=float(swept.work(u,h)))
                pair.append(stats); fields.append(u)
            path = a.folder/(name+"_"+label+".npz")
            np.savez_compressed(path, wound=fields[0], rest=fields[1])
            row = dict(case=name, step=h, jacobian=jac, mixed_blocks=mixed_blocks, crystal_abs_energy_moment=abs_moment,
                       common_crystal_strain_energy_change=common_crystal, pair=pair, file=path.name, sha256=file_sha(path))
            record(a.folder, name+"_"+label, row)
            print(name,label,flush=True)
            rows.append(row)
            mixed.reset()
            moving.reset()
    record(a.folder,"result",dict(input=info,signed=rows))


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ("state","wall","fields","load","axes","zero","evidence","folder"):
        parser.add_argument(key,type=Path)
    parser.add_argument("--step",type=float,default=2e-5)
    parser.add_argument("--tolerance",type=float,default=5e-9)
    run(parser.parse_args())
