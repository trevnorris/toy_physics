"""Carrier transit and angular marker spacing through an archived converting wall."""
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq

from mouth_curved_phase import quadratic
from mouth_order_wave_pulse import Geometry, ROOT, digest, emit


RECORD = Path('research/pde_ledger_v3/_measurements/mouth_conversion_contact_transfer.json')


class WallTrace:
    def __init__(self, mesh, multiple):
        record = json.loads(RECORD.read_text())
        folder = ROOT / f'contact_supply_n{mesh}_01'
        setup, calibration = map(json.loads, (folder / 'stdout').read_text().splitlines()[:2])
        self.geometry = Geometry(mesh, setup['profile'], setup['roof'], setup['far'])
        self.edges, radial, cells = self.geometry.order_blocks[2]
        self.nodes = cells[::len(radial)-1, [0, 3, 6]]
        field = Path(setup['source'])
        order = folder / f'order_{multiple:.5f}.npz'
        self.hashes = {str(path): digest(path) for path in (field, order)}
        expected = record['pinned_inputs_sha256'] | record['generated_fields_sha256']
        if any(expected[path] != sha for path, sha in self.hashes.items()):
            raise ValueError('Field differs from the conversion test')
        with np.load(field) as data:
            np.testing.assert_array_equal(data['points'], self.geometry.points)
            self.psi = data['psi'][self.nodes]
        with np.load(order) as data:
            self.chi = data['chi'][self.nodes]
        self.intake = calibration['intake']
        self.head = multiple * calibration['flat_head']

    def sample(self, values, w):
        if np.any(w < self.edges[0]) or np.any(w > self.edges[-1]):
            raise ValueError('Sample outside guide')
        i = np.clip(np.searchsorted(self.edges, w, side='right')-1, 0, len(self.edges)-2)
        width = np.diff(self.edges)[i]
        basis, deriv = quadratic(2*(w-self.edges[i])/width-1)
        return (basis*values[i]).sum(axis=-1), (deriv*values[i]).sum(axis=-1)*2/width

    def crossing(self, level):
        # All crossings in Q2 elements, including any hidden by equal endpoint signs.
        roots = []
        for left, right, (a, b, c) in zip(self.edges[:-1], self.edges[1:], self.chi, strict=True):
            for t in np.roots([2*(a+c-2*b), 4*b-3*a-c, a-level]):
                if abs(t.imag) < 1e-12 and 0 <= t.real < 1:
                    roots.append(left+(right-left)*t.real)
        if len(roots) != 1:
            raise ValueError(f'Expected one crossing at {level}, found {roots}')
        return roots[0]

    def speed(self, w):
        field, derivative = self.sample(self.psi, w)
        slope = self.geometry.radius(w)[1]
        speed = (derivative/field).imag/(1+slope**2)
        if np.any(speed <= 0):
            raise ValueError('Non-inward wall carrier')
        return speed

    def transit(self, start, stop, degree):
        knots = np.r_[start, self.edges[(start < self.edges) & (self.edges < stop)], stop]
        roots, weights = leggauss(degree)
        widths = np.diff(knots)[:, None]
        sample = knots[:-1, None] + widths*(1+roots)/2
        return np.sum(widths*weights/(2*self.speed(sample)))


def main():
    for mesh in (34, 48):
        for multiple in (0., .5):
            path = WallTrace(mesh, multiple)
            start = path.crossing(.95)
            initial_radius = path.geometry.radius(start)[0]
            readings = []
            for level in (.95, .85, .5, .45, .1, .01):
                end = path.crossing(level)
                transit = path.transit(start, end, 8)
                radius = path.geometry.radius(end)[0]
                readings.append(dict(order=level, depth=end, elapsed=transit,
                                     angular_spacing_ratio=radius/initial_radius))
                np.testing.assert_allclose(path.sample(path.chi, end)[0], level, atol=1e-12)
                np.testing.assert_allclose(transit, path.transit(start, end, 4), atol=1e-10, rtol=1e-10)
            emit(mesh=mesh, head_fraction=multiple, head=path.head, intake=path.intake,
                 field_sha256=path.hashes, readings=readings)


if __name__ == '__main__':
    main()
