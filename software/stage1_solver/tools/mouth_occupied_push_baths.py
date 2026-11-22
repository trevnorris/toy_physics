"""Extend both stationary reservoirs while preserving the original mouth and Q2 mesh."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from mouth_curved_phase import RATE
from mouth_curved_phase_stability import half_plane_bound
from mouth_flared_receiver import EXIT, EXIT_RADIUS, FlaredDomain, emit
from mouth_free_contact import ARCHIVE

BASE_FAR = 9.0
BASE_ROOF = 10.0


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def extended(edges, end):
    """Retain all old cells and add cells no wider than the final old one."""
    if end < edges[-1]:
        raise ValueError("The expanded boundary lies inside the archived mesh")
    if end == edges[-1]:
        return edges
    count = int(np.ceil((end - edges[-1]) / (edges[-1] - edges[-2])))
    return np.r_[edges, np.linspace(edges[-1], end, count + 1)[1:]]


class NestedBath(FlaredDomain):
    def assemble(self, grids):
        angular, stretch = grids[0]
        axial = grids[2][0]
        old_feed = self.rc + (BASE_FAR - self.rc) * stretch
        old_roof = EXIT + (BASE_ROOF - EXIT) * stretch
        old_outer = EXIT_RADIUS + (BASE_FAR - EXIT_RADIUS) * stretch
        feed = extended(old_feed, self.far)
        roof = extended(old_roof, self.roof)
        outer_radius = extended(old_outer, self.far)
        # Block four's chart labels radii by a fraction of the enlarged bath.
        outer_fraction = (outer_radius - EXIT_RADIUS) / (self.far - EXIT_RADIUS)
        super().assemble([(angular, stretch), (feed, stretch), (axial, stretch),
                          (roof, stretch), (roof, outer_fraction)])


def load_bath(folder):
    outcome = json.loads((folder / 'outcome.json').read_text())
    if (outcome['exitCode'] or not outcome['limitsVerified'] or
            outcome['childOutcome']['exitCode'] or outcome['childOutcome']['guardReason']):
        raise ValueError(f'Unsuccessful stationary bath preparation at {folder}')
    path = folder / 'field.npz'
    record, = [json.loads(line) for line in (folder / 'stdout').read_text().splitlines()]
    if record['field_sha256'] != digest(path):
        raise ValueError(f'Changed stationary bath preparation {path}')
    dom = NestedBath(record['n'], record['profile'], record['roof'], record['far'])
    with np.load(path) as data:
        np.testing.assert_array_equal(data['points'], dom.points)
        return dom, data['psi']


def build(args):
    reference = Path(f'_scratch/s11c/mouth_flared_cubic_n{args.n}') / 'field.npz'
    reference_hash = digest(reference)
    if reference_hash != json.loads(ARCHIVE.read_text())['output_field_sha256'][str(reference)]:
        raise ValueError('The baseline mesh reference changed')
    dom = NestedBath(args.n, 'cubic', args.roof, args.far)
    with np.load(reference) as data:
        baseline_points = data['points']
    distance, point_ids = cKDTree(dom.points).query(baseline_points)
    if distance.max() > 2e-11 or len(np.unique(point_ids)) != len(baseline_points):
        raise ValueError('The enlarged mesh lost a baseline point')
    field, state = dom.solve(RATE)
    certificate = half_plane_bound(dom, field, state['bias'])
    np.testing.assert_allclose(field[dom.feed], 1., atol=1e-14, rtol=0.)
    np.testing.assert_allclose(field[dom.bulk], np.exp(1j * state['bias']), atol=1e-14, rtol=0.)
    sections = [dom.section(field, z) for z in (.638, .993, 1.417, 2.083, 2.793, 2.97)]
    if (certificate['half_plane_lower'] <= 0 or abs(state['feed'] - RATE) > 3e-9 or
            abs(state['bulk'] - RATE) > 3e-9 or
            max(abs(s['integrated_jw'] - RATE) for s in sections) > .0025):
        raise ArithmeticError('Stationary intake and branch verification failed')
    np.savez_compressed(args.output, points=dom.points, psi=field)
    emit(test='nested_bath', n=args.n, profile=dom.profile, far=dom.far, roof=dom.roof,
         baseline_far=BASE_FAR, baseline_roof=BASE_ROOF, reference_sha256=reference_hash,
         preserved_baseline_nodes=len(baseline_points), nodes=len(dom.points),
         max_old_node_distance=float(distance.max()), field_sha256=digest(args.output),
         sections=sections, **state, **certificate)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, required=True)
    parser.add_argument('--far', type=float, required=True)
    parser.add_argument('--roof', type=float, required=True)
    parser.add_argument('--output', type=Path, required=True)
    build(parser.parse_args())
