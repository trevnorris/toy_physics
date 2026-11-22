"""Virtual tangential steps on saved 4D half-cells, with both periodic faces included."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_bulk_isotropic_4d import radial
from mouth_bulk_relaxed_shape import PeriodicMaterial
from mouth_bulk_temperature_readout import guarded
from mouth_bulk_thermal_cell import ThermalCell, library
from neutral_patch_sites import emit


BATCH = Path('_scratch/s11c/thermal87_batch02')
REFERENCE = Path('research/pde_ledger_v3/_measurements/mouth_bulk_temperature_steps.json')
STEP = .0002  # Relative tangential displacement, in neutral bond lengths.
PHASES = ('prepare', 'undriven')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_state(reference, rung, history, phase):
    name = f'r{rung}h{history}_{phase}'
    folder = BATCH / name
    rows = guarded(folder)
    prior = reference['runs'][name]
    if sha(folder/'stdout') != prior['stdout_sha256']:
        raise ValueError(f'Original production log changed: {folder}')
    path = folder/'native/end.npz'
    with np.load(path, allow_pickle=False) as data:
        pos, box = data['pos'].copy(), data['box'].copy()
        metadata = int(data['rung']), int(data['history']), str(data['phase']), float(data['temperature'])
    setup = rows[0]
    if metadata != (rung, history, phase, setup['temperature']) or setup['phase'] != phase:
        raise ValueError(f'State does not match its original log: {path}')
    last = rows[-1]
    if last['frame'] != 48 or len(pos) != 128:
        raise ValueError('Incomplete neutral sample')
    return path, pos, box, setup, last


def virtual_step(lib, pos, box, cut, start_labels, origin):
    inv = np.linalg.inv(box)
    normal = inv[:, 0] / np.linalg.norm(inv[:, 0])
    tangent = box[1] / np.linalg.norm(box[1])
    area = np.linalg.det(box)*np.linalg.norm(inv[:, 0])
    np.testing.assert_allclose(tangent@normal, 0., atol=1e-13)
    np.testing.assert_allclose(area, np.sqrt(np.linalg.det(box[1:]@box[1:].T)), rtol=1e-13)
    y = pos@inv
    coord = y[:, 0]-origin
    spatial = (coord % 1. >= .5).astype(int)
    native = ThermalCell(lib, box, len(pos), cut)
    energy, grad, _ = native.read(pos)
    # Site-by-site winding tests the full native image sum independent of the
    # central-difference displacements. A rigid translation has no work or load.
    winding = np.random.default_rng(714).integers(-3, 4, size=pos.shape)
    rewound = pos+winding@box
    np.testing.assert_array_equal(((rewound@inv)[:, 0]-origin) % 1. >= .5, spatial)
    shift_energy, shift_grad, _ = native.read(pos+.37*tangent)
    rew_energy, rew_grad, _ = native.read(rewound)
    np.testing.assert_allclose([shift_energy, rew_energy], energy, rtol=0, atol=2e-9)
    np.testing.assert_allclose([shift_grad, rew_grad], np.stack([grad,grad]), rtol=0, atol=5e-9)
    rigid_load = np.sum(grad@tangent)/(2*area)
    np.testing.assert_allclose(rigid_load, 0., atol=2e-11)

    # Explicit all-image Python enumeration, independent of the native site
    # gradient. A long bond can cross both faces, with opposite directed jumps.
    i, j, images = PeriodicMaterial(len(pos), cut).pairs(y, box)
    d = (y[j]-y[i]+images)@box
    r = np.linalg.norm(d, axis=1)
    _, first, second = radial(r, cut)
    cosine = (d@tangent)/r
    pair_drag = first*cosine
    pair_stiffness = (second-first/r)*cosine**2+first/r
    u, v = coord[i], coord[j]+images[:, 0]
    crossings = np.stack((np.floor(u)-np.floor(v), np.floor(v-.5)-np.floor(u-.5)))
    np.testing.assert_array_equal(crossings.sum(axis=0), spatial[j]-spatial[i])
    # Group displacement changes by exactly q between the two groups; subtract
    # its mean to prevent any center-of-mass translation when counts differ.
    def probe(labels):
        direction = (labels-labels.mean())[:, None]*tangent
        jump = labels[j]-labels[i]
        baseline = np.sum(grad*direction)/(2*area)
        direct = np.dot(pair_drag, jump)/(2*area)
        exact = np.dot(pair_stiffness, jump**2)/(2*area)
        np.testing.assert_allclose(direct, baseline, atol=2e-10)
        increments, force_slopes, energy_slopes = [], [], []
        for step in (STEP, STEP/2):
            es, inc = [], []
            for offset in (step, -step):
                moved = pos+offset*direction
                np.testing.assert_allclose((moved@inv)[:, 0], y[:, 0], atol=2e-13)
                e, g, _ = native.read(moved)
                es.append(e)
                inc.append(np.sum((g-grad)*direction)/(2*area))
            increments.append(inc)
            slope = (inc[0]-inc[1])/(2*step)
            e_slope = ((es[0]-energy)+(es[1]-energy))/(2*area*step**2)
            derivative = (es[0]-es[1])/(4*area*step)
            np.testing.assert_allclose(derivative, baseline, atol=1.e-4, rtol=2e-5)
            np.testing.assert_allclose([slope,e_slope], exact, atol=.0015, rtol=3e-5)
            force_slopes.append(slope)
            energy_slopes.append(e_slope)
        return dict(baseline=baseline, signed_force_increments=increments, full=force_slopes[0],
                    half=force_slopes[1], analytic=exact, energy_full=energy_slopes[0],
                    energy_half=energy_slopes[1])

    face = probe(spatial)
    face_curvatures = crossings@(pair_stiffness*(spatial[j]-spatial[i]))/area
    np.testing.assert_allclose(face_curvatures.mean(), face['analytic'], atol=2e-11)
    face_preloads = crossings@pair_drag/area
    np.testing.assert_allclose(face_preloads.mean(), face['baseline'], atol=2e-11)
    labels = spatial if start_labels is None else start_labels
    tagged = face if np.array_equal(labels,spatial) else probe(labels)
    return spatial, energy, dict(area_per_face=area, normal=normal.tolist(), tangent=tangent.tolist(),
        origin_fraction=origin, upper_count=int(spatial.sum()), changed_side=int(np.count_nonzero(labels!=spatial)),
        face_preloads=face_preloads.tolist(), face_curvatures=face_curvatures.tolist(), spatial=face,
        original_population=tagged, common_force_per_two_areas=rigid_load,
        max_winding_or_translation_force_error=max(abs(shift_grad-grad).max(),abs(rew_grad-grad).max()))


def run(work):
    reference = json.loads(REFERENCE.read_text())
    lib = library(work)
    for rung in (0,4):
        for history in range(3):
            original = None
            for phase in PHASES:
                path, pos, box, setup, last = load_state(reference, rung, history, phase)
                fractions = pos@np.linalg.inv(box)
                if phase == 'prepare':
                    start, original_hash = fractions.copy(), sha(path)
                elif original_hash != setup['source_sha256']:
                    raise ValueError('Undriven endpoint is from a different starting state')
                # Frames remove affine cell motion and rigid center drift using
                # the original trajectory's unwrapped coordinates.
                origin = float(np.mean(fractions[:,0]-start[:,0]))
                spatial, energy, result = virtual_step(lib,pos,box,setup['cutoff'],original,origin)
                np.testing.assert_allclose(energy/len(pos), last['readings'][0]['energy_per_site'], atol=2e-10)
                if phase == 'prepare':
                    original = spatial
                emit(test='instantaneous_cut', temperature=setup['temperature'], history=history,
                     phase=phase, state=str(path), sha256=sha(path), measured_density=len(pos)/np.linalg.det(box),
                     original_neighbors=last['initial'], missing_neighbors=last['original_missing'], **result)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work', required=True, type=Path)
    run(parser.parse_args().work)
