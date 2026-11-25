"""Rematch the native R4 order and carriers at fixed step-108 labels and shape."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


FIELD_TOL = 1e-8
FIELDS = ('order_projected', 'order_kkt', 'amplitude_raw', 'phase_raw', 'carrier_metric')
FORMER = (36353, 45532, 62811)
STRENGTH = {'primary': 1e-8, 'secondary': 1e-10}


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def pin_sources(study):
    """The original evaluator supplies the unchanged mesh, solver and crystal kernels."""
    from mouth_r4_newton_stop import PIN_SHA
    prior = study.parent/'r4_newton_held01'
    pin = prior/'pin/software/stage1_solver/tools'
    manifest = prior/'pinned.json'
    if digest(manifest) != PIN_SHA or any(digest(pin/name) != key
            for name, key in json.loads(manifest.read_bytes()).items()):
        raise ValueError('Step-108 evaluator pin has changed')
    overlay = Path(__file__).resolve().parent
    sys.path.insert(0, str(pin))
    # Only the new checkpoint reader, integration rule and field evaluator
    # override the archived implementation in the narrow run-specific overlay.
    sys.path.insert(0, str(overlay))
    return prior


def field_state(study, case):
    from mouth_r4_newton_stop import read_state
    prior = pin_sources(study)
    from mouth_exact_checkpoint_control import CHECKPOINT_SHA
    from mouth_r4_smooth_fields import GaugeGeometry
    from mouth_r4_smooth_held import SmoothHeld
    from mouth_r4_eulerian_solid import SpatialCrystal

    folder = study/case
    driver = SmoothHeld(study.parent, folder, 1e-8)
    record, fields, labels = read_state(prior/'primary', 'last_active', driver, CHECKPOINT_SHA)
    driver.epsilon = STRENGTH[case]
    trial = driver.trial()
    trial.smooth_order_geometry = GaugeGeometry(trial, aligned=True)
    return driver, trial, SpatialCrystal(trial.b), record, fields, labels


def front(trial, fields):
    from mouth_r4_eulerian_release import order_grid
    from mouth_r4_fold_order import all_low_after
    return all_low_after(trial.model.z[1], order_grid(trial, fields).transpose(0, 2, 1), 0.)


def material_nodes(trial, joined):
    sites = trial.number[trial.b.ids.transpose(0, 2, 1)]
    return np.argwhere(sites == joined).tolist()


def save_json(path, data):
    with path.open('x') as stream:
        json.dump(data, stream, allow_nan=False, indent=2)
        stream.write('\n')


def rematch(study, case):
    driver, trial, solid, record, fields, labels = field_state(study, case)
    from mouth_exact_checkpoint_control import CHECKPOINT_SHA
    from mouth_exchange_carrier import BIAS
    from mouth_r4_newton_stop import PIN_SHA
    from mouth_r4_slaved_melt import FLOW
    from mouth_r4_smooth_fields import coherent, read_fields

    folder = study/case
    path = folder/'rematched.npz'
    if path.exists() or path.with_suffix('.json').exists():
        raise FileExistsError(path)
    labels.flags.writeable = False
    fields.flags.writeable = False
    before, _, initial = read_fields(driver, trial, solid, labels, fields)
    u, _, _, history = coherent(driver, trial, solid, labels, fields,
                                field_target=FIELD_TOL, raw_target=FIELD_TOL)
    physical, load, final = read_fields(driver, trial, solid, labels, u)
    start, end = (v[2, trial.first] for v in (fields, u))
    if (not np.isfinite(u).all() or np.any((end < 0) | (end > 1))
            or not np.array_equal(u[2], end[trial.number])
            or np.any(u[1, trial.model.feed] != 0) or np.any(u[1, trial.model.bulk] != BIAS)):
        raise ArithmeticError('The field solution violates the bounded or registered field domain')
    current = trial.measure(u, load)
    new = np.flatnonzero((start == 0) & (end > 0) & trial.free)
    lost = np.flatnonzero((start > 0) & (end == 0) & trial.free)
    def describe(i):
        return dict(node=int(i), before=float(start[i]), after=float(end[i]),
                    physical_depth=float(trial.depth[i]), material_grid=material_nodes(trial, i))
    with path.open('xb') as stream:
        np.savez_compressed(stream, field=u, labels=labels, x=driver.x, radius=4, winding=1,
                            epsilon=driver.epsilon, iteration=driver.iteration, outer=108,
                            parent_sha256=CHECKPOINT_SHA, source_sha256=record['sha256'],
                            stencil_points=39)
    state_digest = digest(path)
    with np.load(path, allow_pickle=False) as saved:
        if (not np.array_equal(saved['field'], u) or not np.array_equal(saved['labels'], labels)
                or not np.array_equal(saved['x'], driver.x) or float(saved['epsilon']) != driver.epsilon):
            raise ArithmeticError('Stored held state failed its field, label, or geometry round-trip')
    frozen = json.loads((study.parents[2]/'research/pde_ledger_v3/_measurements/mouth_r4_edge_static.json').read_bytes())
    frozen_gauge = frozen['comparisons'][str(driver.epsilon)]['new_gauge']
    result = dict(file=path.name, sha256=state_digest, source_sha256=record['sha256'], outer=108,
                  source_commit='0afc9e1e', evaluator_manifest_sha256=PIN_SHA,
                  sources={p.name:digest(p) for p in map(Path, (__file__,
                      sys.modules['mouth_r4_smooth_fields'].__file__, sys.modules['mouth_r4_aligned_stencil'].__file__))},
                  epsilon=driver.epsilon, stencil_points=39, fixed_radius=4, winding=1,
                  phase_drop=BIAS, x=driver.x.tolist(), initial=initial, final=final,
                  field_steps=history['field_iterations'], work=history['order_step_energy'],
                  before_total=before+initial['gauge_energy'], total=physical+final['gauge_energy'],
                  initial_gauge_minus_frozen=initial['gauge_energy']-frozen_gauge,
                  former=[describe(i) for i in FORMER], newly_ordered=[describe(i) for i in new],
                  newly_zero=[describe(i) for i in lost], front_before=front(trial, fields),
                  front_after=front(trial, u), current={k:float(current[k]) for k in FLOW},
                  max_field_change=dict(zip(('amplitude','phase','order'),
                      np.max(np.abs(u-fields), axis=1).tolist())), tolerance=FIELD_TOL)
    ends = [final[k] for k in FIELDS]
    result['field_matched'] = (all(np.isfinite(ends)) and max(ends) <= FIELD_TOL
                                and abs(result['initial_gauge_minus_frozen']) <= 2e-19)
    save_json(path.with_suffix('.json'), result)
    print(json.dumps(result, allow_nan=False), flush=True)
    if not result['field_matched']:
        raise ArithmeticError('A held field residual or frozen gauge control failed')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study', type=Path)
    parser.add_argument('case', choices=STRENGTH)
    args = parser.parse_args()
    rematch(args.study.resolve(), args.case)
