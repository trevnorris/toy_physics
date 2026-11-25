"""Measure the preserved step-108 R4 Newton state using its frozen evaluator."""
import argparse
from fractions import Fraction
import hashlib
from itertools import product
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import mpmath as mp
import numpy as np


PIN_SHA = '22a06d53f4b64bd6d6be19fc2c8c2f8d038909ebd5f7357c9a327941d386820b'
STATES = {
    'last_active': (108, 'a2210361aa5b38712de7db618b5220985e12f39f3a089351dfd3767bef0a7178'),
    'rematched': (100, '63a2de2ae45c6c2b197573ac70b2f035142f538d59908bb068ceda87b0570542'),
}


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read_state(folder, name, driver, parent):
    step, expected = STATES[name]
    path = folder / (name + '.npz')
    record = json.loads(path.with_suffix('.json').read_bytes())
    if (record['outer'], record['file'], record['sha256'], sha(path)) != (step, path.name, expected, expected):
        raise ValueError(f'Saved {name} step or hash changed')
    with np.load(path, allow_pickle=False) as saved:
        if ((int(saved['radius']), int(saved['winding']), int(saved['iteration']),
             float(saved['epsilon']), str(saved['parent_sha256'])) != (4, 1, driver.iteration, driver.epsilon, parent)
                or not np.array_equal(saved['x'], driver.x)):
            raise ValueError(f'Unexpected {name} geometry, gauge or checkpoint ancestry')
        return record, saved['field'], saved['labels']


def exact_det(axes, labels, cell, fractions):
    """Determinant of the trilinear chart at a cell-local point, in exact binary-input arithmetic."""
    widths = [Fraction(float(a[i+1])) - Fraction(float(a[i])) for a, i in zip(axes, cell)]
    local = list(map(Fraction, fractions))
    gradient = [[Fraction(0) for _ in range(3)] for _ in range(3)]
    for corner in product(range(2), repeat=3):
        node = tuple(i+j for i, j in zip(cell, corner))
        values = [Fraction(float(v)) for v in labels[node]]
        basis = [f if side else 1-f for f, side in zip(local, corner)]
        for j in range(3):
            derivative = (2*corner[j]-1)*basis[(j+1)%3]*basis[(j+2)%3]/widths[j]
            for i in range(3):
                gradient[i][j] += values[i]*derivative
    a, b, c = gradient
    return sum(a[i]*(b[(i+1)%3]*c[(i+2)%3]-b[(i+2)%3]*c[(i+1)%3]) for i in range(3))


def readback(study):
    manifest = study/'pinned.json'
    pin = study/'pin/software/stage1_solver/tools'
    if sha(manifest) != PIN_SHA or any(sha(pin/k) != v for k, v in json.loads(manifest.read_bytes()).items()):
        raise ValueError('Primary evaluator source has changed')
    sys.path.insert(0, str(pin))
    from mouth_exact_checkpoint_control import CHECKPOINT_SHA
    from mouth_exchange_carrier import BIAS
    from mouth_joint_endpoint_audit import geometry
    from mouth_r4_eulerian_corner_audit import point
    from mouth_r4_eulerian_release import order_grid
    from mouth_r4_eulerian_solid import SpatialCrystal
    from mouth_r4_fold_order import all_low_after
    from mouth_r4_newton_held import build_joint
    from mouth_r4_refinement_review import _collar_reference
    from mouth_r4_slaved_melt import FLOW
    from mouth_r4_smooth_fields import read_fields
    from mouth_r4_smooth_held import SmoothHeld

    folder = study/'primary'
    issue = json.loads((folder/'ISSUE.json').read_bytes())
    rejected = issue['witness']['tests'][-1]
    cell, local = rejected['site']['cell'], rejected['site']['chart_fractions']
    if (issue['error'] != 'InfeasibleLabels' or cell != [25, 11, 10] or local != [0, 0, 0.3984375]
            or rejected['site']['chart_det'] >= 0):
        raise ValueError('Unexpected terminal witness')
    with TemporaryDirectory(dir=study) as work:
        driver = SmoothHeld(study.parent, Path(work), 1e-8)
        active, driver.fields, driver.labels = read_state(folder, 'last_active', driver, CHECKPOINT_SHA)
        rematch, old_fields, old_labels = read_state(folder, 'rematched', driver, CHECKPOINT_SHA)
        if not np.array_equal(old_fields, driver.fields) or np.array_equal(old_labels, driver.labels):
            raise ValueError('Expected changed crystal labels and fixed step-100 fields')
        trial = driver.trial()
        driver.current_trial = trial
        shape = geometry(trial.shape, trial.polar_shift)
        if any(record[k] != v for record in (active, rematch) for k, v in shape.items()):
            raise ValueError('Reconstructed checkpoint geometry differs')
        if (np.any(driver.fields[1,trial.model.feed] != 0) or np.any(driver.fields[1,trial.model.bulk] != BIAS)
                or not np.array_equal(driver.fields[2,trial.first][trial.number], driver.fields[2])):
            raise ValueError('Reservoir phase or order registry differs')

        solid = SpatialCrystal(trial.b)
        joint, (_, _, material) = build_joint(driver, trial, solid, registry_update=False)
        physical, load, field = read_fields(driver, trial, solid, driver.labels, driver.fields)
        total = physical + material['gauge_energy']
        flows = trial.measure(driver.fields, load)
        order = order_grid(trial, driver.fields)
        front = all_low_after(trial.model.z[1], order.transpose(0,2,1), 0.)
        axes = joint.opt.axes
        det = exact_det(axes, driver.labels, cell, local)
        with mp.workdps(80):
            accepted = point(trial.b, axes, driver.labels, order, cell, local, _collar_reference)
        corners = order[tuple(slice(i, i+2) for i in cell)].ravel().tolist()
        bound = material['chart']['chart_det_global_lower']
        if (float(det) != accepted['chart_det'] or float(det) < bound
                or material['chart']['chart_lowest_cell'] != cell):
            raise ArithmeticError('Independent trilinear determinant or production enclosure disagrees')
        keys = ('elastic_energy', 'gauge_energy', 'min_Jinv_global', 'min_Jinv_active',
                'max_free_label_force', 'gradient_l2', 'max_unit_melt_residual')
        differences = {key: material[key]-active['material'][key] for key in keys}
        differences['chart_det_global_lower'] = bound-active['material']['chart']['chart_det_global_lower']
        if (any(abs(v) > 1e-11*max(1., abs(active['material'][k])) for k,v in differences.items() if k in keys)
                or abs(differences['chart_det_global_lower']) > 1e-13 or abs(total-active['energy']) > 64*np.spacing(total)
                or abs(field['gauge_energy']-material['gauge_energy']) > 1e-12*material['gauge_energy']):
            raise ArithmeticError('Frozen readback disagrees with accepted state')
        return dict(accepted_step=active['outer'], last_field_rematch=rematch['outer'],
            state_sha256={k: v[1] for k, v in STATES.items()}, source_manifest_sha256=PIN_SHA,
            parent_sha256=CHECKPOINT_SHA, radius=4, winding=1, epsilon=driver.epsilon, phase_drop=BIAS,
            shape=driver.x.tolist(), depth=shape['depth'], exit_radius=shape['exit_radius'],
            physical_fields_equal_last_rematch=True, total_energy=total, physical_energy=physical,
            material={k:material[k] for k in (*keys, 'chart', 'stationary', 'worst_free_active', 'worst_global', 'worst_active')},
            total_minus_logged=total-active['energy'], material_minus_logged=differences,
            field=field, current={k:float(flows[k]) for k in FLOW}, sampled_full_melt_front=front,
            limiting_cell=dict(cell=cell, local=local, chart=accepted['physical_chart'],
                physical=accepted['physical'], chi=accepted['chi'], vertex_chi=corners,
                accepted_chart_det_exact=str(det), accepted_chart_det=float(det), accepted_physical_Jinv=accepted['Jinv'],
                accepted_det_minus_global_lower=float(det-Fraction(bound)), rejected_trial=rejected,
                accepted_minus_rejected_chart_det=float(det-Fraction(rejected['site']['chart_det']))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = readback(args.study.resolve())
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(result, allow_nan=False))
