"""Measure one-sided physical order work on every free zero-order node of step 108.

The absolute energy here uses the model's default summation. Each signed result
is a local long-double energy change. The companion static stencil audit uses
the archived production summation for the absolute-energy invariance check.
"""
import argparse
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np

from mouth_r4_newton_stop import PIN_SHA, read_state, sha


def evaluate(study, output):
    manifest = study / 'pinned.json'
    pin = study / 'pin/software/stage1_solver/tools'
    if sha(manifest) != PIN_SHA or any(sha(pin / k) != v for k, v in json.loads(manifest.read_bytes()).items()):
        raise ValueError('Changed source for the held run')
    sys.path.insert(0, str(pin))
    from mouth_continuous_depth_mesh import axes_for_depth
    from mouth_continuous_shape_relax import ShapeTrial
    from mouth_exact_checkpoint_control import CHECKPOINT_SHA
    from mouth_r4_eulerian_solid import SpatialCrystal
    from mouth_r4_smooth_fields import material_load
    from mouth_r4_smooth_held import SmoothHeld
    from mouth_r4_smooth_order_boundary import PhysicalOrderEnergy, gauge_coefficients

    with TemporaryDirectory(dir=output) as work:
        driver = SmoothHeld(study.parent, Path(work), 1e-8)
        record, fields, labels = read_state(study / 'primary', 'last_active', driver, CHECKPOINT_SHA)
        shape = driver.shape(driver.x)
        trial = ShapeTrial(driver.base, driver.prep, driver.x, 8, 24, 6, 3, 2, precondition=False,
                           axes=axes_for_depth(shape.length, driver.reference_depth))
        trial.model.material_load, _ = material_load(trial, SpatialCrystal(trial.b), labels, fields)
        energy, gradient = trial.model.evaluate(fields)
        native_gradient = trial.pull.T @ gradient[2]
        chi = fields[2, trial.first]
        if not np.array_equal(chi[trial.number], fields[2]):
            raise ArithmeticError('Broken order registry')
        ids = np.flatnonzero((chi == 0) & trial.free)
        oracle = PhysicalOrderEnergy(trial, fields)
        increment = np.longdouble(1e-9)
        work = np.empty(len(ids), dtype=np.longdouble)
        slope = np.empty_like(work)
        for i, node in enumerate(ids):
            local = oracle.at(node)
            slope[i] = local(0)
            work[i] = increment * local(increment)
        mismatch = np.max(np.abs(slope - native_gradient[ids]))
        if not np.all(np.isfinite(work)) or mismatch > 5e-10:
            raise ArithmeticError(f'Physical work failed its slope control: {mismatch}')
        coeff, sites, _ = gauge_coefficients(trial, labels, fields)
        a, b, c = coeff[:, ids]
        unit = increment * (a + increment * (b + increment*c))
        if np.any(unit > 0):
            raise ArithmeticError('Increasing bounded order increased a positive weighted gauge')
        bad = np.flatnonzero(work <= 0)
        node_locations = {}
        for idx in bad:
            i = ids[idx]
            node_locations[str(i)] = dict(material_grid=np.argwhere(sites == i).tolist(),
                                          work=float(work[idx]), derivative=float(slope[idx]),
                                          old_stencil_unit_gauge_work=float(unit[idx]))
        pos = work > 0
        old = {}
        for epsilon in (1e-8, 1e-10):
            gauge = np.longdouble(epsilon) * unit
            total = work + gauge
            old[str(epsilon)] = dict(negative=int(np.sum(total<0)), zero=int(np.sum(total==0)),
                                     min_total=float(total.min()),
                                     worst_node=int(ids[np.argmin(total)]),
                                     max_abs_gauge_over_positive_physical=float(np.max(np.abs(gauge[pos])/work[pos])))
        record_out = dict(state_sha256=record['sha256'], source_manifest_sha256=PIN_SHA,
                          frozen_physical_energy=float(energy), increment=float(increment),
                          free_zero=len(ids), positive=int(np.sum(pos)), negative=int(np.sum(work<0)),
                          zero=int(np.sum(work==0)), min_physical_work=float(work.min()),
                          worst_node=int(ids[np.argmin(work)]), derivative_max_error=float(mismatch),
                          nonpositive_nodes=node_locations, vertex_gauge=old)
        with (output / 'physical.npz').open('xb') as f:
            np.savez_compressed(f, ids=ids, physical=work, slope=slope, old_unit_gauge=unit,
                                sha=record['sha256'])
        with (output / 'physical.json').open('x') as stream:
            json.dump(record_out, stream, indent=2, allow_nan=False)
            stream.write('\n')
        print(json.dumps(record_out, allow_nan=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    evaluate(args.study.resolve(), args.output.resolve())
