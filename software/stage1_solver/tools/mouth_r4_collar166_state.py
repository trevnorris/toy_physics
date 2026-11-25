"""Frozen and field-rematched R4 readouts for the step-166 collar reference comparison."""
import argparse
from pathlib import Path
import sys

import numpy as np

from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
from mouth_continuous_depth_mesh import axes_for_depth
from mouth_continuous_shape_relax import ShapeTrial
from mouth_r4_aligned_crystal import AlignedCrystal, held_forces
from mouth_r4_aligned_rematch import FIELD_TOL, FIELDS, STRENGTH, front
from mouth_r4_eulerian_release import order_grid
from mouth_r4_eulerian_solid import SpatialCrystal
from mouth_r4_segmented_held import sign_sweep
from mouth_r4_slaved_melt import FLOW
from mouth_r4_smooth_fields import GaugeGeometry, coherent, material_load, read_fields
from mouth_r4_smooth_held import SmoothHeld
from mouth_r4_stage2_outcome import verify_guard
from mouth_whole_r4_relax import atomic_json, save_fields


SOURCE = '4b84a0d2ad53889a1d619e68960df70233ceb589fc4fe6c336c6e59e832a2694'


def load_held(study, output, epsilon, smooth=True, *, solve_fields=False):
    if epsilon not in STRENGTH.values():
        raise ValueError('Expected a target held gauge strength')
    source = study.parent / 'r4_segmented_held01/primary/step_00166.npz'
    record = read(source.with_suffix('.json'))
    if record['outer'] != 166 or not digest(source) == SOURCE == record['sha256']:
        raise ValueError('Step-166 source digest differs')
    driver = SmoothHeld(study.parent, output, epsilon)
    with np.load(source, allow_pickle=False) as state:
        if (int(state['iteration']) != record['iteration'] or int(state['radius']) != 4 or int(state['winding']) != 1
                or str(state['parent_sha256']) != CHECKPOINT_SHA or not np.array_equal(driver.x, state['x'])):
            raise ValueError('Changed held shape or preparation')
        driver.labels, driver.fields = state['labels'], state['field']
    trial = ShapeTrial(driver.base, driver.prep, driver.x, 8, 24, 6, 3, 2,
                       precondition=solve_fields,
                       axes=axes_for_depth(driver.shape(driver.x).length, driver.reference_depth))
    trial.model.extended_energy_sum = True
    driver.current_trial = trial
    trial.smooth_order_geometry = GaugeGeometry(trial, aligned=True)
    if smooth:
        from mouth_r4_smooth_collar_reference import enable_smooth_collar
        enable_smooth_collar(trial.b)
    if not np.array_equal(driver.fields[2, trial.first][trial.number], driver.fields[2]):
        raise ValueError('Saved order is discontinuous at the bonded faces')
    return driver, trial, SpatialCrystal(trial.b), order_grid(trial, driver.fields)


def readout(driver, trial, solid):
    labels, fields = driver.labels, driver.fields
    joint = AlignedCrystal(driver, trial, solid, labels, order_grid(trial, fields))
    _, _, material = joint.evaluate(labels)
    if trial.model.precondition:
        physical, load, field = read_fields(driver, trial, solid, labels, fields)
        if abs(material['gauge_energy'] - field['gauge_energy']) > 1e-15:
            raise ArithmeticError('Field and label gauges disagree')
    else:
        load, exchange = material_load(trial, solid, labels, fields)
        trial.model.material_load = load
        physical, _ = trial.model.evaluate(fields)
        field = dict(shared_work_error=exchange)
    measured = trial.measure(fields, load)
    return dict(total_energy=float(physical + material['gauge_energy']), physical_energy=float(physical),
                crystal_energy=material['elastic_energy'], gauge_energy=material['gauge_energy'],
                peak_active=material['max_free_label_force'], peak_node=material['worst_free_active'],
                field=field, chart=material['chart'], min_Jinv=material['min_Jinv_global'],
                front=front(trial, fields), current={k: float(measured[k]) for k in FLOW})


def run(study, output, case, operation):
    output.mkdir(exist_ok=False, parents=True)
    eps = STRENGTH[case]
    smooth = operation != 'legacy'
    driver, trial, solid, _ = load_held(study, output, eps, smooth, solve_fields=operation=='rematch')
    result = dict(command=sys.orig_argv, source_sha256=SOURCE, script_sha256=digest(Path(__file__)),
                  radius=4, outer=166, epsilon=eps, operation=operation, frozen=readout(driver, trial, solid))
    if operation == 'legacy':
        forces = held_forces(driver, trial, solid, driver.labels, driver.fields, output)
        result['six'] = forces['six_total']
        result['modes'] = [row['mode'] for row in forces['rows']]
    elif operation == 'rematch':
        driver.labels.flags.writeable = False
        driver.fields, _, _, history = coherent(driver, trial, solid, driver.labels, driver.fields,
                                                field_target=FIELD_TOL, raw_target=FIELD_TOL)
        result['rematched'] = matched = readout(driver, trial, solid)
        if not all(np.isfinite(matched['field'][k]) and matched['field'][k] <= FIELD_TOL for k in FIELDS):
            raise ArithmeticError(f'Field rematch incomplete: {matched["field"]}')
        result['field_work'] = history['order_step_energy']
        state = output / 'rematched.npz'
        save_fields(state, labels=driver.labels, field=driver.fields, x=driver.x, radius=4, winding=1,
                    epsilon=eps, outer=166, parent_sha256=CHECKPOINT_SHA, initial_sha256=SOURCE)
        result['state_sha256'] = digest(state)
        result['sign'] = sign_sweep(driver, trial, output / 'order_sign.json', matched)
    result['passed'] = True
    atomic_json(output / 'result.json', result)
    print(operation, case, result['frozen']['crystal_energy'], result.get('six', ''), flush=True)


def static_comparison(study, output):
    cases = {}
    for case, epsilon in STRENGTH.items():
        original = (study / 'legacy_primary_static' if case == 'primary'
                    else study / 'frozen_controls/legacy_secondary')
        smooth = study / 'frozen_controls' / ('static_' + case)
        signed = study / 'derivatives_lean' / (case + '_shape')
        for folder in (original, smooth, signed):
            verify_guard(folder.with_name(folder.name + '_guard'), 0)
        c0, c2, work = (read(folder / 'result.json') for folder in (original, smooth, signed))
        if (not all(row['passed'] and row['epsilon'] == epsilon for row in (c0, c2, work))
                or not c0['source_sha256'] == c2['source_sha256'] == work['state_sha256'] == SOURCE
                or c0['modes'] != [r['mode'] for r in work['rows']]):
            raise ValueError('Inconsistent fixed-label readout')
        a, b = c0['frozen'], c2['frozen']
        keys = ('total_energy', 'physical_energy', 'crystal_energy', 'gauge_energy', 'peak_active')
        changes = {key: dict(linear=a[key], cubic=b[key], delta=b[key]-a[key]) for key in keys}
        currents = {k: b['current'][k]-a['current'][k] for k in FLOW}
        forces = [dict(mode=name, linear=old, cubic=new['analytic'], delta=new['analytic']-old,
                       central=new['central']) for name, old, new in zip(c0['modes'], c0['six'], work['rows'])]
        if any(currents.values()) or changes['gauge_energy']['delta'] or a['front'] != b['front']:
            raise ArithmeticError('Fixed-label currents, gauge or order front changed')
        cases[case] = dict(epsilon=epsilon, scalars=changes, shape_forces=forces,
                           shape_tolerance=work['threshold'], flow=b['current'], flow_delta=currents,
                           front=b['front'], chart=b['chart'], min_sampled_Jinv=b['min_Jinv'],
                           old_peak=a['peak_node'], new_peak=b['peak_node'],
                           evidence={str(d.relative_to(study)):digest(d/'result.json') for d in
                                     (original, smooth, signed)})
    atomic_json(output, dict(command=sys.orig_argv, source_commit='ec9fa56b', state_sha256=SOURCE,
                             radius=4, step=166, labels_and_coherent_fields_held=True, cases=cases))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('case', choices=[*STRENGTH, 'both'])
    parser.add_argument('operation', choices=['legacy', 'static', 'rematch', 'compare'])
    opts = parser.parse_args()
    if opts.operation == 'compare' and opts.case == 'both':
        static_comparison(opts.study.resolve(), opts.output.resolve())
    elif opts.case in STRENGTH and opts.operation != 'compare':
        run(opts.study.resolve(), opts.output.resolve(), opts.case, opts.operation)
    else:
        parser.error('Use both with compare, or one target gauge with legacy, static or rematch')
