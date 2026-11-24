"""Compare independent fixed-radius held equilibria from the corner projection test."""
import argparse
from pathlib import Path

import numpy as np

from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
from mouth_r4_eulerian_relax_labels import FORCE_TOL, absmax
from mouth_r4_stage2_outcome import verify_guard
from mouth_whole_r4_relax import FIELD_TOL, atomic_json


def load_held(study, name, epsilon):
    case = study / name
    verify_guard(study / (name + '_guard'), 0)
    end = read(case / 'HELD.json')
    state = read(case / 'checkpoint_00010.json')
    data = case / state['file']
    if not (end['saved'] == state['sha256'] == digest(data)
            and state['epsilon'] == end['epsilon'] == epsilon
            and state['phase_drop'] == .4510084282857874):
        raise ValueError(f'Unmatched endpoint or gauge at {name}')
    with np.load(data) as saved:
        if (int(saved['radius']), int(saved['winding']), str(saved['parent_sha256'])) != (4, 1, CHECKPOINT_SHA):
            raise ValueError('Wrong material or geometry')
        shape = saved['x']
        order = saved['field'][2]
    m, f, current = state['material'], state['field'], state['current']
    signed, half = (np.array(state[k]) for k in ('shape_force', 'half_shape_forces'))
    checks = dict(
        all_six=signed.shape == half.shape == (6,) and absmax(signed-half) <= 1.e-4,
        crystal=m['stationary'] and max(m['max_free_label_force'], m['max_free_elastic_force']) <= FORCE_TOL,
        gauge=m['max_unit_melt_residual'] <= 3.e-8 and m['max_dead_elastic_force'] <= 1.e-12,
        orientation=min(m['min_Jinv_active'], m['min_Jinv_global'], m['chart']['chart_det_global_lower']) > 0.,
        fields=max(f[k] for k in ('carrier_metric', 'order_kkt', 'order_projected')) < FIELD_TOL
            and max(f[k] for k in ('amplitude_raw', 'phase_raw')) <= FORCE_TOL,
        healthy_steps=state['smallest_accepted_active_fraction'] is None
            or state['smallest_accepted_active_fraction'] > 1.e-4,
        currents=max(abs(current[k]) for k in ('port_error', 'disk_minus_feed_current')) <= 3.e-5)
    if not all(checks.values()):
        raise ValueError(f'Unqualified held state at {name}: {checks}')
    return state, shape, order


def compare(study):
    a, xa, ca = load_held(study, 'primary', 1.e-8)
    b, xb, cb = load_held(study, 'secondary', 1.e-10)
    if not np.array_equal(xa, xb):
        raise ArithmeticError('The paired states have different held shapes')
    tests = {}

    def test(key, delta, limit):
        tests[key] = dict(signed_primary_minus_secondary=delta, absolute_limit=limit,
                          passed=bool(delta is not None and np.max(np.abs(delta)) <= limit))

    test('total_energy', a['energy']-b['energy'], 1.e-6)
    test('physical_energy', (a['energy']-a['material']['gauge_energy'])
                            -(b['energy']-b['material']['gauge_energy']), 1.e-6)
    test('crystal_energy', a['material']['elastic_energy']-b['material']['elastic_energy'], 1.e-6)
    test('force', (np.array(a['shape_force'])-b['shape_force']).tolist(), 1.e-4)
    test('half_step_force', (np.array(a['half_shape_forces'])-b['half_shape_forces']).tolist(), 1.e-4)
    front = [s['sampled_melted_front'] for s in (a,b)]
    fd = 0. if front == [None,None] else None if None in front else front[0]-front[1]
    test('melted_front', fd, 1.e-4)
    test('largest_order_difference', float((ca-cb).flat[np.argmax(np.abs(ca-cb))]), 1.e-4)
    for key, tol in [('feed_current',3.e-5), ('bulk_current',3.e-5),
                     ('weak_disk_current',3.e-5), ('dense_total_subsonic_share',1.e-4)]:
        test(key, a['current'][key]-b['current'][key], tol)
    result = dict(passed=all(test['passed'] for test in tests.values()), comparisons=tests,
                  primary_sha256=a['sha256'], secondary_sha256=b['sha256'], shape=xa.tolist())
    atomic_json(study/'paired_held.json', result)
    return result['passed']


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study', type=Path)
    raise SystemExit(0 if compare(parser.parse_args().study) else 1)
